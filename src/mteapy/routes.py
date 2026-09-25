"""Persistence for enumerated alternate-route results.

This is the boundary between the two sides of context-aware task scoring:
`mteapy.enumeration` (slow, solver-dependent) writes here once per
model/task-list version; `mteapy.context_scoring` (fast, the default path
most users hit) only ever reads from here, never re-enumerates.

The schema matches, field for field, the database already used to store the
Human-GEM v2.0.1 route-enumeration results this package's design was
validated against, so an existing precomputed database can be read (and
extended with new tasks/models) without any migration.
"""

from __future__ import annotations

import hashlib
import sqlite3
from datetime import datetime, timezone

SCHEMA = """
CREATE TABLE IF NOT EXISTS models (
    model_id INTEGER PRIMARY KEY,
    path TEXT NOT NULL,
    sha256 TEXT NOT NULL UNIQUE,
    n_reactions INTEGER,
    n_genes INTEGER,
    solver TEXT,
    solver_version TEXT,
    created_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS tasks (
    source TEXT NOT NULL,
    task_id TEXT NOT NULL,
    description TEXT,
    PRIMARY KEY (source, task_id)
);
CREATE TABLE IF NOT EXISTS enumeration_runs (
    source TEXT NOT NULL,
    task_id TEXT NOT NULL,
    model_id INTEGER NOT NULL REFERENCES models(model_id),
    status TEXT NOT NULL,
    n_routes INTEGER NOT NULL,
    max_routes INTEGER NOT NULL,
    hit_cap INTEGER NOT NULL,
    truncated INTEGER NOT NULL,
    cumulative_time_seconds REAL NOT NULL,
    last_run_at TEXT NOT NULL,
    PRIMARY KEY (source, task_id, model_id),
    FOREIGN KEY (source, task_id) REFERENCES tasks(source, task_id)
);
CREATE TABLE IF NOT EXISTS routes (
    route_id INTEGER PRIMARY KEY,
    source TEXT NOT NULL,
    task_id TEXT NOT NULL,
    model_id INTEGER NOT NULL,
    reaction_set_hash TEXT NOT NULL,
    n_reactions INTEGER NOT NULL,
    first_seen_at TEXT NOT NULL,
    UNIQUE (source, task_id, model_id, reaction_set_hash)
);
CREATE TABLE IF NOT EXISTS route_reactions (
    route_id INTEGER NOT NULL REFERENCES routes(route_id),
    reaction_id TEXT NOT NULL,
    PRIMARY KEY (route_id, reaction_id)
);
CREATE INDEX IF NOT EXISTS idx_routes_task ON routes(source, task_id, model_id);
CREATE INDEX IF NOT EXISTS idx_route_reactions_reaction ON route_reactions(reaction_id);
"""


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def reaction_set_hash(reactions) -> str:
    """Content hash identifying a route by its reaction set, order-independent."""
    return hashlib.sha256("|".join(sorted(reactions)).encode()).hexdigest()


def model_file_sha256(path: str) -> str:
    """Hash a model file's raw bytes, for `register_model`'s dedup key."""
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def connect(path: str) -> sqlite3.Connection:
    """Open (creating if needed) a routes database at `path` with the schema applied."""
    conn = sqlite3.connect(path)
    conn.executescript(SCHEMA)
    return conn


def register_model(
    conn: sqlite3.Connection,
    path: str,
    sha256: str,
    n_reactions: int | None = None,
    n_genes: int | None = None,
    solver: str | None = None,
    solver_version: str | None = None,
) -> int:
    """Insert `models` row if this exact content hash hasn't been seen, and return its model_id either way."""
    row = conn.execute("SELECT model_id FROM models WHERE sha256 = ?", (sha256,)).fetchone()
    if row is not None:
        return row[0]
    cur = conn.execute(
        "INSERT INTO models (path, sha256, n_reactions, n_genes, solver, solver_version, created_at) "
        "VALUES (?, ?, ?, ?, ?, ?, ?)",
        (path, sha256, n_reactions, n_genes, solver, solver_version, _now()),
    )
    conn.commit()
    return cur.lastrowid


def register_task(conn: sqlite3.Connection, source: str, task_id: str, description: str = "") -> None:
    conn.execute(
        "INSERT INTO tasks (source, task_id, description) VALUES (?, ?, ?) "
        "ON CONFLICT (source, task_id) DO UPDATE SET description = excluded.description",
        (source, task_id, description),
    )
    conn.commit()


def record_enumeration_result(
    conn: sqlite3.Connection,
    source: str,
    task_id: str,
    model_id: int,
    routes: list[frozenset[str]],
    status: str,
    max_routes: int,
    hit_cap: bool,
    truncated: bool,
    elapsed_seconds: float,
) -> list[int]:
    """Persist one task's enumeration output: new distinct routes (deduped by
    reaction-set hash against whatever this (source, task_id, model) already
    has) plus an upserted `enumeration_runs` summary row.

    `routes` should already come from `mteapy.enumeration`, in the order
    found (index 0 is the pFBA reference route). Returns the route_id
    assigned to each entry of `routes`, in the same order -- for a route
    that already existed (matching reaction-set hash), that existing
    route_id, not a new one.
    """
    route_ids = []
    for reactions in routes:
        h = reaction_set_hash(reactions)
        existing = conn.execute(
            "SELECT route_id FROM routes WHERE source = ? AND task_id = ? AND model_id = ? AND reaction_set_hash = ?",
            (source, task_id, model_id, h),
        ).fetchone()
        if existing is not None:
            route_ids.append(existing[0])
            continue
        cur = conn.execute(
            "INSERT INTO routes (source, task_id, model_id, reaction_set_hash, n_reactions, first_seen_at) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            (source, task_id, model_id, h, len(reactions), _now()),
        )
        route_id = cur.lastrowid
        conn.executemany(
            "INSERT OR IGNORE INTO route_reactions (route_id, reaction_id) VALUES (?, ?)",
            [(route_id, rid) for rid in reactions],
        )
        route_ids.append(route_id)

    conn.execute(
        "INSERT INTO enumeration_runs "
        "(source, task_id, model_id, status, n_routes, max_routes, hit_cap, truncated, cumulative_time_seconds, last_run_at) "
        "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?) "
        "ON CONFLICT (source, task_id, model_id) DO UPDATE SET "
        "status = excluded.status, n_routes = excluded.n_routes, max_routes = excluded.max_routes, "
        "hit_cap = excluded.hit_cap, truncated = excluded.truncated, "
        "cumulative_time_seconds = enumeration_runs.cumulative_time_seconds + excluded.cumulative_time_seconds, "
        "last_run_at = excluded.last_run_at",
        (source, task_id, model_id, status, len(routes), max_routes, int(hit_cap), int(truncated), elapsed_seconds, _now()),
    )
    conn.commit()
    return route_ids


def load_task_routes(conn: sqlite3.Connection, source: str, task_id: str, model_id: int) -> dict[int, frozenset[str]]:
    """route_id -> reaction set, for one (source, task_id, model)."""
    route_rows = conn.execute(
        "SELECT route_id FROM routes WHERE source = ? AND task_id = ? AND model_id = ?",
        (source, task_id, model_id),
    ).fetchall()
    result = {}
    for (route_id,) in route_rows:
        reactions = frozenset(
            r[0] for r in conn.execute(
                "SELECT reaction_id FROM route_reactions WHERE route_id = ?", (route_id,)
            ).fetchall()
        )
        result[route_id] = reactions
    return result


def load_multiroute_tasks(
    conn: sqlite3.Connection, model_id: int, min_routes: int = 2
) -> dict[tuple[str, str], dict[int, frozenset[str]]]:
    """(source, task_id) -> {route_id: reaction set}, only for tasks with at least `min_routes` distinct routes.

    This is the entry point context-aware scoring actually uses: single-route
    tasks have no topological choice to make (nothing for the "topology"
    dimension to distinguish), so callers scoring for topological diversity
    should start here rather than `load_task_routes` task-by-task.
    """
    rows = conn.execute(
        "SELECT source, task_id, route_id FROM routes WHERE model_id = ?", (model_id,)
    ).fetchall()
    by_task: dict[tuple[str, str], list[int]] = {}
    for source, task_id, route_id in rows:
        by_task.setdefault((source, task_id), []).append(route_id)

    result = {}
    for key, route_ids in by_task.items():
        if len(route_ids) < min_routes:
            continue
        result[key] = {
            route_id: frozenset(
                r[0] for r in conn.execute(
                    "SELECT reaction_id FROM route_reactions WHERE route_id = ?", (route_id,)
                ).fetchall()
            )
            for route_id in route_ids
        }
    return result


def latest_model_id(conn: sqlite3.Connection) -> int | None:
    row = conn.execute("SELECT model_id FROM models ORDER BY model_id DESC LIMIT 1").fetchone()
    return row[0] if row else None
