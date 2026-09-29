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
    definition_hash TEXT,
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

# Added after the schema above was first deployed (an existing database, on
# disk before this column existed, needs an ALTER TABLE -- `CREATE TABLE IF
# NOT EXISTS` alone never adds a column to an already-existing table). NULL
# means "never solved for", distinct from a real (possibly zero) flux.
_FLUX_COLUMN_MIGRATION = "ALTER TABLE route_reactions ADD COLUMN flux REAL"

# Same reasoning: `tasks.definition_hash` (see mteapy.tasks.task_definition_hash)
# is used to guard route-enumeration resume against a since-edited task
# definition; NULL means "unknown" (e.g. a task registered before this
# column existed), which callers should treat as "can't verify, don't resume".
_DEFINITION_HASH_COLUMN_MIGRATION = "ALTER TABLE tasks ADD COLUMN definition_hash TEXT"


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
    existing_columns = {row[1] for row in conn.execute("PRAGMA table_info(route_reactions)")}
    if "flux" not in existing_columns:
        conn.execute(_FLUX_COLUMN_MIGRATION)
        conn.commit()
    existing_task_columns = {row[1] for row in conn.execute("PRAGMA table_info(tasks)")}
    if "definition_hash" not in existing_task_columns:
        conn.execute(_DEFINITION_HASH_COLUMN_MIGRATION)
        conn.commit()
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


def register_task(
    conn: sqlite3.Connection, source: str, task_id: str, description: str = "",
    definition_hash: str | None = None,
) -> None:
    """Insert or update a task's metadata.

    `definition_hash` (see `mteapy.tasks.task_definition_hash`) is optional:
    passing it always updates the stored value, but omitting it (leaving
    the default `None`) never clobbers a hash a previous call already set --
    a caller that doesn't know/care about resume-safety can still call this
    freely without erasing another caller's hash.
    """
    conn.execute(
        "INSERT INTO tasks (source, task_id, description, definition_hash) VALUES (?, ?, ?, ?) "
        "ON CONFLICT (source, task_id) DO UPDATE SET "
        "description = excluded.description, "
        "definition_hash = COALESCE(excluded.definition_hash, tasks.definition_hash)",
        (source, task_id, description, definition_hash),
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
    fluxes: list[dict[str, float]] | None = None,
) -> list[int]:
    """Persist one task's enumeration output: new distinct routes (deduped by
    reaction-set hash against whatever this (source, task_id, model) already
    has) plus an upserted `enumeration_runs` summary row.

    `routes` should already come from `mteapy.enumeration`, in the order
    found (index 0 is the pFBA reference route). Returns the route_id
    assigned to each entry of `routes`, in the same order -- for a route
    that already existed (matching reaction-set hash), that existing
    route_id, not a new one.

    `fluxes`, if given, must be the same length as `routes`: each entry is
    that route's `RouteResult.fluxes` (reaction_id -> solved flux), stored
    alongside the route so a later caller never has to re-solve the LP just
    to draw the route's network diagram. Omitting it leaves `flux` NULL for
    any newly-inserted route (unchanged behavior for callers that only care
    about topology); a route that already existed gets its flux backfilled
    from `fluxes` even though its reaction rows aren't otherwise touched.
    """
    if fluxes is not None and len(fluxes) != len(routes):
        raise ValueError(f"fluxes must have the same length as routes ({len(routes)}), got {len(fluxes)}")

    route_ids = []
    for i, reactions in enumerate(routes):
        h = reaction_set_hash(reactions)
        existing = conn.execute(
            "SELECT route_id FROM routes WHERE source = ? AND task_id = ? AND model_id = ? AND reaction_set_hash = ?",
            (source, task_id, model_id, h),
        ).fetchone()
        route_flux = fluxes[i] if fluxes is not None else None

        if existing is not None:
            route_id = existing[0]
            if route_flux:
                # Only touch reactions actually present in route_flux -- a
                # caller-supplied partial dict must never null out a
                # previously-stored value for a reaction it simply didn't
                # include, matching save_route_fluxes's semantics.
                conn.executemany(
                    "UPDATE route_reactions SET flux = ? WHERE route_id = ? AND reaction_id = ?",
                    [(flux, route_id, rid) for rid, flux in route_flux.items()],
                )
        else:
            cur = conn.execute(
                "INSERT INTO routes (source, task_id, model_id, reaction_set_hash, n_reactions, first_seen_at) "
                "VALUES (?, ?, ?, ?, ?, ?)",
                (source, task_id, model_id, h, len(reactions), _now()),
            )
            route_id = cur.lastrowid
            if route_flux is not None:
                conn.executemany(
                    "INSERT INTO route_reactions (route_id, reaction_id, flux) VALUES (?, ?, ?)",
                    [(route_id, rid, route_flux.get(rid)) for rid in reactions],
                )
            else:
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
        # len(set(route_ids)), not len(routes): if a caller ever passed two
        # entries that hash to the same route (both resolve to the same
        # route_id), n_routes must reflect the true distinct count, not the
        # raw input length.
        (source, task_id, model_id, status, len(set(route_ids)), max_routes, int(hit_cap), int(truncated), elapsed_seconds, _now()),
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


def load_route_fluxes(conn: sqlite3.Connection, route_id: int) -> dict[str, float]:
    """reaction_id -> flux for one route, restricted to reactions that
    actually have a persisted value (NULL rows -- never solved for -- are
    excluded, not returned as 0.0).

    A caller wanting a route's full flux vector should check the result
    against the route's known reaction set: fewer entries than that means
    "not yet computed for this route", not "these reactions carry zero
    flux" -- fall back to `mteapy.network.compute_route_fluxes` and then
    `save_route_fluxes` to persist it for next time.
    """
    rows = conn.execute(
        "SELECT reaction_id, flux FROM route_reactions WHERE route_id = ? AND flux IS NOT NULL",
        (route_id,),
    ).fetchall()
    return dict(rows)


def save_route_fluxes(conn: sqlite3.Connection, route_id: int, fluxes: dict[str, float]) -> None:
    """Persist a solved flux value for each of a route's reactions.

    Meant to be called the first (and only) time a route's flux is needed
    -- e.g. lazily, from a visualizer, right after
    `mteapy.network.compute_route_fluxes` -- so that LP solve never has to
    happen again for this route.
    """
    conn.executemany(
        "UPDATE route_reactions SET flux = ? WHERE route_id = ? AND reaction_id = ?",
        [(flux, route_id, rid) for rid, flux in fluxes.items()],
    )
    conn.commit()


def latest_model_id(conn: sqlite3.Connection) -> int | None:
    row = conn.execute("SELECT model_id FROM models ORDER BY model_id DESC LIMIT 1").fetchone()
    return row[0] if row else None


def list_tasks(conn: sqlite3.Connection, model_id: int) -> list[dict]:
    """Every (source, task_id) with at least one enumerated route for `model_id`.

    Returns dicts with `source`, `task_id`, `description` (from the `tasks`
    table -- blank if that task was never registered with one) and
    `n_routes` (this model's route count) -- the listing a task-picker UI
    needs, without loading every route's own reaction set.

    Ordered numerically for purely-numeric task ids (every task list this
    project actually uses), with any non-numeric id (task lists are not
    required to use numeric ids -- see `mteapy.tasks`'s module docstring)
    sorted after those, alphabetically -- a plain `CAST(... AS INTEGER)`
    would silently collapse every non-numeric id to 0, clustering them all
    together in an arbitrary order instead.
    """
    rows = conn.execute(
        """
        SELECT r.source, r.task_id, COALESCE(t.description, ''), COUNT(*)
        FROM routes r
        LEFT JOIN tasks t ON t.source = r.source AND t.task_id = r.task_id
        WHERE r.model_id = ?
        GROUP BY r.source, r.task_id
        ORDER BY (r.task_id GLOB '[0-9]*') DESC, CAST(r.task_id AS INTEGER), r.task_id
        """,
        (model_id,),
    ).fetchall()
    return [
        {"source": source, "task_id": task_id, "description": description, "n_routes": n_routes}
        for source, task_id, description, n_routes in rows
    ]


def get_task_definition_hash(conn: sqlite3.Connection, source: str, task_id: str) -> str | None:
    """The `mteapy.tasks.task_definition_hash` stored for (source, task_id),
    or None if the task was never registered with one (either never
    registered at all, or registered before this column existed)."""
    row = conn.execute(
        "SELECT definition_hash FROM tasks WHERE source = ? AND task_id = ?", (source, task_id),
    ).fetchone()
    return row[0] if row else None


def get_enumeration_status(conn: sqlite3.Connection, source: str, task_id: str, model_id: int) -> dict | None:
    """The latest enumeration_runs summary for (source, task_id, model_id),
    or None if it has never been enumerated against this model at all.

    Includes a computed `is_exhaustive`: True only if the last run ended
    'optimal' without being capped at `max_routes` or externally truncated
    -- i.e. the solver's next cut attempt came back infeasible, proving no
    further alternate route exists. A resume should skip any task where
    this is already True (there is nothing more to find) and otherwise
    seed from its currently-stored routes (see `load_task_routes` and
    `mteapy.enumeration.iter_alternate_routes`'s `seed_routes`).
    """
    row = conn.execute(
        "SELECT status, n_routes, max_routes, hit_cap, truncated, cumulative_time_seconds, last_run_at "
        "FROM enumeration_runs WHERE source = ? AND task_id = ? AND model_id = ?",
        (source, task_id, model_id),
    ).fetchone()
    if row is None:
        return None
    status, n_routes, max_routes, hit_cap, truncated, cumulative_time_seconds, last_run_at = row
    return {
        "status": status,
        "n_routes": n_routes,
        "max_routes": max_routes,
        "hit_cap": bool(hit_cap),
        "truncated": bool(truncated),
        "cumulative_time_seconds": cumulative_time_seconds,
        "last_run_at": last_run_at,
        "is_exhaustive": status == "optimal" and not hit_cap and not truncated,
    }


def reset_source(conn: sqlite3.Connection, source: str, model_id: int) -> None:
    """Delete all routes/route_reactions/enumeration_runs for (source, model_id).

    Leaves the `tasks` rows themselves (description, definition_hash) and
    the `models` row alone -- only this source's enumeration *results* for
    this model are discarded, so re-running from scratch doesn't need to
    re-supply task descriptions, and other sources/models sharing the same
    `tasks`/`models` rows are unaffected.
    """
    route_ids = [
        row[0] for row in conn.execute(
            "SELECT route_id FROM routes WHERE source = ? AND model_id = ?", (source, model_id),
        ).fetchall()
    ]
    conn.executemany("DELETE FROM route_reactions WHERE route_id = ?", [(rid,) for rid in route_ids])
    conn.execute("DELETE FROM routes WHERE source = ? AND model_id = ?", (source, model_id))
    conn.execute("DELETE FROM enumeration_runs WHERE source = ? AND model_id = ?", (source, model_id))
    conn.commit()
