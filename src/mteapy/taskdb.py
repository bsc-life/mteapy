"""Self-contained task + route database (schema v2).

One database per model. It holds the model's entities (reactions,
metabolites), the task lists defined against that model, every task's full
definition (so the original task file is only an import source), and the
enumerated alternate routes with their solved fluxes. Metabolites and
reactions are stored once and referenced by integer key everywhere else,
which both shrinks the route tables and lets the importer verify that every
task references only things that exist in the model.

The model XML itself is not embedded: `models.sha256` pins which file the
database was built against.

Task definitions round-trip: `load_task` rebuilds a `MetabolicTask` whose
`task_definition_hash` equals the one stored at import time.

A task that cannot be imported cleanly (an unresolvable metabolite or
reaction, a repeated IN/OUT metabolite, or a task file COMMENTS cell
starting with ``INVALID:``) is still stored, with ``valid = 0`` and an
``import_note``, but without definition rows -- the database stays a
faithful copy of the list while consumers can skip what cannot be scored.
"""

from __future__ import annotations

import hashlib
import json
import sqlite3
import warnings
from datetime import datetime, timezone

from cobra.core import Model

from mteapy.task_model import (
    _EQU_ARROW_PATTERN, _EQU_TERM_SPLIT_PATTERN, _parse_equation_term, build_metabolite_lookup,
)
from mteapy.tasks import (
    BoundedMetabolite, ChangedBound, EquationConstraint, MetabolicTask, parse_task_file, task_definition_hash,
)

SCHEMA_VERSION = 2
INVALID_MARKER = "INVALID:"

SCHEMA = """
CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);

CREATE TABLE models (
    model_id INTEGER PRIMARY KEY,
    name TEXT NOT NULL,
    version TEXT,
    sha256 TEXT NOT NULL UNIQUE,
    n_reactions INTEGER,
    n_genes INTEGER,
    origin_repo TEXT,
    origin_ref TEXT,
    created_at TEXT NOT NULL
);
CREATE TABLE reactions (
    rxn_pk INTEGER PRIMARY KEY,
    model_id INTEGER NOT NULL REFERENCES models(model_id),
    rxn_id TEXT NOT NULL,
    UNIQUE (model_id, rxn_id)
);
CREATE TABLE metabolites (
    met_pk INTEGER PRIMARY KEY,
    model_id INTEGER NOT NULL REFERENCES models(model_id),
    met_id TEXT NOT NULL,
    name TEXT NOT NULL,
    compartment TEXT NOT NULL,
    UNIQUE (model_id, met_id)
);

CREATE TABLE task_lists (
    task_list_id INTEGER PRIMARY KEY,
    name TEXT NOT NULL UNIQUE,
    model_id INTEGER NOT NULL REFERENCES models(model_id),
    description TEXT,
    origin_file TEXT,
    origin_sha256 TEXT,
    origin_repo TEXT,
    origin_ref TEXT,
    imported_at TEXT NOT NULL
);
CREATE TABLE tasks (
    task_pk INTEGER PRIMARY KEY,
    task_list_id INTEGER NOT NULL REFERENCES task_lists(task_list_id),
    task_id TEXT NOT NULL,
    description TEXT,
    system TEXT,
    subsystem TEXT,
    should_fail INTEGER NOT NULL DEFAULT 0,
    print_flux INTEGER NOT NULL DEFAULT 0,
    comments TEXT,
    annotations TEXT,
    definition_hash TEXT NOT NULL,
    valid INTEGER NOT NULL DEFAULT 1,
    import_note TEXT,
    UNIQUE (task_list_id, task_id)
);
CREATE TABLE task_inputs (
    task_pk INTEGER NOT NULL REFERENCES tasks(task_pk),
    ord INTEGER NOT NULL,
    met_pk INTEGER NOT NULL REFERENCES metabolites(met_pk),
    lower_bound REAL NOT NULL,
    upper_bound REAL NOT NULL,
    PRIMARY KEY (task_pk, ord)
);
CREATE TABLE task_outputs (
    task_pk INTEGER NOT NULL REFERENCES tasks(task_pk),
    ord INTEGER NOT NULL,
    met_pk INTEGER NOT NULL REFERENCES metabolites(met_pk),
    lower_bound REAL NOT NULL,
    upper_bound REAL NOT NULL,
    PRIMARY KEY (task_pk, ord)
);
CREATE TABLE task_equations (
    equation_pk INTEGER PRIMARY KEY,
    task_pk INTEGER NOT NULL REFERENCES tasks(task_pk),
    ord INTEGER NOT NULL,
    equation TEXT NOT NULL,
    lower_bound REAL NOT NULL,
    upper_bound REAL NOT NULL
);
CREATE TABLE task_equation_terms (
    equation_pk INTEGER NOT NULL REFERENCES task_equations(equation_pk),
    met_pk INTEGER NOT NULL REFERENCES metabolites(met_pk),
    coefficient REAL NOT NULL
);
CREATE TABLE task_changed_bounds (
    task_pk INTEGER NOT NULL REFERENCES tasks(task_pk),
    ord INTEGER NOT NULL,
    rxn_pk INTEGER NOT NULL REFERENCES reactions(rxn_pk),
    lower_bound REAL NOT NULL,
    upper_bound REAL NOT NULL,
    PRIMARY KEY (task_pk, ord)
);

CREATE TABLE enumeration_runs (
    task_pk INTEGER PRIMARY KEY REFERENCES tasks(task_pk),
    status TEXT NOT NULL,
    n_routes INTEGER NOT NULL,
    max_routes INTEGER NOT NULL,
    hit_cap INTEGER NOT NULL,
    truncated INTEGER NOT NULL,
    cumulative_time_seconds REAL NOT NULL,
    solver TEXT,
    solver_version TEXT,
    last_run_at TEXT NOT NULL
);
CREATE TABLE routes (
    route_id INTEGER PRIMARY KEY,
    task_pk INTEGER NOT NULL REFERENCES tasks(task_pk),
    reaction_set_hash TEXT NOT NULL,
    n_reactions INTEGER NOT NULL,
    first_seen_at TEXT NOT NULL,
    UNIQUE (task_pk, reaction_set_hash)
);
CREATE TABLE route_reactions (
    route_id INTEGER NOT NULL REFERENCES routes(route_id),
    rxn_pk INTEGER NOT NULL REFERENCES reactions(rxn_pk),
    flux REAL NOT NULL,
    PRIMARY KEY (route_id, rxn_pk)
) WITHOUT ROWID;
"""


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def reaction_set_hash(reactions) -> str:
    return hashlib.sha256("\n".join(sorted(reactions)).encode("utf-8")).hexdigest()


def create_db(path: str) -> sqlite3.Connection:
    """Create a new, empty schema-v2 database at `path` (must not exist)."""
    import os

    if os.path.exists(path):
        raise FileExistsError(path)
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA foreign_keys = ON")
    conn.executescript(SCHEMA)
    conn.execute("INSERT INTO meta VALUES ('schema_version', ?)", (str(SCHEMA_VERSION),))
    conn.commit()
    return conn


def register_model(
    conn: sqlite3.Connection, model: Model, name: str, sha256: str, version: str | None = None,
    origin_repo: str | None = None, origin_ref: str | None = None,
) -> int:
    """Insert the model row plus every reaction and metabolite; returns model_id."""
    cur = conn.execute(
        "INSERT INTO models (name, version, sha256, n_reactions, n_genes, origin_repo, origin_ref, created_at) "
        "VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
        (name, version, sha256, len(model.reactions), len(model.genes), origin_repo, origin_ref, _now()),
    )
    model_id = cur.lastrowid
    conn.executemany(
        "INSERT INTO reactions (model_id, rxn_id) VALUES (?, ?)", [(model_id, r.id) for r in model.reactions]
    )
    conn.executemany(
        "INSERT INTO metabolites (model_id, met_id, name, compartment) VALUES (?, ?, ?, ?)",
        [(model_id, m.id, m.name, m.compartment) for m in model.metabolites],
    )
    conn.commit()
    return model_id


def reaction_pks(conn: sqlite3.Connection, model_id: int) -> dict[str, int]:
    return dict(conn.execute("SELECT rxn_id, rxn_pk FROM reactions WHERE model_id = ?", (model_id,)))


def _metabolite_pks(conn: sqlite3.Connection, model_id: int) -> dict[str, int]:
    return dict(conn.execute("SELECT met_id, met_pk FROM metabolites WHERE model_id = ?", (model_id,)))


def _equation_terms(equation: str) -> list[tuple[str, float]]:
    """(``name[compartment]`` token, signed coefficient) per term, mirroring
    `mteapy.task_model._add_equation_reaction` (reactants negative)."""
    arrow = _EQU_ARROW_PATTERN.search(equation)
    if not arrow:
        raise ValueError(f"Equation {equation!r} does not contain '=>' or '<=>'")
    terms: dict[str, float] = {}
    for side, sign in ((equation[: arrow.start()], -1.0), (equation[arrow.end():], 1.0)):
        for term in _EQU_TERM_SPLIT_PATTERN.split(side):
            if term.strip():
                coef, token = _parse_equation_term(term.strip())
                terms[token] = terms.get(token, 0.0) + sign * coef
    return list(terms.items())


def import_task_list(
    conn: sqlite3.Connection, model: Model, model_id: int, name: str, task_file: str,
    description: str | None = None, origin_sha256: str | None = None, origin_repo: str | None = None,
    origin_ref: str | None = None,
) -> dict[str, int]:
    """Parse `task_file` and store the task list `name` against `model_id`.

    Returns ``{task_id: task_pk}``. Tasks that fail validation are stored
    with ``valid = 0`` (see the module docstring), never dropped.
    """
    lookup = build_metabolite_lookup(model)
    met_pks = _metabolite_pks(conn, model_id)
    rxn_pks = reaction_pks(conn, model_id)

    cur = conn.execute(
        "INSERT INTO task_lists (name, model_id, description, origin_file, origin_sha256, origin_repo, origin_ref, "
        "imported_at) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
        (name, model_id, description, task_file, origin_sha256, origin_repo, origin_ref, _now()),
    )
    task_list_id = cur.lastrowid

    def met_pk(token: str) -> int:
        met_id = lookup.get(token.strip())
        if met_id is None:
            raise KeyError(f"metabolite {token!r} not in model")
        return met_pks[met_id]

    result: dict[str, int] = {}
    for task in parse_task_file(task_file):
        problems: list[str] = []
        if task.comments.strip().startswith(INVALID_MARKER):
            problems.append(task.comments.strip())

        rows_in, rows_out, rows_eq, rows_ch = [], [], [], []
        for bucket, rows in ((task.inputs, rows_in), (task.outputs, rows_out)):
            seen = set()
            for k, bm in enumerate(bucket):
                try:
                    pk = met_pk(bm.metabolite)
                except KeyError as exc:
                    problems.append(str(exc))
                    continue
                if pk in seen:
                    problems.append(f"duplicate metabolite {bm.metabolite!r}")
                seen.add(pk)
                rows.append((k, pk, bm.lower_bound, bm.upper_bound))
        for k, eq in enumerate(task.equations):
            try:
                terms = [(met_pk(tok), coef) for tok, coef in _equation_terms(eq.equation)]
            except (KeyError, ValueError) as exc:
                problems.append(str(exc))
                continue
            rows_eq.append((k, eq, terms))
        for k, cb in enumerate(task.changed_bounds):
            if cb.reaction_id not in rxn_pks:
                problems.append(f"reaction {cb.reaction_id!r} not in model")
                continue
            rows_ch.append((k, rxn_pks[cb.reaction_id], cb.lower_bound, cb.upper_bound))

        valid = not problems
        cur = conn.execute(
            "INSERT INTO tasks (task_list_id, task_id, description, system, subsystem, should_fail, print_flux, "
            "comments, annotations, definition_hash, valid, import_note) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (task_list_id, task.id, task.description, task.system, task.subsystem, int(task.should_fail),
             int(task.print_flux), task.comments, json.dumps(task.annotations, sort_keys=True),
             task_definition_hash(task), int(valid), "; ".join(problems) or None),
        )
        task_pk = cur.lastrowid
        result[task.id] = task_pk
        if not valid:
            continue
        conn.executemany("INSERT INTO task_inputs VALUES (?, ?, ?, ?, ?)", [(task_pk, *r) for r in rows_in])
        conn.executemany("INSERT INTO task_outputs VALUES (?, ?, ?, ?, ?)", [(task_pk, *r) for r in rows_out])
        for k, eq, terms in rows_eq:
            eq_pk = conn.execute(
                "INSERT INTO task_equations (task_pk, ord, equation, lower_bound, upper_bound) VALUES (?, ?, ?, ?, ?)",
                (task_pk, k, eq.equation, eq.lower_bound, eq.upper_bound),
            ).lastrowid
            conn.executemany(
                "INSERT INTO task_equation_terms VALUES (?, ?, ?)", [(eq_pk, pk, coef) for pk, coef in terms]
            )
        conn.executemany(
            "INSERT INTO task_changed_bounds VALUES (?, ?, ?, ?, ?)", [(task_pk, *r) for r in rows_ch]
        )
    conn.commit()
    return result


def load_task(conn: sqlite3.Connection, task_pk: int) -> MetabolicTask:
    """Rebuild the `MetabolicTask` stored under `task_pk`.

    Raises `ValueError` for a task flagged ``valid = 0`` (it has no
    definition rows to rebuild from).
    """
    row = conn.execute(
        "SELECT task_id, description, system, subsystem, should_fail, print_flux, comments, annotations, valid, "
        "import_note FROM tasks WHERE task_pk = ?", (task_pk,)
    ).fetchone()
    if row is None:
        raise KeyError(task_pk)
    task_id, description, system, subsystem, should_fail, print_flux, comments, annotations, valid, note = row
    if not valid:
        raise ValueError(f"task {task_id!r} is flagged invalid: {note}")

    def bounded(table: str) -> list[BoundedMetabolite]:
        return [
            BoundedMetabolite(f"{name}[{comp}]", lb, ub)
            for name, comp, lb, ub in conn.execute(
                f"SELECT m.name, m.compartment, t.lower_bound, t.upper_bound FROM {table} t "
                "JOIN metabolites m USING (met_pk) WHERE t.task_pk = ? ORDER BY t.ord", (task_pk,)
            )
        ]

    return MetabolicTask(
        id=task_id, description=description or "", should_fail=bool(should_fail), system=system or "",
        subsystem=subsystem or "", inputs=bounded("task_inputs"), outputs=bounded("task_outputs"),
        equations=[
            EquationConstraint(eq, lb, ub) for eq, lb, ub in conn.execute(
                "SELECT equation, lower_bound, upper_bound FROM task_equations WHERE task_pk = ? ORDER BY ord",
                (task_pk,))
        ],
        changed_bounds=[
            ChangedBound(rxn_id, lb, ub) for rxn_id, lb, ub in conn.execute(
                "SELECT r.rxn_id, c.lower_bound, c.upper_bound FROM task_changed_bounds c "
                "JOIN reactions r USING (rxn_pk) WHERE c.task_pk = ? ORDER BY c.ord", (task_pk,))
        ],
        print_flux=bool(print_flux), comments=comments or "", annotations=json.loads(annotations or "{}"),
    )


# ---------------------------------------------------------------------------
# Opening, lookup and querying
# ---------------------------------------------------------------------------

def model_file_sha256(path: str) -> str:
    return _file_sha256(path)


def task_source_sha256(path: str) -> str:
    return _file_sha256(path)


def _file_sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def connect(path: str, check_same_thread: bool = True) -> sqlite3.Connection:
    """Open an existing schema-v2 database (foreign keys enforced).

    Refuses a missing file (to create one use
    `create_db`) and any database that is not schema v2, with a message
    pointing at `mteapy.cmds.migrate_db` for legacy ones.

    `check_same_thread=False` lets a server share one connection across its
    worker threads; the caller then owns serializing access.
    """
    import os

    if not os.path.exists(path):
        raise FileNotFoundError(f"routes database not found: {path}")
    conn = sqlite3.connect(path, check_same_thread=check_same_thread)
    conn.execute("PRAGMA foreign_keys = ON")
    try:
        version = conn.execute("SELECT value FROM meta WHERE key = 'schema_version'").fetchone()
    except sqlite3.OperationalError:
        version = None
    if version is None or int(version[0]) != SCHEMA_VERSION:
        conn.close()
        raise ValueError(
            f"{path} is not a schema-v{SCHEMA_VERSION} routes database "
            f"(legacy databases can be converted with `python -m mteapy.cmds.migrate_db`)"
        )
    return conn


def get_model(conn: sqlite3.Connection, sha256: str) -> dict | None:
    row = conn.execute(
        "SELECT model_id, name, version, sha256, n_reactions, n_genes FROM models WHERE sha256 = ?", (sha256,)
    ).fetchone()
    return dict(zip(("model_id", "name", "version", "sha256", "n_reactions", "n_genes"), row)) if row else None


def verify_model_file(conn: sqlite3.Connection, model_id: int, path: str) -> None:
    """Raise `ValueError` unless the file at `path` is the model `model_id` was built from."""
    expected = conn.execute("SELECT sha256 FROM models WHERE model_id = ?", (model_id,)).fetchone()[0]
    if model_file_sha256(path) != expected:
        raise ValueError(
            f"{path} does not match the model this database was built against (sha256 {expected[:12]}...); "
            "its reaction ids may not line up with the stored routes"
        )


def task_list_routes_fingerprint(conn: sqlite3.Connection, task_list: str) -> str:
    """sha256 over every (task id, route reaction-set hash) of `task_list`.

    Saved results record it so that loading them can tell whether the routes
    they were scored against are still exactly the ones in this database."""
    task_list_id = get_task_list(conn, task_list)["task_list_id"]
    rows = conn.execute(
        "SELECT t.task_id, r.reaction_set_hash FROM routes r JOIN tasks t USING (task_pk) "
        "WHERE t.task_list_id = ?", (task_list_id,)
    ).fetchall()
    h = hashlib.sha256()
    for task_id, route_hash in sorted(rows):
        h.update(f"{task_id}\t{route_hash}\n".encode())
    return h.hexdigest()


def list_task_lists(conn: sqlite3.Connection) -> list[dict]:
    rows = conn.execute(
        "SELECT tl.name, tl.model_id, m.name, tl.description, "
        "(SELECT COUNT(*) FROM tasks t WHERE t.task_list_id = tl.task_list_id), "
        "(SELECT COUNT(DISTINCT r.task_pk) FROM routes r JOIN tasks t USING (task_pk) "
        " WHERE t.task_list_id = tl.task_list_id) "
        "FROM task_lists tl JOIN models m USING (model_id) ORDER BY tl.name"
    ).fetchall()
    return [dict(zip(("name", "model_id", "model_name", "description", "n_tasks", "n_tasks_with_routes"), r))
            for r in rows]


def get_task_list(conn: sqlite3.Connection, name: str) -> dict:
    """The task list `name`; raises `KeyError` listing the available names if it doesn't exist."""
    row = conn.execute(
        "SELECT task_list_id, name, model_id, description, origin_file, origin_sha256, origin_repo, origin_ref "
        "FROM task_lists WHERE name = ?", (name,)
    ).fetchone()
    if row is None:
        available = [r[0] for r in conn.execute("SELECT name FROM task_lists ORDER BY name")]
        raise KeyError(f"unknown task list {name!r}; available: {available}")
    return dict(zip(("task_list_id", "name", "model_id", "description", "origin_file", "origin_sha256",
                     "origin_repo", "origin_ref"), row))


def get_task_pk(conn: sqlite3.Connection, task_list: str, task_id: str) -> int | None:
    row = conn.execute(
        "SELECT t.task_pk FROM tasks t JOIN task_lists tl USING (task_list_id) WHERE tl.name = ? AND t.task_id = ?",
        (task_list, task_id),
    ).fetchone()
    return row[0] if row else None


def load_task_by_id(conn: sqlite3.Connection, task_list: str, task_id: str) -> MetabolicTask:
    pk = get_task_pk(conn, task_list, task_id)
    if pk is None:
        raise KeyError(f"task {task_id!r} not in task list {task_list!r}")
    return load_task(conn, pk)


def load_task_routes(conn: sqlite3.Connection, task_pk: int) -> dict[int, frozenset[str]]:
    """route_id -> reaction-id set, for one task."""
    routes: dict[int, set[str]] = {}
    for route_id, rxn_id in conn.execute(
        "SELECT r.route_id, x.rxn_id FROM routes r JOIN route_reactions rr USING (route_id) "
        "JOIN reactions x USING (rxn_pk) WHERE r.task_pk = ? ORDER BY r.route_id", (task_pk,)
    ):
        routes.setdefault(route_id, set()).add(rxn_id)
    return {rid: frozenset(rx) for rid, rx in routes.items()}


def load_task_list_routes(
    conn: sqlite3.Connection, task_list: str, min_routes: int = 1
) -> dict[str, dict[int, frozenset[str]]]:
    """task_id -> {route_id: reaction set} for every task of `task_list` with at
    least `min_routes` routes (single-route tasks included by default: there is
    simply nothing for the topology dimension to distinguish)."""
    get_task_list(conn, task_list)
    by_task: dict[str, dict[int, set[str]]] = {}
    for task_id, route_id, rxn_id in conn.execute(
        "SELECT t.task_id, r.route_id, x.rxn_id FROM tasks t JOIN task_lists tl USING (task_list_id) "
        "JOIN routes r USING (task_pk) JOIN route_reactions rr USING (route_id) JOIN reactions x USING (rxn_pk) "
        "WHERE tl.name = ?", (task_list,)
    ):
        by_task.setdefault(task_id, {}).setdefault(route_id, set()).add(rxn_id)
    return {
        task_id: {rid: frozenset(rx) for rid, rx in routes.items()}
        for task_id, routes in by_task.items() if len(routes) >= min_routes
    }


def load_route_fluxes(conn: sqlite3.Connection, route_id: int) -> dict[str, float]:
    """reaction_id -> solved signed flux for one route (always complete: flux is NOT NULL)."""
    return dict(conn.execute(
        "SELECT x.rxn_id, rr.flux FROM route_reactions rr JOIN reactions x USING (rxn_pk) WHERE rr.route_id = ?",
        (route_id,),
    ))


def list_tasks(conn: sqlite3.Connection, task_list: str | None = None) -> list[dict]:
    """Every task (of `task_list`, or all lists) with its route count and validity,
    ordered by list then numerically by id (non-numeric ids last)."""
    sql = (
        "SELECT tl.name, t.task_id, t.description, t.valid, t.import_note, "
        "(SELECT COUNT(*) FROM routes r WHERE r.task_pk = t.task_pk) "
        "FROM tasks t JOIN task_lists tl USING (task_list_id)"
    )
    params: tuple = ()
    if task_list is not None:
        get_task_list(conn, task_list)
        sql += " WHERE tl.name = ?"
        params = (task_list,)
    rows = [
        {"task_list": n, "task_id": tid, "description": d or "", "valid": bool(v), "import_note": note,
         "n_routes": nr}
        for n, tid, d, v, note, nr in conn.execute(sql, params)
    ]
    rows.sort(key=lambda r: (r["task_list"], not r["task_id"].isdigit(),
                             int(r["task_id"]) if r["task_id"].isdigit() else 0, r["task_id"]))
    return rows


# ---------------------------------------------------------------------------
# Recording enumeration results
# ---------------------------------------------------------------------------

def get_enumeration_status(conn: sqlite3.Connection, task_pk: int) -> dict | None:
    """The task's latest enumeration summary, or None if never enumerated.

    ``is_exhaustive`` is true only when the solver proved no further
    alternate route exists (status optimal, cap not hit, not truncated).
    """
    row = conn.execute(
        "SELECT status, n_routes, max_routes, hit_cap, truncated, cumulative_time_seconds, solver, solver_version "
        "FROM enumeration_runs WHERE task_pk = ?", (task_pk,)
    ).fetchone()
    if row is None:
        return None
    status, n_routes, max_routes, hit_cap, truncated, secs, solver, solver_version = row
    return {
        "status": status, "n_routes": n_routes, "max_routes": max_routes, "hit_cap": bool(hit_cap),
        "truncated": bool(truncated), "cumulative_time_seconds": secs, "solver": solver,
        "solver_version": solver_version,
        "is_exhaustive": status == "optimal" and not hit_cap and not truncated,
    }


def record_enumeration_result(
    conn: sqlite3.Connection,
    task_pk: int,
    routes: list[frozenset[str]],
    fluxes: list[dict[str, float]],
    status: str,
    max_routes: int,
    hit_cap: bool,
    truncated: bool,
    elapsed_seconds: float,
    solver: str | None = None,
    solver_version: str | None = None,
) -> list[int]:
    """Persist one task's enumeration output: new distinct routes (deduplicated
    by reaction-set hash against what the task already has) plus an upserted
    `enumeration_runs` summary row. Returns each entry's route_id, in order.

    `fluxes[i]` is route `i`'s solved flux for every reaction in its support
    (`mteapy.enumeration.RouteResult.fluxes`); it is required, because the
    schema stores a flux for every route reaction -- the viewer never has to
    solve anything. A reaction id absent from the model, or a flux dict that
    doesn't cover its route, raises `ValueError` before anything is written.
    """
    if len(fluxes) != len(routes):
        raise ValueError(f"fluxes must have the same length as routes ({len(routes)}), got {len(fluxes)}")
    model_id = conn.execute(
        "SELECT tl.model_id FROM tasks t JOIN task_lists tl USING (task_list_id) WHERE t.task_pk = ?", (task_pk,)
    ).fetchone()[0]
    rxn_pks = reaction_pks(conn, model_id)
    for reactions, flux in zip(routes, fluxes):
        unknown = [r for r in reactions if r not in rxn_pks]
        if unknown:
            raise ValueError(f"reactions not in the model: {sorted(unknown)[:5]}")
        missing = [r for r in reactions if flux.get(r) is None]
        if missing:
            raise ValueError(f"no flux given for route reactions: {sorted(missing)[:5]}")

    # A repeated support within one batch means the enumerator kept re-finding
    # the same route (typically solver noise leaving a cut unapplied), so the
    # cap was not really reached and the task is not proven exhaustive either:
    # record it as truncated so a resume revisits it.
    n_repeats = len(routes) - len({reaction_set_hash(r) for r in routes})
    if n_repeats:
        warnings.warn(
            f"task_pk {task_pk}: {n_repeats} repeated route support(s) dropped; recording the run as truncated",
            RuntimeWarning, stacklevel=2,
        )
        hit_cap, truncated = False, True

    route_ids = []
    for reactions, flux in zip(routes, fluxes):
        h = reaction_set_hash(reactions)
        row = conn.execute(
            "SELECT route_id FROM routes WHERE task_pk = ? AND reaction_set_hash = ?", (task_pk, h)
        ).fetchone()
        if row is not None:
            route_ids.append(row[0])
            continue
        route_id = conn.execute(
            "INSERT INTO routes (task_pk, reaction_set_hash, n_reactions, first_seen_at) VALUES (?, ?, ?, ?)",
            (task_pk, h, len(reactions), _now()),
        ).lastrowid
        conn.executemany(
            "INSERT INTO route_reactions (route_id, rxn_pk, flux) VALUES (?, ?, ?)",
            [(route_id, rxn_pks[r], float(flux[r])) for r in reactions],
        )
        route_ids.append(route_id)

    conn.execute(
        "INSERT INTO enumeration_runs (task_pk, status, n_routes, max_routes, hit_cap, truncated, "
        "cumulative_time_seconds, solver, solver_version, last_run_at) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?) "
        "ON CONFLICT (task_pk) DO UPDATE SET status = excluded.status, n_routes = excluded.n_routes, "
        "max_routes = excluded.max_routes, hit_cap = excluded.hit_cap, truncated = excluded.truncated, "
        "cumulative_time_seconds = enumeration_runs.cumulative_time_seconds + excluded.cumulative_time_seconds, "
        "solver = excluded.solver, solver_version = excluded.solver_version, last_run_at = excluded.last_run_at",
        (task_pk, status, len(set(route_ids)), max_routes, int(hit_cap), int(truncated), elapsed_seconds, solver,
         solver_version, _now()),
    )
    conn.commit()
    return route_ids


def reset_task_list_routes(conn: sqlite3.Connection, task_list: str) -> None:
    """Delete every route and enumeration summary of `task_list` (its tasks stay)."""
    task_list_id = get_task_list(conn, task_list)["task_list_id"]
    pks = "(SELECT task_pk FROM tasks WHERE task_list_id = ?)"
    conn.execute(f"DELETE FROM route_reactions WHERE route_id IN (SELECT route_id FROM routes WHERE task_pk IN {pks})",
                 (task_list_id,))
    conn.execute(f"DELETE FROM routes WHERE task_pk IN {pks}", (task_list_id,))
    conn.execute(f"DELETE FROM enumeration_runs WHERE task_pk IN {pks}", (task_list_id,))
    conn.commit()


def reset_task_routes(conn: sqlite3.Connection, task_pk: int) -> None:
    """Delete one task's routes and enumeration summary (the task itself stays)."""
    conn.execute("DELETE FROM route_reactions WHERE route_id IN (SELECT route_id FROM routes WHERE task_pk = ?)",
                 (task_pk,))
    conn.execute("DELETE FROM routes WHERE task_pk = ?", (task_pk,))
    conn.execute("DELETE FROM enumeration_runs WHERE task_pk = ?", (task_pk,))
    conn.commit()
