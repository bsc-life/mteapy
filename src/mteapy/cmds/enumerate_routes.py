"""`run-mtea tasks import` / `run-mtea tasks enumerate-routes` -- build and
extend the self-contained task + route database (`mteapy.taskdb`) that
context-aware scoring reads from (see also `mteapy.enumeration`,
`mteapy.context_scoring`).

No expression/sample data is involved at all. `import` stores a RAVEN-style
task list, validated against the model, as a named task list; 
`enumerate-routes` then finds, for each valid task of that list, every
alternate reaction set (up to --max-routes) that can accomplish it, together
with the solved flux of each route, and persists them. This is
dataset-agnostic structural precomputation, run once per (model, task list)
and reused by every sample later scored against it -- which is why these are
`tasks` commands, not `analyze` ones.

A task list's definitions are immutable once imported. To change a task,
import the edited file under a new list name (or reset and re-import), so
stored routes can never silently belong to an older definition.

--mode selects how to handle a task list that already has results:
  (unset)  "plain" -- always fully re-enumerate every task from scratch,
           relying on record_enumeration_result's dedup-by-hash to avoid
           duplicate rows. This is what a first-ever run does anyway.
  reset    Wipe all existing routes/enumeration summaries for this task
           list first (taskdb.reset_task_list_routes), then enumerate
           everything from scratch. The task definitions are kept. Prompts
           for confirmation first (this permanently discards enumeration
           results), unless --yes.
  resume   Skip any task already exhaustively enumerated (see
           taskdb.get_enumeration_status); for a task that previously hit its
           max_routes cap, seed the MILP with its already-known routes
           (mteapy.enumeration's `seed_routes`) and search only for
           genuinely new ones beyond those.
"""
from __future__ import annotations

import os
import sqlite3
import subprocess
import sys
import time

from cobra.io import read_sbml_model

from mteapy import registry, taskdb
from mteapy.enumeration import enumerate_alternate_routes
from mteapy.task_model import build_metabolite_lookup


def resolve_paths(model: str | None, db: str | None, model_file: str | None) -> tuple[str, str]:
    """(db_path, model_path) for a command: `--model` (a registered model key
    or a model folder; default: the bundled Human-GEM) supplies both, and an
    explicit `--db` / `--model-file` overrides its half. The registry is only
    consulted when something is missing, so building a brand-new database
    from explicit files needs no registered model. Exits with a clear
    message instead of a traceback when the model cannot be resolved."""
    if db and model_file:
        return db, model_file
    try:
        entry = registry.resolve_model(model)
    except (KeyError, ValueError) as exc:
        print(f"ERROR: {exc.args[0]}", file=sys.stderr)
        sys.exit(2)
    return db or entry.db_path, model_file or entry.model_path


def _git_provenance(file_path: str) -> tuple[str | None, str | None]:
    """Best-effort (origin_repo, origin_ref) for the git repo containing
    `file_path`, or (None, None) if it's not inside one (or git isn't
    available). Used to record which exact commit of the task list's
    source repo a route database's tasks came from."""
    directory = os.path.dirname(os.path.abspath(file_path))
    try:
        remote = subprocess.run(
            ["git", "-C", directory, "remote", "get-url", "origin"],
            capture_output=True, text=True, check=True, timeout=5,
        ).stdout.strip()
        ref = subprocess.run(
            ["git", "-C", directory, "rev-parse", "HEAD"],
            capture_output=True, text=True, check=True, timeout=5,
        ).stdout.strip()
        return remote, ref
    except Exception:  # noqa: BLE001 -- provenance is best-effort, never fatal
        return None, None


def _confirm_reset(conn, task_list: str, skip_confirmation: bool) -> bool:
    """Prints how many routes --mode reset would permanently delete and asks
    for confirmation, unless --yes was passed. Returns True to proceed."""
    n_routes = conn.execute(
        "SELECT COUNT(*) FROM routes r JOIN tasks t USING (task_pk) JOIN task_lists tl USING (task_list_id) "
        "WHERE tl.name = ?", (task_list,),
    ).fetchone()[0]
    if n_routes == 0:
        return True
    print(f"WARNING: --mode reset will permanently delete {n_routes} existing route(s) of task list {task_list!r}.")
    if skip_confirmation:
        return True
    if not sys.stdin.isatty():
        print("Refusing to reset non-interactively without --yes.", flush=True)
        return False
    reply = input("Type 'yes' to continue: ").strip().lower()
    return reply == "yes"


def _solver_name(model) -> str:
    return model.solver.interface.__name__.rsplit(".", 1)[-1].removesuffix("_interface")


def run_import(args) -> None:
    """`run-mtea tasks import`: store a task file as a named task list."""
    db_path, model_path = resolve_paths(args.model, args.db, args.model_file)

    print("Loading model...", flush=True)
    model = read_sbml_model(model_path)
    sha = taskdb.model_file_sha256(model_path)

    if os.path.exists(db_path):
        conn = taskdb.connect(db_path)
        known = taskdb.get_model(conn, sha)
    else:
        conn, known = taskdb.create_db(db_path), None
        print(f"created {db_path}")
    if known is None:
        if not args.model_name:
            print("This model is not yet registered in the database: pass --model-name (and optionally "
                  "--model-version).", file=sys.stderr)
            sys.exit(2)
        model_id = taskdb.register_model(conn, model, args.model_name, sha, args.model_version)
    else:
        model_id = known["model_id"]

    origin_repo, origin_ref = _git_provenance(args.task_file)
    try:
        tasks = taskdb.import_task_list(
            conn, model, model_id, args.task_list, args.task_file, description=args.description,
            origin_sha256=taskdb.task_source_sha256(args.task_file), origin_repo=origin_repo, origin_ref=origin_ref,
        )
    except sqlite3.IntegrityError:
        print(f"Task list {args.task_list!r} already exists in {db_path}; task lists are immutable -- import "
              f"under a new name.", file=sys.stderr)
        sys.exit(2)
    invalid = conn.execute(
        "SELECT t.task_id, t.import_note FROM tasks t JOIN task_lists tl USING (task_list_id) "
        "WHERE tl.name = ? AND t.valid = 0", (args.task_list,)).fetchall()
    print(f"Imported {len(tasks)} tasks as {args.task_list!r}; {len(invalid)} flagged invalid.")
    for task_id, note in invalid:
        print(f"  task {task_id}: {note}")


def run(args) -> None:
    """`run-mtea tasks enumerate-routes`."""
    task_list = args.task_list
    db_path, model_path = resolve_paths(args.model, args.db, args.model_file)
    only_ids = set(args.only.split(",")) if args.only else None

    try:
        conn = taskdb.connect(db_path)
        tl = taskdb.get_task_list(conn, task_list)
        taskdb.verify_model_file(conn, tl["model_id"], model_path)
    except (FileNotFoundError, KeyError, ValueError) as exc:
        print(f"ERROR: {exc.args[0]}", file=sys.stderr)
        sys.exit(2)

    print("Loading model...", flush=True)
    model = read_sbml_model(model_path)
    if args.solver:
        model.solver = args.solver
    solver = _solver_name(model)
    print(f"solver: {solver}", flush=True)
    lookup = build_metabolite_lookup(model)

    if args.mode == "reset":
        if not _confirm_reset(conn, task_list, args.skip_confirmation):
            print("Aborted.", flush=True)
            return
        print(f"--mode reset: clearing existing routes of task list {task_list!r}...", flush=True)
        taskdb.reset_task_list_routes(conn, task_list)

    rows = [r for r in taskdb.list_tasks(conn, task_list) if r["valid"]]
    skipped_invalid = [r["task_id"] for r in taskdb.list_tasks(conn, task_list) if not r["valid"]]
    if skipped_invalid:
        print(f"skipping {len(skipped_invalid)} task(s) flagged invalid: {skipped_invalid}", flush=True)
    if only_ids is not None:
        rows = [r for r in rows if r["task_id"] in only_ids]
    else:
        rows = [r for r in rows if not r["task_id"].isdigit() or int(r["task_id"]) >= args.start_at]
        if args.limit:
            rows = rows[: args.limit]
    print(f"{len(rows)} tasks to process from task list {task_list!r} (mode={args.mode or 'plain'})", flush=True)

    t0 = time.time()
    for i, row in enumerate(rows, 1):
        task_pk = taskdb.get_task_pk(conn, task_list, row["task_id"])
        task = taskdb.load_task(conn, task_pk)
        seed_reaction_sets: list[frozenset[str]] = []
        seed_fluxes: list[dict[str, float]] = []

        if args.mode == "resume":
            status = taskdb.get_enumeration_status(conn, task_pk)
            if status and status["is_exhaustive"]:
                print(f"[{i}/{len(rows)}] task {task.id}: already exhaustive ({status['n_routes']} routes), skipping", flush=True)
                continue
            if status:
                existing = taskdb.load_task_routes(conn, task_pk)
                seed_reaction_sets = list(existing.values())
                seed_fluxes = [taskdb.load_route_fluxes(conn, rid) for rid in existing.keys()]

        start = time.time()
        stop_info: dict = {}
        try:
            new_results = enumerate_alternate_routes(
                model, task, max_routes=args.max_routes, met_lookup=lookup, seed_routes=seed_reaction_sets,
                stop_info=stop_info,
            )
        except Exception as exc:  # noqa: BLE001
            elapsed = time.time() - start
            print(f"[{i}/{len(rows)}] task {task.id} ERROR ({elapsed:.1f}s): {exc}", flush=True)
            taskdb.record_enumeration_result(
                conn, task_pk, seed_reaction_sets, seed_fluxes, status="error", max_routes=args.max_routes,
                hit_cap=False, truncated=bool(seed_reaction_sets), elapsed_seconds=elapsed, solver=solver,
            )
            continue

        elapsed = time.time() - start
        all_reaction_sets = seed_reaction_sets + [r.reactions for r in new_results]
        all_fluxes = seed_fluxes + [r.fluxes for r in new_results]
        reason = stop_info.get("reason")

        if not all_reaction_sets:
            print(f"[{i}/{len(rows)}] task {task.id} infeasible ({elapsed:.1f}s)", flush=True)
            taskdb.record_enumeration_result(
                conn, task_pk, [], [], status="infeasible", max_routes=args.max_routes, hit_cap=False,
                truncated=False, elapsed_seconds=elapsed, solver=solver,
            )
            continue

        hit_cap = reason == "max_routes_reached"
        truncated = reason == "degenerate_duplicate"
        taskdb.record_enumeration_result(
            conn, task_pk, all_reaction_sets, all_fluxes, status="optimal", max_routes=args.max_routes,
            hit_cap=hit_cap, truncated=truncated, elapsed_seconds=elapsed, solver=solver,
        )
        degenerate_note = " [stopped: numerical degeneracy]" if truncated else ""
        sizes = [len(r) for r in all_reaction_sets]
        new_note = f", {len(new_results)} new" if seed_reaction_sets else ""
        print(f"[{i}/{len(rows)}] task {task.id} ({task.description}): "
              f"{len(all_reaction_sets)} route(s){new_note}, sizes {sizes} ({elapsed:.1f}s){degenerate_note}", flush=True)

    total = time.time() - t0
    print(f"Done: {len(rows)} tasks in {total:.1f}s ({total / max(len(rows), 1):.1f}s/task avg).", flush=True)
