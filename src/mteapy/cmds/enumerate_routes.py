"""`run-mtea tasks enumerate-routes` -- build/extend the alternate-route
database context-aware scoring reads from (see `mteapy.enumeration`,
`mteapy.routes`, `mteapy.context_scoring`).

No expression/sample data is involved at all: this enumerates, for a given
model + RAVEN-style task list, every alternate reaction set (up to
--max-routes) that can accomplish each task, and persists them. It's
dataset-agnostic structural precomputation, run once per (model,
task-list) pair and reused by every sample later scored against it --
which is why it's a `tasks` command, not an `analyze` one.

--mode selects how to handle a (source, model) that already has results:
  (unset)  "plain" -- always fully re-enumerate every task from scratch,
           relying on record_enumeration_result's dedup-by-hash to avoid
           duplicate rows. This is what a first-ever run does anyway.
  reset    Wipe all existing routes/enumeration_runs for this
           (source, model_id) first (mteapy.routes.reset_source), then
           enumerate everything from scratch. Task metadata (description,
           definition_hash) is kept. Prompts for confirmation first (this
           permanently discards enumeration results), unless --yes.
  resume   Skip any task already exhaustively enumerated (see
           mteapy.routes.get_enumeration_status); for a task that
           previously hit its max_routes cap, seed the MILP with its
           already-known routes (mteapy.enumeration's `seed_routes`) and
           search only for genuinely new ones beyond those. Refuses (loudly,
           per-task) to resume a task whose definition has changed since
           the seed routes were recorded (mteapy.tasks.task_definition_hash)
           -- silently seeding cuts from a stale definition would produce
           wrong results with no error.
"""
from __future__ import annotations

import os
import subprocess
import sys
import time

from cobra.io import read_sbml_model

from mteapy.enumeration import enumerate_alternate_routes
from mteapy.routes import (
    connect, get_enumeration_status, get_task_definition_hash, load_route_fluxes, load_task_routes,
    model_file_sha256, record_enumeration_result, register_model, register_task, register_task_source,
    reset_source, task_source_sha256,
)
from mteapy.task_model import build_metabolite_lookup
from mteapy.tasks import parse_task_file, task_definition_hash

MTEAPY_DATA = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "data")
DEFAULT_MODEL_PATH = os.path.join(MTEAPY_DATA, "HumanGEM_v201.xml")
DEFAULT_DB_PATH = os.path.join(MTEAPY_DATA, "routes_human2.db")


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


def _confirm_reset(conn, source: str, model_id: int, skip_confirmation: bool) -> bool:
    """Prints how many routes --mode reset would permanently delete and asks
    for confirmation, unless --yes was passed. Returns True to proceed."""
    n_routes = conn.execute(
        "SELECT COUNT(*) FROM routes WHERE source = ? AND model_id = ?", (source, model_id),
    ).fetchone()[0]
    if n_routes == 0:
        return True
    print(f"WARNING: --mode reset will permanently delete {n_routes} existing route(s) "
          f"for source={source!r}, model_id={model_id}.")
    if skip_confirmation:
        return True
    if not sys.stdin.isatty():
        print("Refusing to reset non-interactively without --yes.", flush=True)
        return False
    reply = input("Type 'yes' to continue: ").strip().lower()
    return reply == "yes"


def run(args) -> None:
    source = args.source
    task_file = args.task_file
    model_path = args.model_file or DEFAULT_MODEL_PATH
    db_path = args.db or DEFAULT_DB_PATH
    only_ids = set(args.only.split(",")) if args.only else None

    print("Loading model...", flush=True)
    model = read_sbml_model(model_path)
    if args.solver:
        model.solver = args.solver
    print(f"solver: {model.solver.interface.__name__}", flush=True)
    lookup = build_metabolite_lookup(model)

    conn = connect(db_path)
    model_id = register_model(
        conn, model_path, sha256=model_file_sha256(model_path),
        n_reactions=len(model.reactions), n_genes=len(model.genes),
    )
    print(f"model_id={model_id} ({len(model.reactions)} reactions, {len(model.genes)} genes)", flush=True)

    origin_repo, origin_ref = _git_provenance(task_file)
    source_changed = register_task_source(
        conn, source, task_file, sha256=task_source_sha256(task_file),
        origin_repo=origin_repo, origin_ref=origin_ref,
    )
    if source_changed:
        print(f"WARNING: the task-list file for source={source!r} has changed (different content hash) "
              f"since it was last recorded -- existing routes may no longer match the current task "
              f"definitions. Consider --mode reset if that's not intended.", flush=True)

    if args.mode == "reset":
        if not _confirm_reset(conn, source, model_id, args.skip_confirmation):
            print("Aborted.", flush=True)
            return
        print(f"--mode reset: clearing existing routes for source={source!r}, model_id={model_id}...", flush=True)
        reset_source(conn, source, model_id)

    tasks = parse_task_file(task_file)
    if only_ids is not None:
        tasks = [t for t in tasks if t.id in only_ids]
    else:
        tasks = [t for t in tasks if int(t.id) >= args.start_at]
        if args.limit:
            tasks = tasks[: args.limit]
    print(f"{len(tasks)} tasks to process from {task_file} (mode={args.mode or 'plain'})", flush=True)

    t0 = time.time()
    for i, task in enumerate(tasks, 1):
        definition_hash = task_definition_hash(task)
        seed_reaction_sets: list[frozenset[str]] = []
        seed_fluxes: list[dict[str, float]] = []

        if args.mode == "resume":
            status = get_enumeration_status(conn, source, task.id, model_id)
            if status and status["is_exhaustive"]:
                register_task(conn, source, task.id, task.description, definition_hash=definition_hash)
                print(f"[{i}/{len(tasks)}] task {task.id}: already exhaustive ({status['n_routes']} routes), skipping", flush=True)
                continue
            if status:
                stored_hash = get_task_definition_hash(conn, source, task.id)
                if stored_hash is not None and stored_hash != definition_hash:
                    print(f"[{i}/{len(tasks)}] task {task.id}: SKIPPED -- task definition changed since its "
                          f"stored routes were found (hash {stored_hash[:8]} -> {definition_hash[:8]}); "
                          f"refusing to seed stale cuts. Use --mode reset to redo this source from scratch.", flush=True)
                    continue
                existing = load_task_routes(conn, source, task.id, model_id)
                seed_reaction_sets = list(existing.values())
                seed_fluxes = [load_route_fluxes(conn, rid) for rid in existing.keys()]

        register_task(conn, source, task.id, task.description, definition_hash=definition_hash)
        start = time.time()
        stop_info: dict = {}
        try:
            new_results = enumerate_alternate_routes(
                model, task, max_routes=args.max_routes, met_lookup=lookup, seed_routes=seed_reaction_sets,
                stop_info=stop_info,
            )
        except Exception as exc:  # noqa: BLE001
            elapsed = time.time() - start
            print(f"[{i}/{len(tasks)}] task {task.id} ERROR ({elapsed:.1f}s): {exc}", flush=True)
            if seed_reaction_sets:
                record_enumeration_result(
                    conn, source, task.id, model_id, seed_reaction_sets, status="error",
                    max_routes=args.max_routes, hit_cap=False, truncated=True, elapsed_seconds=elapsed,
                    fluxes=seed_fluxes,
                )
            else:
                record_enumeration_result(
                    conn, source, task.id, model_id, [], status="error",
                    max_routes=args.max_routes, hit_cap=False, truncated=False, elapsed_seconds=elapsed,
                )
            continue

        elapsed = time.time() - start
        all_reaction_sets = seed_reaction_sets + [r.reactions for r in new_results]
        all_fluxes = seed_fluxes + [r.fluxes for r in new_results]
        reason = stop_info.get("reason")

        if not all_reaction_sets:
            print(f"[{i}/{len(tasks)}] task {task.id} infeasible ({elapsed:.1f}s)", flush=True)
            record_enumeration_result(
                conn, source, task.id, model_id, [], status="infeasible",
                max_routes=args.max_routes, hit_cap=False, truncated=False, elapsed_seconds=elapsed,
            )
            continue

        hit_cap = reason == "max_routes_reached"
        truncated = reason == "degenerate_duplicate"
        record_enumeration_result(
            conn, source, task.id, model_id, all_reaction_sets, status="optimal",
            max_routes=args.max_routes, hit_cap=hit_cap, truncated=truncated, elapsed_seconds=elapsed,
            fluxes=all_fluxes,
        )
        degenerate_note = " [stopped: numerical degeneracy]" if truncated else ""
        sizes = [len(r) for r in all_reaction_sets]
        new_note = f", {len(new_results)} new" if seed_reaction_sets else ""
        print(f"[{i}/{len(tasks)}] task {task.id} ({task.description}): "
              f"{len(all_reaction_sets)} route(s){new_note}, sizes {sizes} ({elapsed:.1f}s){degenerate_note}", flush=True)

    total = time.time() - t0
    print(f"Done: {len(tasks)} tasks in {total:.1f}s ({total / max(len(tasks), 1):.1f}s/task avg).", flush=True)
