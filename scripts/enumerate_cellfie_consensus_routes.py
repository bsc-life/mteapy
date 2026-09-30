"""Enumerate alternate routes for the CellFie-consensus task list
(metabolicTasks_CellfieConsensus.txt -- the task list used in the AGS-TIDE
paper, github.com/bsc-life/ags-paper) against our current Human-GEM v2.0.1
model, storing results into the same routes_human2.db under
source="cellfie_consensus" (or --source), alongside the existing
source="full" (metabolicTasks_Full.txt) routes already there.

Deliberately NOT re-deriving the paper's own older Human-GEM build (the one
bundled in mteapy's own package data, ~13085 reactions) -- using our
current model version here means downstream scores won't exactly reproduce
the paper's original published TIDE numbers, which is expected and
accepted: the goal is a usable, topology-aware route database for this
task list on the model version we're already using everywhere else in
this project, not a historical reproduction. Each route's flux vector is
captured and persisted alongside it (mteapy.enumeration.RouteResult /
mteapy.routes.record_enumeration_result's `fluxes` argument), so a later
network view never needs to re-solve for it.

--mode selects how to handle a (source, model) that already has results:
  (unset)  "plain" -- always fully re-enumerate every task from scratch,
           relying on record_enumeration_result's dedup-by-hash to avoid
           duplicate rows. This is what a first-ever run does anyway.
  reset    Wipe all existing routes/enumeration_runs for this
           (source, model_id) first (mteapy.routes.reset_source), then
           enumerate everything from scratch. Task metadata (description,
           definition_hash) is kept.
  resume   Skip any task already exhaustively enumerated (see
           mteapy.routes.get_enumeration_status); for a task that
           previously hit its max_routes cap, seed the MILP with its
           already-known routes (mteapy.enumeration's `seed_routes`) and
           search only for genuinely new ones beyond those. Refuses (loudly,
           per-task) to resume a task whose definition has changed since
           the seed routes were recorded (mteapy.tasks.task_definition_hash)
           -- an edited task list is a real possibility in this project
           (see the Human-GEM curation-fix work), and silently seeding cuts
           from a stale definition would produce wrong results with no error.

Supports --limit N and --start-at TASK_ID too, for timing/partial runs
independent of --mode.
"""
import argparse
import os
import subprocess
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

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))       # mteapy/scripts/
MTEAPY_DATA = os.path.join(os.path.dirname(SCRIPT_DIR), "src", "mteapy", "data")
# Human-GEM is a sibling checkout of the whole workspace, not something
# bundled into mteapy -- see PROVENANCE.md for why the model/routes DB are
# a frozen copy inside mteapy while the task-list source stays external.
WORKSPACE = os.path.dirname(os.path.dirname(SCRIPT_DIR))

MODEL_PATH = os.path.join(MTEAPY_DATA, "HumanGEM_v201.xml")
DB_PATH = os.path.join(MTEAPY_DATA, "routes_human2.db")
TASK_FILE = os.path.join(WORKSPACE, "Human-GEM", "data", "metabolicTasks", "metabolicTasks_CellfieConsensus.txt")
SOURCE = "cellfie_consensus"
MAX_ROUTES = 10


def _git_provenance(file_path: str) -> tuple[str | None, str | None]:
    """Best-effort (origin_repo, origin_ref) for the git repo containing
    `file_path`, or (None, None) if it's not inside one (or git isn't
    available). Used to record which exact commit of the task list's
    source repo (e.g. Human-GEM) a route database's tasks came from."""
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=None, help="Only process the first N tasks (for timing).")
    parser.add_argument("--start-at", type=int, default=0, help="Skip tasks with id < this.")
    parser.add_argument("--solver", type=str, default=None, help="COBRApy solver name (e.g. 'gurobi', 'glpk'). Default: whatever cobra picks.")
    parser.add_argument("--source", type=str, default=SOURCE, help="source label to store routes under (lets a second run, e.g. a different solver, be stored separately for comparison).")
    parser.add_argument("--mode", choices=["reset", "resume"], default=None,
                        help="'reset': wipe existing routes for this (source, model) first. "
                             "'resume': skip already-exhaustive tasks, continue capped ones from their stored routes. "
                             "Default: plain full re-run of every task.")
    parser.add_argument("--task-file", type=str, default=TASK_FILE, help="Path to a RAVEN-style task list file.")
    parser.add_argument("--only", type=str, default=None,
                        help="Comma-separated task ids to process, ignoring --limit/--start-at (for targeted repairs).")
    args = parser.parse_args()
    source = args.source
    task_file = args.task_file
    only_ids = set(args.only.split(",")) if args.only else None

    print("Loading model...", flush=True)
    model = read_sbml_model(MODEL_PATH)
    if args.solver:
        model.solver = args.solver
    print(f"solver: {model.solver.interface.__name__}", flush=True)
    lookup = build_metabolite_lookup(model)

    conn = connect(DB_PATH)
    model_id = register_model(
        conn, MODEL_PATH, sha256=model_file_sha256(MODEL_PATH),
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
                model, task, max_routes=MAX_ROUTES, met_lookup=lookup, seed_routes=seed_reaction_sets,
                stop_info=stop_info,
            )
        except Exception as exc:  # noqa: BLE001
            elapsed = time.time() - start
            print(f"[{i}/{len(tasks)}] task {task.id} ERROR ({elapsed:.1f}s): {exc}", flush=True)
            if seed_reaction_sets:
                # Preserve the already-good seed data and its accurate count
                # rather than clobbering it with an empty-routes/error record.
                record_enumeration_result(
                    conn, source, task.id, model_id, seed_reaction_sets, status="error",
                    max_routes=MAX_ROUTES, hit_cap=False, truncated=True, elapsed_seconds=elapsed,
                    fluxes=seed_fluxes,
                )
            else:
                record_enumeration_result(
                    conn, source, task.id, model_id, [], status="error",
                    max_routes=MAX_ROUTES, hit_cap=False, truncated=False, elapsed_seconds=elapsed,
                )
            continue

        elapsed = time.time() - start
        all_reaction_sets = seed_reaction_sets + [r.reactions for r in new_results]
        all_fluxes = seed_fluxes + [r.fluxes for r in new_results]
        reason = stop_info.get("reason")

        if not all_reaction_sets:
            # No seeds and no routes found at all -- the task itself is infeasible.
            print(f"[{i}/{len(tasks)}] task {task.id} infeasible ({elapsed:.1f}s)", flush=True)
            record_enumeration_result(
                conn, source, task.id, model_id, [], status="infeasible",
                max_routes=MAX_ROUTES, hit_cap=False, truncated=False, elapsed_seconds=elapsed,
            )
            continue

        # hit_cap=True only when *this* search genuinely ran out of budget
        # (there may be more); truncated=True when it stopped early for the
        # numerical degeneracy reason (mteapy.enumeration's stop_info) --
        # neither "more definitely exist" nor "proven exhaustive", so a
        # future --mode resume must not treat this task as settled either way.
        hit_cap = reason == "max_routes_reached"
        truncated = reason == "degenerate_duplicate"
        record_enumeration_result(
            conn, source, task.id, model_id, all_reaction_sets, status="optimal",
            max_routes=MAX_ROUTES, hit_cap=hit_cap, truncated=truncated, elapsed_seconds=elapsed,
            fluxes=all_fluxes,
        )
        degenerate_note = " [stopped: numerical degeneracy]" if truncated else ""
        sizes = [len(r) for r in all_reaction_sets]
        new_note = f", {len(new_results)} new" if seed_reaction_sets else ""
        print(f"[{i}/{len(tasks)}] task {task.id} ({task.description}): "
              f"{len(all_reaction_sets)} route(s){new_note}, sizes {sizes} ({elapsed:.1f}s){degenerate_note}", flush=True)

    total = time.time() - t0
    print(f"Done: {len(tasks)} tasks in {total:.1f}s ({total / max(len(tasks), 1):.1f}s/task avg).", flush=True)


if __name__ == "__main__":
    main()
