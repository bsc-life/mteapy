"""Import route-enumeration results produced elsewhere (e.g. the MN5 greasy
jobs' one-JSON-per-task output) into a schema-v2 database.

Those JSON files carry each task's route *supports* only. The database
stores a solved flux for every route reaction, so this solves each route's
minimum-total-flux LP on a route-restricted submodel
(`mteapy.network.compute_route_fluxes_submodel`) before recording it --
cheap (about a second per route), parallel across routes.

    python -m mteapy.cmds.import_routes DB MODEL.xml --task-list HumanGEM-Full \\
        --json-dir results/human2_greasy_out_fix1 --solver-label cplex --processes 8

Each JSON file is ``{"task_id", "status", "n_routes", "max_routes", "hit_cap",
"truncated", "time_seconds", "routes": [[reaction ids], ...]}``. A task that
already has routes in the database is skipped unless --replace is given
(which deletes only that task's routes first).
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from multiprocessing import Pool

from cobra.io import read_sbml_model

from mteapy import taskdb
from mteapy.network import compute_route_fluxes_submodel
from mteapy.task_model import build_metabolite_lookup

_state: dict = {}


def _init(model_path, db_path, task_list, solver):
    _state["model"] = read_sbml_model(model_path)
    if solver:
        _state["model"].solver = solver
    _state["lookup"] = build_metabolite_lookup(_state["model"])
    _state["conn"] = taskdb.connect(db_path)
    _state["task_list"] = task_list
    _state["tasks"] = {}


def _solve(job):
    task_id, k, reactions = job
    try:
        task = _state["tasks"].get(task_id)
        if task is None:
            task = _state["tasks"][task_id] = taskdb.load_task_by_id(_state["conn"], _state["task_list"], task_id)
        fluxes = compute_route_fluxes_submodel(_state["model"], task, set(reactions), _state["lookup"])
        if any(v != v for v in fluxes.values()):
            return task_id, k, None, "no optimal solution"
        return task_id, k, fluxes, None
    except Exception as exc:  # noqa: BLE001 -- one bad route must not abort the batch
        return task_id, k, None, f"{type(exc).__name__}: {exc}"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("db")
    ap.add_argument("model")
    ap.add_argument("--task-list", required=True)
    ap.add_argument("--json-dir", required=True)
    ap.add_argument("--solver-label", default=None, help="solver that produced the routes (recorded as provenance)")
    ap.add_argument("--solver", default=None, help="cobra solver for the flux solves (default: cobra's default)")
    ap.add_argument("--processes", type=int, default=1)
    ap.add_argument("--replace", action="store_true")
    args = ap.parse_args(argv)

    conn = taskdb.connect(args.db)
    tl = taskdb.get_task_list(conn, args.task_list)
    taskdb.verify_model_file(conn, tl["model_id"], args.model)

    payloads, jobs = {}, []
    for path in sorted(glob.glob(os.path.join(args.json_dir, "*.json"))):
        data = json.load(open(path))
        task_id = str(data["task_id"])
        pk = taskdb.get_task_pk(conn, args.task_list, task_id)
        if pk is None:
            print(f"skip {path}: task {task_id!r} not in task list {args.task_list!r}")
            continue
        if not conn.execute("SELECT valid FROM tasks WHERE task_pk = ?", (pk,)).fetchone()[0]:
            print(f"skip task {task_id}: flagged invalid")
            continue
        if taskdb.load_task_routes(conn, pk) and not args.replace:
            print(f"skip task {task_id}: already has routes (use --replace)")
            continue
        payloads[task_id] = (pk, data)
        jobs += [(task_id, k, r) for k, r in enumerate(data["routes"])]
    print(f"{len(payloads)} tasks, {len(jobs)} routes to solve", flush=True)

    flux_by = {}
    failed = []
    with Pool(args.processes, initializer=_init, initargs=(args.model, args.db, args.task_list, args.solver)) as pool:
        for i, (task_id, k, fluxes, err) in enumerate(pool.imap_unordered(_solve, jobs, chunksize=4), 1):
            if err:
                failed.append((task_id, k, err))
            else:
                flux_by[(task_id, k)] = fluxes
            if i % 200 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)} solved, {len(failed)} failed", flush=True)
    if failed:
        for task_id, k, err in failed[:20]:
            print(f"FAILED task {task_id} route {k}: {err}", file=sys.stderr)
        print("Nothing written for tasks with failed routes.", file=sys.stderr)
    bad_tasks = {t for t, _, _ in failed}

    for task_id, (pk, data) in payloads.items():
        if task_id in bad_tasks:
            continue
        if args.replace:
            taskdb.reset_task_routes(conn, pk)
        routes = [frozenset(r) for r in data["routes"]]
        taskdb.record_enumeration_result(
            conn, pk, routes, [flux_by[(task_id, k)] for k in range(len(routes))], status=data["status"],
            max_routes=data["max_routes"], hit_cap=data["hit_cap"], truncated=data["truncated"],
            elapsed_seconds=data["time_seconds"], solver=args.solver_label,
        )
        print(f"task {task_id}: {len(routes)} routes recorded")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
