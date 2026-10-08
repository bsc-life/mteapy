"""Build a schema-v2 database (`mteapy.taskdb`) from a legacy routes database.

    python -m mteapy.cmds.migrate_db OLD.db MODEL.xml NEW.db \\
        --model-name Human-GEM --model-version 2.0.1 \\
        --task-list Full:full:TASKS_Full.txt:cplex[:LEGACY_FILE] \\
        --task-list CellFie:cellfie:TASKS_Cellfie.txt:gurobi

Each ``--task-list`` is ``NEW_NAME:OLD_SOURCE:TASK_FILE:SOLVER`` with an
optional fifth field, LEGACY_FILE: the task file the legacy routes were
actually enumerated against (default: TASK_FILE). Its sha256 must equal the
one the legacy database recorded for that source.

Tasks are re-imported from TASK_FILE (validated against the model). A legacy
route is carried over only if the definition it was enumerated for (the
hash the legacy database stored, or -- for tasks registered before hashes
were stored -- the hash from LEGACY_FILE) equals the task's current
definition hash, the task is valid, and the route is non-empty. Everything
dropped is counted by reason.
"""

from __future__ import annotations

import argparse
import sqlite3
import sys
from collections import Counter

from cobra.io import read_sbml_model

from mteapy import taskdb
from mteapy.taskdb import model_file_sha256, task_source_sha256
from mteapy.tasks import parse_task_file, task_definition_hash


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("old_db")
    ap.add_argument("model")
    ap.add_argument("new_db")
    ap.add_argument("--model-name", required=True)
    ap.add_argument("--model-version")
    ap.add_argument("--model-origin-repo")
    ap.add_argument("--model-origin-ref")
    ap.add_argument("--task-list", action="append", required=True, metavar="NAME:SOURCE:FILE:SOLVER")
    ap.add_argument("--task-origin-repo")
    ap.add_argument("--task-origin-ref")
    args = ap.parse_args(argv)

    old = sqlite3.connect(f"file:{args.old_db}?mode=ro", uri=True)
    model = read_sbml_model(args.model)
    old_model_hash = old.execute("SELECT sha256 FROM models").fetchone()[0]
    if model_file_sha256(args.model) != old_model_hash:
        print("model file hash does not match the legacy database's model", file=sys.stderr)
        return 2

    new = taskdb.create_db(args.new_db)
    model_id = taskdb.register_model(
        new, model, args.model_name, model_file_sha256(args.model), args.model_version,
        args.model_origin_repo, args.model_origin_ref,
    )
    rxn_pks = taskdb.reaction_pks(new, model_id)
    old_model_id = old.execute("SELECT model_id FROM models").fetchone()[0]

    for spec in args.task_list:
        name, source, task_file, solver, *rest = spec.split(":")
        legacy_file = rest[0] if rest else task_file
        recorded_sha = old.execute("SELECT sha256 FROM task_sources WHERE source = ?", (source,)).fetchone()[0]
        if task_source_sha256(legacy_file) != recorded_sha:
            print(f"[{name}] legacy task file {legacy_file} does not match the sha256 the legacy "
                  f"database recorded for source {source!r}", file=sys.stderr)
            return 2
        legacy_hash = {t.id: task_definition_hash(t) for t in parse_task_file(legacy_file)}
        tasks = taskdb.import_task_list(
            new, model, model_id, name, task_file, origin_sha256=task_source_sha256(task_file),
            origin_repo=args.task_origin_repo, origin_ref=args.task_origin_ref,
        )
        valid = {tid for tid, pk in tasks.items()
                 if new.execute("SELECT valid FROM tasks WHERE task_pk = ?", (pk,)).fetchone()[0]}
        print(f"[{name}] {len(tasks)} tasks imported, {len(tasks) - len(valid)} flagged invalid")

        old_hash = {
            tid: stored or legacy_hash.get(tid)
            for tid, stored in old.execute("SELECT task_id, definition_hash FROM tasks WHERE source = ?", (source,))
        }
        dropped, kept_tasks = Counter(), set()
        for route_id, task_id, rhash, first_seen in old.execute(
            "SELECT route_id, task_id, reaction_set_hash, first_seen_at FROM routes WHERE source = ? AND model_id = ?",
            (source, old_model_id),
        ).fetchall():
            if task_id not in tasks:
                dropped["task not in new list"] += 1
                continue
            task_pk = tasks[task_id]
            new_hash = new.execute("SELECT definition_hash FROM tasks WHERE task_pk = ?", (task_pk,)).fetchone()[0]
            if old_hash.get(task_id) != new_hash:
                dropped["task definition changed"] += 1
                continue
            if task_id not in valid:
                dropped["task invalid"] += 1
                continue
            rows = old.execute(
                "SELECT reaction_id, flux FROM route_reactions WHERE route_id = ?", (route_id,)
            ).fetchall()
            if not rows:
                dropped["empty route"] += 1
                continue
            if any(flux is None for _, flux in rows):
                dropped["route without fluxes"] += 1
                continue
            new_route = new.execute(
                "INSERT INTO routes (task_pk, reaction_set_hash, n_reactions, first_seen_at) VALUES (?, ?, ?, ?)",
                (task_pk, rhash, len(rows), first_seen),
            ).lastrowid
            new.executemany(
                "INSERT INTO route_reactions VALUES (?, ?, ?)", [(new_route, rxn_pks[r], f) for r, f in rows]
            )
            kept_tasks.add(task_id)

        for (task_id, status, n_routes, max_routes, hit_cap, truncated, secs, last_run) in old.execute(
            "SELECT task_id, status, n_routes, max_routes, hit_cap, truncated, cumulative_time_seconds, last_run_at "
            "FROM enumeration_runs WHERE source = ? AND model_id = ?", (source, old_model_id),
        ).fetchall():
            if task_id not in kept_tasks:
                continue
            new.execute(
                "INSERT INTO enumeration_runs VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (tasks[task_id], status, n_routes, max_routes, hit_cap, truncated, secs, solver, None, last_run),
            )
        new.commit()
        n_routes = new.execute(
            "SELECT COUNT(*) FROM routes JOIN tasks USING (task_pk) JOIN task_lists USING (task_list_id) "
            "WHERE task_lists.name = ?", (name,)).fetchone()[0]
        print(f"[{name}] routes carried over: {n_routes}; tasks with routes: {len(kept_tasks)}; "
              f"dropped: {dict(dropped)}")

    new.execute("VACUUM")
    return 0


if __name__ == "__main__":
    sys.exit(main())
