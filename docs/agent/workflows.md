# Workflows

Concrete steps for tasks that come up repeatedly. All paths below are
relative to this repo's root unless stated otherwise.

## Dev environment setup (from scratch)

This repo has its own independent venv at `venv/` (not shared with any
other project). If it's missing or broken (see `debugging.md`'s venv
entry), rebuild it:

```sh
# Match whatever Python version your solver installs were actually built
# for -- 3.10 last time, because of a prebuilt CPLEX egg.
python3.10 -m venv venv
source venv/bin/activate
python -m pip install --upgrade pip

# requirements.txt already pins cobra-netgraph via git+ssh (not a local
# path) -- this alone is a real end-to-end install, pulling that
# dependency fresh from its own repo.
python -m pip install -r requirements.txt

# This package itself, editable.
python -m pip install -e . --no-deps

# Solver + optional extras.
python -m pip install gurobipy python-libsbml

# CPLEX bindings, if MN5-parity local testing is needed (not required
# for day-to-day work, which uses Gurobi locally per the solver-choice rule):
cd /opt/ibm/ILOG/CPLEX_Studio221/cplex/python/3.10/x86-64_linux  # adjust version
python setup.py install
cd -

# Optional, for interactive/notebook-based work against this venv:
python -m pip install jupyterlab ipykernel
python -m ipykernel install --user --name mteapy --display-name "mteapy (venv)"
```

Verify: `python -c "import cobra, mteapy, mteapy.routes, cobra_netgraph,
gurobipy; print('ok')"`.

## Running tests

```sh
source venv/bin/activate
pytest -q
```

152 tests expected to pass as of the current feature branch state.

## Running the CLI

```sh
# TAS, context-aware, min aggregation (default), full task list:
run-mtea analyze TAS expr.tsv --gene_col geneID \
    --mapping-strategy context-aware --routes-source full

# TAS, classic (single fixed reaction set per task):
run-mtea analyze TAS expr.tsv --gene_col geneID --mapping-strategy classic

# CellFie, context-aware:
run-mtea analyze CellFie expr.tsv --gene_col geneID \
    --mapping-strategy context-aware --routes-source full

# TIDE (always needs a differential-expression file, not raw expression):
run-mtea analyze TIDE dea.tsv --lfc_col log2FoldChange
```

`--routes-source` must be an actual `source` registered in the routes
database -- use `--routes-source full` or `--routes-source
cellfie_consensus_gurobi` for the two currently-registered task-list
families (see `architecture.md`).

## Enumerating routes for a task list

```sh
run-mtea tasks enumerate-routes path/to/task_list.txt \
    --source <label> --task-list {full,cellfie} --solver gurobi \
    --max-routes 100 --mode {reset,resume}
```

`--mode resume` skips tasks already proven exhaustive and seeds new
search from existing routes for the rest; refuses to resume a task whose
definition changed since those routes were found. `--mode reset` wipes
existing routes for that `(source, model)` first -- it prompts for
confirmation unless `--yes` is passed, since it's destructive.

## Running the webapp

```sh
cd webapp
source ../venv/bin/activate
python -m uvicorn server:app --port 8765
```

Then open `http://127.0.0.1:8765/` (or `ssh -L 8765:localhost:8765
<host>` from a remote machine, then open it locally). The server loads
the model and routes database once at startup, which takes a few
seconds -- wait for `Ready: model_id=...` in its log before expecting
`/api/tasks` to respond.

## Adding a new `analyze` method (the TAS precedent)

1. Add a `compute_<Name>` function, following `compute_TAS`'s signature
   shape (`expr_data`, `task_structure`, `model`, your method's own
   parameters, `mapping_strategy`, `tasks_routes`) -- reuse
   `score_tasks_matrix` for the context-aware path rather than
   reimplementing route/complex scoring.
2. Add a `<Name>_parser` in `parser.py`, mirroring `CellFie_parser`'s
   structure; call `_add_mapping_strategy_args(<Name>_parser)` for the
   shared classic/context-aware flags.
3. Add an `elif args.analyze_command == "<Name>":` branch in
   `cmds/run_mtea.py`'s `main()`, mirroring the CellFie branch (file
   checks, gene-id validation, loading the right model/task_structure
   per mapping strategy, calling your `compute_` function, saving
   results via `add_task_metadata`).
4. Update the `else:` branch's usage message and `README.md`'s framework
   table + CLI usage line to list the new method.
5. Add tests in `tests/test_context_scoring.py` (or wherever the
   underlying function lives) covering both mapping strategies and any
   new error paths, using the `toy_model` fixture for fast,
   hand-constructed-network cases.

## Git LFS

`src/mteapy/data/routes_human2.db` is LFS-tracked. On a fresh clone:
`git lfs install` (once per machine) before cloning, or `git lfs pull`
after an already-done clone that only got the pointer stub.
