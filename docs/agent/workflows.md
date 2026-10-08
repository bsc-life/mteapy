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

182 tests expected to pass as of the current feature branch state.

## Running the CLI

```sh
# TAS, context-aware, min aggregation (default), full task list:
run-mtea analyze TAS expr.tsv --gene_col geneID \
    --mapping-strategy context-aware --task-list HumanGEM-Full

# TAS, classic (single fixed reaction set per task):
run-mtea analyze TAS expr.tsv --gene_col geneID --mapping-strategy classic

# CellFie, context-aware:
run-mtea analyze CellFie expr.tsv --gene_col geneID \
    --mapping-strategy context-aware --task-list HumanGEM-Full

# TIDE (always needs a differential-expression file, not raw expression):
run-mtea analyze TIDE dea.tsv --lfc_col log2FoldChange
```

`--task-list` must be a task list stored in the database
(`HumanGEM-Full` is the default; `CellFie` is the other currently stored).
`--routes-db`/`--routes-model-file` select another database/model; the
model file's sha256 must match the one the database was built against.

## Models: listing, adding, and the manifest

```sh
run-mtea models list          # every model found, its task lists, and any problems
run-mtea analyze TAS expr.tsv ... --mapping-strategy context-aware --model Human-GEM-2.0.1
```

Models are found under `$MTEAPY_MODELS` (os.pathsep-separated) if set, else
`~/.mteapy/models` plus the bundled `src/mteapy/data/models/`. `--model` takes a
key or a model folder; `--routes-db` / `--routes-model-file` override its halves.

A model is a folder -- see `mteapy/registry.py`'s docstring for the full
`manifest.json` format:

```
MyModel/
  manifest.json      name, version, model/database files, sha256s, task lists, annotations
  model.xml          the SBML model
  routes.db          schema-v2 task/route database (built with `tasks import` + `enumerate-routes`)
  tasks/*.txt        the task files the database was imported from
  annotations/*.tsv  optional BiGG/EC tables (metabolites.tsv, reactions.tsv)
```

To add a model: build `routes.db` (next section), copy the model and task files
in, write the manifest (sha256 of the model file and each task file; the task-list
names must match the ones imported), then `run-mtea models list` -- it reports
any mismatch (modified model file, database built against another model, a
task list whose hash differs from the one the database was imported from).
New packaged files need a `[tool.setuptools.package-data]` pattern in
`pyproject.toml` (the bundled folder already has one); verify with
`pip wheel . --no-deps --no-build-isolation`.

## Building / extending the task + route database

```sh
# 1. Store a task list (validated against the model; creates the DB if needed).
run-mtea tasks import path/to/task_list.txt --task-list HumanGEM-Full \
    --model-name Human-GEM --model-version 2.0.1 --db new.db

# 2. Enumerate routes (with fluxes) for every valid task of that list.
run-mtea tasks enumerate-routes --task-list HumanGEM-Full --db new.db \
    --solver gurobi --max-routes 100 --mode {reset,resume}
```

`--mode resume` skips tasks already proven exhaustive and seeds new
search from existing routes for the rest. `--mode reset` wipes that
task list's routes first and prompts unless `--yes`. Tasks flagged
invalid at import are skipped. Task lists are immutable: to change a
task, edit the file and import it under a new list name.

Results computed elsewhere (e.g. the MN5 greasy jobs' one-JSON-per-task
output, which has route supports only):

```sh
python -m mteapy.cmds.import_routes new.db HumanGEM_v201.xml \
    --task-list HumanGEM-Full --json-dir results/<dir> \
    --solver-label cplex --processes 8
```

Converting a legacy (schema v1) database: `python -m
mteapy.cmds.migrate_db OLD.db MODEL.xml NEW.db --model-name ...
--task-list NAME:OLD_SOURCE:TASK_FILE:SOLVER[:LEGACY_TASK_FILE]` -- see
its docstring (routes are carried over only if their task definition
hasn't changed).

## Running the webapp

```sh
cd webapp
source ../venv/bin/activate
python -m uvicorn server:app --port 8765
```

`MTEAPY_MODELS=<model folder or dir of them>` makes it look for models elsewhere than the bundled
ones. Then open `http://127.0.0.1:8765/` (or `ssh -L 8765:localhost:8765
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

`src/mteapy/data/models/*/routes.db` is LFS-tracked. On a fresh clone:
`git lfs install` (once per machine) before cloning, or `git lfs pull`
after an already-done clone that only got the pointer stub.
