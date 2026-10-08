# Architecture

## What this library actually does

Standard metabolic task scoring (CellFie, TIDE) does two things that
throw away real biology: it picks **one fixed reaction set** per task
(usually via one parsimonious-flux solve), and it collapses a reaction's
GPR rule to **one flattened number** (AND = min, OR = max over gene
values). mteapy keeps both of those as the default ("classic") behavior
for compatibility, but adds an orthogonal axis that scores what classic
mode discards:

- **Topology**: a task can admit several, equally optimal *routes*
  (distinct reaction sets) -- `mteapy.enumeration` finds all of them via
  the recursive MILP integer-cut method of Lee, Phalakornkule, Domach &
  Grossmann (2000): solve for the minimum-total-flux reference route,
  then repeatedly force at least one previously-used reaction off and
  re-solve, until infeasible or a cap is hit.
- **Regulation**: within one *fixed* reaction, an OR in its GPR is really
  offering a choice of which candidate enzyme complex realizes it.
  `mteapy.complexes.get_enzymes` decomposes a GPR into its candidate
  complexes; `mteapy.context_scoring.assess_reaction` scores each
  complex (min of its member genes -- the AND requirement) and picks the
  best-supported one per sample.

`mapping_strategy` ("classic" | "context-aware") is the CLI-facing name
for this axis, and it's orthogonal to *which scoring convention* is used
(CellFie / TIDE / TAS) -- see `_add_mapping_strategy_args` in `parser.py`,
shared by all three.

## Package layout

```
src/mteapy/
  registry.py        Model registry: finds model folders by their
                      manifest.json, integrity-checks them (file hashes,
                      database built against this model), resolves a key or
                      path to a `ModelEntry`. Nothing model-specific is
                      hard-coded elsewhere.
  methods.py         Declarative scoring-method specs (parameters, defaults,
                      choices, help) -- the CLI and the GUI both render from
                      it (TAS now).
  taskdb.py          Self-contained task + route database (schema v2): model
                      entities, task lists/tasks, routes + fluxes. See
                      "Task + route database" below.
  enumeration.py      Alternate-route enumeration (the recursive MILP).
  complexes.py        GPR -> candidate-complex decomposition.
  context_scoring.py  The context-aware scoring engine: assess_reaction,
                      score_task, score_tasks_matrix, compute_TAS.
                      TaskActivityReport is the core return type -- TAS
                      is named after it on purpose.
  cellfie.py          CellFie: calculate_GAL/RAL, calculate_CellFie_scores
                      (classic) and _context_aware (via score_tasks_matrix),
                      compute_CellFie (the dispatcher both CLI paths use).
  tide.py             TIDE / TIDE-essential, same classic/context-aware split.
  task_model.py       Builds a per-task constrained model copy from a
                      MetabolicTask (IN/OUT/EQU -> pseudo reactions).
  tasks.py            RAVEN-style task-list file parsing.
  network.py          Route -> bipartite graph for the webapp's viewer
                      (BiGG/EC annotation, currency-metabolite handling,
                      flux-direction edges).
  utils.py            map_gpr (the plain, non-context-aware GPR projection
                      TAS's reaction-level math reduces to), check_ensemblid,
                      add_task_metadata, parallel-worker helpers.
  parser.py           All CLI argument definitions (argparse).
  cmds/
    run_mtea.py        `run-mtea analyze {TIDE,CellFie,TAS}` dispatch.
    enumerate_routes.py `run-mtea tasks enumerate-routes` -- builds/extends
                        the task/route database (`tasks import`, `tasks enumerate-routes`).
  data/models/HumanGEM/  The bundled model folder: manifest.json,
                      HumanGEM_v201.xml, routes.db (LFS), tasks/*.txt,
                      annotations/*.tsv.
  data/               Also HumanGEM.xml.gz (older, see version-mismatch note
                      below), task_structure_matrix.tsv, task_metadata.tsv,
                      essential-genes tables -- the classic/TIDE path's inputs.
webapp/
  server.py           FastAPI app: task listing, dataset upload, live
                      scoring, route-network JSON for the frontend.
  static/              Vanilla JS/D3/dagre frontend (no build step).
tests/                 pytest; toy_model fixture pattern for fast,
                      hand-constructed-network tests alongside a few
                      tests that load the real bundled model.
```

## Task + route database (`mteapy.taskdb`, schema v2)

One self-contained SQLite database per model (`meta.schema_version = 2`;
`taskdb.connect` refuses anything else and points at `cmds/migrate_db.py`).
The model XML is *not* embedded -- `models.sha256` pins the exact file
(`taskdb.verify_model_file`), and every consumer checks it.

- `models`, `reactions(rxn_pk, rxn_id)`, `metabolites(met_pk, met_id, name,
  compartment)` -- the model's entities, stored once. Everything below
  references them by integer key, which both shrinks the route tables and
  lets the importer reject a task that names a metabolite/reaction the
  model doesn't have.
- `task_lists(name, model_id, origin_file/sha256/repo/ref, ...)` -- a named
  task list against a model (`HumanGEM-Full`, `CellFie`). A list's name is
  the single handle callers select by; there is no separate "source" vs
  "task_list" notion any more (that split caused the double-counting bug
  in `decisions.md`).
- `tasks(task_pk, task_list_id, task_id, description, system, subsystem,
  should_fail, comments, definition_hash, valid, import_note)` plus child
  tables `task_inputs`, `task_outputs`, `task_equations` +
  `task_equation_terms`, `task_changed_bounds`. A task's full definition
  lives here; the original file is only an import source
  (`run-mtea tasks import`). `taskdb.load_task` rebuilds a
  `MetabolicTask` whose `task_definition_hash` equals the stored one
  (tested). **Task lists are immutable** once imported: to change a task,
  import the edited file under a new name.
- A task that fails validation (unresolvable metabolite/reaction,
  duplicate IN/OUT, or a task-file COMMENTS cell starting `INVALID:`) is
  stored with `valid = 0` and an `import_note`, without definition rows.
  Consumers skip it; the DB stays a faithful copy of the list.
- `enumeration_runs(task_pk, status, n_routes, max_routes, hit_cap,
  truncated, cumulative_time_seconds, solver, solver_version)` -- latest
  attempt per task; solver provenance lives here, per task.
  `is_exhaustive` (computed) is `status=='optimal' and not hit_cap and
  not truncated`.
- `routes(route_id, task_pk, reaction_set_hash, n_reactions)` /
  `route_reactions(route_id, rxn_pk, flux NOT NULL)` (`WITHOUT ROWID`) --
  deduplicated by reaction-set hash, every route with a solved signed flux
  for every reaction, so no viewer ever solves an LP.
  `record_enumeration_result` *requires* the fluxes
  (`mteapy.enumeration.RouteResult.fluxes`) and refuses a partial dict.

Currently two task lists against Human-GEM v2.0.1: `HumanGEM-Full` (257
tasks) and `CellFie` (193). Importing results produced elsewhere (e.g.
MN5 greasy JSON, which carries supports only) goes through
`cmds/import_routes.py`, which solves the fluxes on a route-restricted
submodel (`network.compute_route_fluxes_submodel`) first.

## A known, accepted model-version mismatch

The **classic** mapping-strategy CLI path loads `data/HumanGEM.xml.gz`
(an older bundled model) and `data/task_structure_matrix.tsv`. The
**context-aware** path loads whatever model the model folder's `routes.db` was
actually built against (checked by sha256: `models/HumanGEM/HumanGEM_v201.xml`, newer). These are
*not* the same model version. This is intentional and already documented
in `enumerate_routes.py`'s own docstring -- don't "fix" it by silently
swapping one file for the other; if it ever needs reconciling, that's a
deliberate decision, not a bug fix.

## Webapp (`webapp/server.py`)

Local-only FastAPI app, not part of the installable package. Key things
that are *not* obvious from a first read:

- Every lookup is keyed by **`(task_list, task_id)`**, never `task_id`
  alone -- different task lists reuse the same numeric task_id for
  unrelated tasks. Task definitions come from the database
  (`taskdb.load_task`), not from task files. A task flagged invalid gets
  a 409 from the network endpoint.
- `/api/tasks/{task_list}/{task_id}/network` makes `dataset_id`/`sample`
  **optional**. With neither given, every route scores against an empty
  signal, which makes every route tie (score 0.0 everywhere) --
  `score_task` returns *all* tied routes as "winning", which is exactly
  the plain "show me the topology, no data" mode. Capped at
  `MAX_TOPOLOGY_PANELS` (15): a task can have ~100 routes and building
  every panel eagerly in one request would hold the single-process
  server's event loop for the duration. Fluxes are plain DB lookups
  (nothing is solved at request time).
- Models come from `mteapy.registry` (a model folder's `manifest.json`,
  integrity-checked; unusable models are *listed* with their `problems`) and
  load lazily into a per-model `ModelContext` (model, DB connection,
  caches). Annotation tables come from the manifest, not a sibling
  checkout. `MTEAPY_MODELS` overrides where models are looked for.
- Scoring is a **Run** (`webapp/runs.py`): `POST /api/runs {model,
  task_list, dataset_id, method, params}` validates against
  `mteapy.methods`, then scores every task of the list against every
  sample in a worker thread (`context_scoring.score_tasks_report`);
  the page polls `GET /api/runs/{id}` for progress and fetches the
  tasks x samples matrix from `/results`. The worker uses its own DB
  connection (the handlers' shared one is not thread-safe) and a
  per-(model, task list) cache of routes + complex cache. The network
  endpoint takes `run_id` + `sample` and scores with *that run's*
  parameters so it always matches the table. Results live in memory
  (last 20 runs); saving/loading them is not built yet.
- Method parameters are never hard-coded in JS: `GET /api/methods` serves the
  spec and the Analysis bar renders its form from it.
- The frontend can hide the data/analysis/legend/task/detail panels
  (state in `localStorage`).
- The frontend paints that no-data mode with a neutral fill, not the
  normal evidence colors -- "no evidence" (red) would otherwise read as
  a negative finding instead of "nothing asked yet".
