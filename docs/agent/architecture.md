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
  routes.py          SQLite routes database: models, tasks, task_sources,
                      enumeration_runs, routes, route_reactions. See
                      "Routes database schema" below.
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
                        the routes database.
  data/               Bundled HumanGEM_v201.xml, routes_human2.db (LFS),
                      HumanGEM.xml.gz (older, see version-mismatch note
                      below), task_structure_matrix.tsv, task_metadata.tsv.
webapp/
  server.py           FastAPI app: task listing, dataset upload, live
                      scoring, route-network JSON for the frontend.
  static/              Vanilla JS/D3/dagre frontend (no build step).
tests/                 pytest; toy_model fixture pattern for fast,
                      hand-constructed-network tests alongside a few
                      tests that load the real bundled model.
```

## Routes database schema (`mteapy.routes`)

- `models(model_id, path, sha256, n_reactions, n_genes, ...)` -- one row
  per distinct model file content (keyed by hash, not path).
- `tasks(source, task_id, description, definition_hash)` -- one row per
  (source, task_id); `definition_hash` lets `--mode resume` detect a
  task definition that changed since routes were last found for it.
- `task_sources(source, task_list, file_path, sha256, origin_repo,
  origin_ref, recorded_at)` -- provenance per `source`, plus `task_list`
  grouping multiple sources under one logical family (e.g. a
  solver-comparison re-run) so a caller can select "the cellfie task
  list" without knowing which specific `source` currently backs it. See
  `get_sources_for_task_list`/`register_task_source`. Optional fields use
  `COALESCE` in their upsert so a caller omitting one never silently
  erases a previously-recorded value -- **follow this pattern for any
  new optional column you add here**, it was a real bug once (see
  `decisions.md`).
- `enumeration_runs(source, task_id, model_id, status, n_routes,
  max_routes, hit_cap, truncated, cumulative_time_seconds, last_run_at)`
  -- one row per (source, task_id, model_id), the latest attempt's
  summary. `is_exhaustive` (computed, not stored) is `status=='optimal'
  and not hit_cap and not truncated` -- true only when the solver proved
  no further alternate route exists.
- `routes(route_id, source, task_id, model_id, reaction_set_hash,
  n_reactions, first_seen_at)` / `route_reactions(route_id, reaction_id,
  flux)` -- the actual enumerated routes, deduplicated by reaction-set
  hash, with each route's solved flux vector persisted so a later caller
  (e.g. the webapp's network view) never has to re-solve the LP.

Two currently-registered `task_list` families, both against Human-GEM
v2.0.1: `full` (`source="full"`, 257 tasks) and `cellfie`
(`source="cellfie_consensus_gurobi"`, 193 tasks -- the non-Gurobi
`cellfie_consensus` variant was retired, see `decisions.md`).

## A known, accepted model-version mismatch

The **classic** mapping-strategy CLI path loads `data/HumanGEM.xml.gz`
(an older bundled model) and `data/task_structure_matrix.tsv`. The
**context-aware** path loads whatever model `routes_human2.db` was
actually enumerated against (`data/HumanGEM_v201.xml`, newer). These are
*not* the same model version. This is intentional and already documented
in `enumerate_routes.py`'s own docstring -- don't "fix" it by silently
swapping one file for the other; if it ever needs reconciling, that's a
deliberate decision, not a bug fix.

## Webapp (`webapp/server.py`)

Local-only FastAPI app, not part of the installable package. Key things
that are *not* obvious from a first read:

- Every lookup is keyed by **`(source, task_id)`**, never `task_id`
  alone -- different sources reuse the same numeric task_id for
  unrelated tasks. Each source's own task-list file is resolved lazily
  via `task_sources` (`get_task_source`), not hardcoded.
- `/api/tasks/{source}/{task_id}/network` makes `dataset_id`/`sample`
  **optional**. With neither given, every route scores against an empty
  signal, which makes every route tie (score 0.0 everywhere) --
  `score_task` returns *all* tied routes as "winning", which is exactly
  the plain "show me the topology, no data" mode. Capped at
  `MAX_TOPOLOGY_PANELS` (15) because a task can have ~100 routes and
  building each one's graph triggers a flux solve on first view;
  rendering all of them eagerly would block the server's single-threaded
  event loop for the whole solve duration.
- The frontend paints that no-data mode with a neutral fill, not the
  normal evidence colors -- "no evidence" (red) would otherwise read as
  a negative finding instead of "nothing asked yet".
