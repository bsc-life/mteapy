# mteapy task visualizer

A local, CellFie/TIDE-style GUI for context-aware metabolic task scoring:
upload any gene expression table, score every task in the routes database
live, and inspect the winning route's network (with GPR/complex evidence,
EC/BiGG ids, flux direction, and task input/output coloring) by clicking a
row.

Not tied to GTEx: the bundled GTEx tissue matrix is offered as a one-click
example dataset, but "Upload expression" accepts any CSV/TSV with a gene-id
column and one or more sample columns.

## Run it

```
cd webapp
../metabolic-task-builder/.venv/bin/python3 -m uvicorn server:app --reload --port 8765
```

Then open http://127.0.0.1:8765/.

The first request that renders a given task's network takes a few seconds
per route (an LP solve over the genome-scale model to recover that route's
flux vector) -- this is cached in memory after that, so re-viewing the same
task (even against a different sample) is near-instant. Re-scoring all 220
tasks against a newly picked sample currently takes several seconds too
(pure arithmetic, but 220 tasks x up to 100 routes each); the whole score
table is recomputed on every sample change.

## What's static vs. per-sample

- **Static (per model+task, computed once and cached for the server's
  lifetime)**: enumerated routes (`data/routes_human2.db`, already built),
  each reaction's GPR-derived candidate complexes, EC/BiGG annotations,
  each route's solved flux vector, and which metabolites are the task's
  declared IN/OUT boundary.
- **Per-sample (cheap, recomputed live)**: each candidate complex's score
  (min over its genes' expression), each reaction's evidence classification,
  each route's aggregate score, and which route(s) win for that sample.

This split is why the tool can score an arbitrary uploaded dataset live
without needing COBRApy or a solver at request time for the scoring step
itself -- only the one-time-per-route flux solve needs them.

## API

- `GET /api/tasks` -- every task with a route count.
- `POST /api/datasets` (multipart file) -- upload an expression table, get
  back a `dataset_id` and its sample columns.
- `GET /api/datasets/{id}` -- sample columns for an existing dataset (used
  for the bundled `gtex-example` dataset without re-uploading).
- `GET /api/datasets/{id}/scores?sample=...` -- score every task against
  one sample column.
- `GET /api/datasets/{id}/network/{task_id}?sample=...` -- the winning
  route(s)' network graph for one task/sample, in the same
  `{task_id, task_description, tissue, panels: [...]}` shape the earlier
  static-Artifact prototype used.

## Frontend: how a network result gets rendered

The rendering code is split the same way the Python side is split between
`cobra_netgraph` (builds the graph) and `mteapy.network` (decorates it with
task-scoring fields), just kept as two files in this same directory rather
than two packages:

- **`static/netgraph-viz.js`** is a generic dagre+D3 renderer. It knows
  nothing about tasks, genes, GPRs, or flux -- it only knows how to lay
  out, draw, zoom/pan, and dispatch clicks on a `{nodes, edges}` graph,
  given plain callbacks that turn a node/edge into a shape, a color, a
  label. It started as a standalone cross-project package
  (`netgraph-viz`), but came back in-repo: a future gap-filling/
  reconstruction visualization tool would need a different-enough
  interaction model (curating a large draft network vs. inspecting one
  small enumerated route) that sharing a renderer across both would mean
  either bloating it with options to please both, or one project fighting
  it. `cobra-netgraph`, by contrast, stayed a shared Python package,
  because building the same bipartite graph from a COBRA model genuinely
  is one shared problem, not two.
- **`static/app.js`** is the only mteapy-specific piece: it fetches a
  result from this server's API, then supplies netgraph-viz.js with the
  callbacks that know what `evidence`, `score`, `io` and `flux` mean.

The actual load path, end to end:

1. Clicking a task row calls `selectTask(taskId)`, which does
   `GET /api/datasets/{id}/network/{taskId}?sample=...` (see API below) and
   stores the JSON response as `currentData` -- `{task_id,
   task_description, tissue, score, is_complete, panels: [{route_id,
   nodes, edges}, ...]}` (more than one panel only when several routes are
   tied for the winning score).
2. `showPanel(0)` picks the first panel and calls `renderPanel(host, panel,
   1000, 800)`.
3. `renderPanel` builds two D3 scales *from that panel's own data* --a
   node-score color scale and a flux-width scale-- then calls
   `NetGraphViz.render(host, panel, {...})`, passing those scales in via
   closures inside the required callbacks (`nodeFill`, `edgeColor`,
   `edgeWidth`, ...). This is why the scales are built fresh per panel
   view: colors/widths are only ever comparable within one rendered route,
   not across the whole task list.
4. netgraph-viz lays the graph out with dagre and draws it; clicking a node
   calls back into `onNodeClick`, which is `renderDetail(d)` in `app.js` --
   this is where the reaction/metabolite detail panel (equation, EC/BiGG
   ids, the GPR/complex bipartite mini-diagram) is built, entirely outside
   netgraph-viz.

If you change how a task's network should look (new field, new coloring
rule), the change belongs in `app.js`'s callbacks, not in
`netgraph-viz.js` -- that file should stay ignorant of what a "reaction" or
an "evidence" value is.

## Known rough edges

- Single-process, in-memory dataset/cache storage -- restarting the server
  drops uploaded datasets (the bundled GTEx example reloads automatically).
- No auth; it's meant to run on `127.0.0.1` for one local user.
- Scoring all 220 tasks per sample-change is not yet optimized (it's
  correct, just not fast) -- fine for occasional interactive use, would
  want batching/caching work before scoring many samples in a loop.
