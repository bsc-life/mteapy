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
../venv/bin/python -m uvicorn server:app --port 8765
```

Then open http://127.0.0.1:8765/.

Models come from the registry (`mteapy.registry`): each model folder's
`manifest.json` names its model file, task/route database and annotation
tables and is integrity-checked (a model that fails is listed as unusable,
with the reasons, instead of being hidden). By default the server looks in
`~/.mteapy/models` and mteapy's bundled models; set `MTEAPY_MODELS` to an
`os.pathsep`-separated list of model folders (or directories of them) to
look elsewhere. A model is loaded the first time it is selected -- a few
seconds -- and kept for the server's lifetime.

Nothing is solved at request time: every route's solved flux is stored in
the database.

## Using it

- **Model / Task list**: pick a model, then one of its task lists. Task
  ids are only unique within a list, so lists are never mixed in one table.
- **Analysis bar**: pick a method (TAS or CellFie), set its parameters (the form is
  generated from `mteapy.methods`, the same definitions the CLI uses), load a
  dataset, press **Run**. Every task of the selected list is scored against
  *every* sample in a background worker (progress bar, Cancel); when it
  finishes the Score column shows the selected sample, and switching sample is
  instant. The network view is scored with the run's own parameters. Changing
  the model, task list or dataset clears the results. CellFie's thresholds
  are computed over the whole loaded dataset (as in the CLI), so a sample's
  CellFie score depends on the other samples.
- **Save / Load results**: *Save results* downloads a self-contained JSON file
  (scores, method + parameters, model key and sha256, a fingerprint of the
  task list's routes, and the scoring signal restricted to the genes the routes
  reach -- the expression for TAS, the gene activity levels for CellFie).
  *Export TSV* is just the task x sample score matrix. *Load results* switches
  the view to the saved model/task list/method/parameters and shows the scores
  and networks without recomputing, but only accepts a file whose model file
  and routes are exactly the ones installed (otherwise it says why). A loaded
  TAS file's expression can be re-run; CellFie's cannot (its signal is the
  thresholded activity, not the expression).
- **View toggles** (top right) and the `x` on a panel show or hide the data bar,
  analysis bar, legend, task panel and detail panel -- the network takes the
  freed space: it is zoomed to the panel's width (its height follows the graph,
  at least the window's remaining height), and re-fits when panels toggle or the
  window resizes. Keyboard: `c`, `a`, `l`, `t`, `d` (ignored while typing in a
  field). The choice is remembered in the browser.

## What's static vs. per-sample

- **Static (per model + task list, computed once and cached for the
  server's lifetime)**: enumerated routes, each reaction's GPR-derived
  candidate complexes, EC/BiGG annotations (Human-GEM only), each route's
  stored flux vector, and which metabolites are the task's declared IN/OUT
  boundary.
- **Per-sample (cheap, recomputed live)**: each candidate complex's score
  (min over its genes' expression), each reaction's evidence classification,
  each route's aggregate score, and which route(s) win for that sample.

## API

- `GET /api/models` -- every discovered model with its task lists
  (`available: false` if its model file wasn't found next to the database).
  Reads the databases only; never loads a model.
- `GET /api/models/{model}/task_lists/{task_list}/tasks` -- the valid tasks
  of a list with route counts (this is what first loads the model).
- `POST /api/datasets` (multipart file) -- upload an expression table, get
  back a `dataset_id` and its sample columns.
- `GET /api/datasets/{id}` -- sample columns for an existing dataset (used
  for the bundled `gtex-example` dataset without re-uploading).
- `GET /api/methods` -- every scoring method with its parameters (name, type,
  default, choices, help).
- `POST /api/runs` `{model, task_list, dataset_id, method, params}` -- validate
  and start a run (202). 400 for a bad method/parameter or a dataset sharing no
  genes with the model; `GET /api/runs/{id}` -- status and progress;
  `DELETE /api/runs/{id}` -- cancel; `GET /api/runs/{id}/results` -- the
  tasks x samples score matrix (with completeness and tie flags).
- `GET /api/models/{model}/task_lists/{task_list}/tasks/{task_id}/network
  [?run_id=...&sample=...]` -- the winning route(s)' network graph for one
  task, scored with that run's parameters (all route variants, unscored, when
  no run is given; a task flagged invalid gets 409).

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
