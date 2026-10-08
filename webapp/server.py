"""Local GUI for context-aware metabolic task scoring (mteapy).

A user uploads their own gene expression data (any CSV/TSV of gene_id + one
or more sample columns), picks a model and one of its task lists, and the
server scores every task of that list against whichever sample they pick
(mteapy.context_scoring, live), while the D3/dagre network view renders
whichever route(s) won for a task the user clicks on.

Models come from the registry (`mteapy.registry`): each model folder's
manifest names its model file, task/route database (`mteapy.taskdb`) and
annotation tables, and is integrity-checked. Each model is loaded lazily,
the first time it is used, into a `ModelContext` that owns its model,
database connection and caches.

No solver is needed: every route's solved flux is stored in the database.
The network view still needs COBRApy and the genome-scale model (gene and
metabolite annotations), which is why this is a local Python app and not
something that runs in a browser. Run with:

    uvicorn server:app --port 8765   (from this directory)

then open http://127.0.0.1:8765/. ``MTEAPY_MODELS`` (an os.pathsep-separated
list of model folders / directories of them) overrides where models are
looked for; the default is ~/.mteapy/models plus the bundled models.
"""

from __future__ import annotations

import importlib.metadata
import io
import json
import os
import sqlite3
import threading
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from cobra.io import read_sbml_model
from fastapi import Body, FastAPI, File, HTTPException, UploadFile
from fastapi.responses import FileResponse, Response
from fastapi.staticfiles import StaticFiles

from mteapy import methods, registry, taskdb
from mteapy.cellfie import calculate_GAL
from mteapy.context_scoring import build_complex_cache, score_task, score_tasks_report
from mteapy.network import build_route_graph, load_bigg_ids, resolve_boundary_ids
from mteapy.task_model import build_metabolite_lookup
from runs import RunManager, check_cancelled

# The GTEx example dataset is a sibling checkout of the whole workspace,
# not bundled into mteapy.
WORKSPACE = Path(__file__).resolve().parent.parent.parent
GTEX_PATH = WORKSPACE / "metabolic-variability" / "data" / "raw" / "gtex" / "GTEx_v10_gene_median_tpm.gct.gz"

app = FastAPI(title="mteapy task visualizer")

MODELS: dict[str, registry.ModelEntry] = {e.key: e for e in registry.discover_models()}
print(f"Models found: {[k + ('' if e.available else ' (unusable)') for k, e in MODELS.items()] or 'none'}")


@dataclass
class ModelContext:
    """Everything the server holds for one loaded model. Caches are keyed by
    (task_list, task_id) -- different task lists reuse the same numeric ids
    for unrelated tasks -- and are independent of any sample, so they're
    computed once and reused for every dataset/sample scored against a task."""

    entry: registry.ModelEntry
    model: object
    met_lookup: dict
    met_bigg: dict | None
    rxn_bigg: dict | None
    # One shared connection, used only from async handlers (never concurrently).
    db: object
    task_pks: dict = field(default_factory=dict)
    routes: dict = field(default_factory=dict)
    complexes: dict = field(default_factory=dict)
    fluxes: dict = field(default_factory=dict)
    boundaries: dict = field(default_factory=dict)
    # Whole-task-list structures a run needs; built on first use, then shared.
    list_routes: dict = field(default_factory=dict)      # task_list -> {task_id: {route_id: reactions}}
    list_complexes: dict = field(default_factory=dict)   # task_list -> complex cache over all its reactions
    list_descriptions: dict = field(default_factory=dict)  # task_list -> {task_id: description}
    lock: threading.Lock = field(default_factory=threading.Lock)


_CONTEXTS: dict[str, ModelContext] = {}


def get_context(model_key: str) -> ModelContext:
    """The loaded context for `model_key`, loading it on first use (a few seconds)."""
    if model_key in _CONTEXTS:
        return _CONTEXTS[model_key]
    entry = MODELS.get(model_key)
    if entry is None:
        raise HTTPException(404, f"Unknown model {model_key!r}; available: {list(MODELS)}")
    if not entry.available:
        raise HTTPException(409, f"Model {model_key!r} cannot be used: "
                                 f"{'; '.join(entry.problems) or 'its model file was not found next to the database'}")
    print(f"Loading model {model_key} ...", flush=True)
    model = read_sbml_model(entry.model_path)
    met_bigg = rxn_bigg = None
    if entry.annotations.get("metabolites") and entry.annotations.get("reactions"):
        met_bigg, rxn_bigg = load_bigg_ids(entry.annotations["metabolites"], entry.annotations["reactions"])
    db = taskdb.connect(entry.db_path, check_same_thread=False)
    _CONTEXTS[model_key] = ctx = ModelContext(entry, model, build_metabolite_lookup(model), met_bigg, rxn_bigg, db)
    print(f"Ready: {model_key}.", flush=True)
    return ctx


def _task_pk(ctx: ModelContext, task_list: str, task_id: str) -> int:
    key = (task_list, task_id)
    if key not in ctx.task_pks:
        pk = taskdb.get_task_pk(ctx.db, task_list, task_id)
        if pk is None:
            raise HTTPException(404, f"Task {task_id!r} not found in task list {task_list!r}")
        ctx.task_pks[key] = pk
    return ctx.task_pks[key]


def get_task(ctx: ModelContext, task_list: str, task_id: str):
    try:
        return taskdb.load_task(ctx.db, _task_pk(ctx, task_list, task_id))
    except ValueError as exc:  # flagged invalid at import: it has no definition to build a network from
        raise HTTPException(409, str(exc)) from None


def get_routes(ctx: ModelContext, task_list: str, task_id: str) -> dict[int, frozenset[str]]:
    key = (task_list, task_id)
    if key not in ctx.routes:
        ctx.routes[key] = taskdb.load_task_routes(ctx.db, _task_pk(ctx, task_list, task_id))
    return ctx.routes[key]


def get_complex_cache(ctx: ModelContext, task_list: str, task_id: str, routes: dict[int, frozenset[str]]) -> dict:
    key = (task_list, task_id)
    if key not in ctx.complexes:
        all_reactions = set().union(*routes.values()) if routes else set()
        ctx.complexes[key] = build_complex_cache(ctx.model, all_reactions)
    return ctx.complexes[key]


def get_flux(ctx: ModelContext, route_id: int) -> dict[str, float]:
    """A route's flux never depends on any sample, and is stored for every
    route in the database, so this is a plain lookup -- nothing is solved."""
    if route_id not in ctx.fluxes:
        ctx.fluxes[route_id] = taskdb.load_route_fluxes(ctx.db, route_id)
    return ctx.fluxes[route_id]


def get_boundary_ids(ctx: ModelContext, task_list: str, task_id: str) -> tuple[set[str], set[str]]:
    key = (task_list, task_id)
    if key not in ctx.boundaries:
        ctx.boundaries[key] = resolve_boundary_ids(get_task(ctx, task_list, task_id), ctx.met_lookup)
    return ctx.boundaries[key]


# In-memory uploaded datasets: dataset_id -> DataFrame (index=gene_id, columns=samples).
_DATASETS: dict[str, pd.DataFrame] = {}


def _load_gtex_example() -> str:
    if not GTEX_PATH.exists():
        return ""
    df = pd.read_csv(GTEX_PATH, sep="\t", skiprows=2, index_col=0)
    df = df.drop(columns=["Description"])
    df.index = df.index.str.replace(r"\.\d+$", "", regex=True)
    df = df.groupby(df.index).mean()
    dataset_id = "gtex-example"
    _DATASETS[dataset_id] = df
    print(f"Loaded bundled GTEx example dataset ({df.shape[1]} tissues).")
    return dataset_id


GTEX_DATASET_ID = _load_gtex_example()


def get_dataset(dataset_id: str) -> pd.DataFrame:
    df = _DATASETS.get(dataset_id)
    if df is None:
        raise HTTPException(404, f"Unknown dataset {dataset_id!r} (upload it first)")
    return df


def _check_task_list(ctx: ModelContext, task_list: str) -> None:
    try:
        taskdb.get_task_list(ctx.db, task_list)
    except KeyError as exc:
        raise HTTPException(404, exc.args[0]) from None


@app.get("/api/models")
async def api_models():
    """Every discovered model with its task lists. Reads the databases only --
    it never loads a model, so listing stays instant."""
    models = []
    for key, entry in MODELS.items():
        task_lists = []
        if os.path.isfile(entry.db_path):
            try:
                conn = taskdb.connect(entry.db_path)
            except (ValueError, sqlite3.DatabaseError):
                conn = None
            if conn is not None:
                try:
                    task_lists = taskdb.list_task_lists(conn)
                finally:
                    conn.close()
        models.append({
            "key": key, "name": entry.name, "version": entry.version, "available": entry.available,
            "problems": list(entry.problems), "description": entry.description, "license": entry.license,
            "gene_id_type": entry.gene_id_type, "loaded": key in _CONTEXTS, "task_lists": task_lists,
        })
    return {"models": models, "gtex_dataset_id": GTEX_DATASET_ID}


@app.get("/api/models/{model}/task_lists/{task_list}/tasks")
async def api_tasks(model: str, task_list: str):
    ctx = get_context(model)
    _check_task_list(ctx, task_list)
    return {"model": model, "task_list": task_list,
            "tasks": [t for t in taskdb.list_tasks(ctx.db, task_list) if t["valid"]]}


@app.post("/api/datasets")
async def api_upload_dataset(file: UploadFile = File(...)):
    """Accept a CSV/TSV of gene expression: first column = gene id (Ensembl
    or matching the model's gene naming), remaining columns = one sample
    each. Any CellFie/TIDE-style flat expression table works."""
    raw = await file.read()
    sep = "\t" if file.filename.endswith((".tsv", ".txt", ".gct")) else ","
    try:
        df = pd.read_csv(io.BytesIO(raw), sep=sep, index_col=0)
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(400, f"Could not parse {file.filename!r}: {exc}") from None
    df = df.select_dtypes(include="number")
    if df.empty:
        raise HTTPException(400, "No numeric sample columns found in the uploaded file")
    dataset_id = str(uuid.uuid4())
    _DATASETS[dataset_id] = df
    return {"dataset_id": dataset_id, "samples": list(df.columns), "n_genes": len(df)}


@app.get("/api/datasets/{dataset_id}")
async def api_dataset_info(dataset_id: str):
    df = get_dataset(dataset_id)
    return {"dataset_id": dataset_id, "samples": list(df.columns), "n_genes": len(df)}


@app.get("/api/methods")
async def api_methods():
    """Every scoring method with its parameters (name, type, default, choices, help)."""
    return {"methods": methods.describe_methods()}


RUNS = RunManager()


def _task_list_structures(ctx: ModelContext, task_list: str):
    """(routes, complex cache, descriptions) for a whole task list, built once per model.
    Runs in the worker thread, so it uses its own read-only connection rather
    than the one the request handlers share."""
    with ctx.lock:
        if task_list not in ctx.list_routes:
            conn = taskdb.connect(ctx.entry.db_path)
            try:
                routes = taskdb.load_task_list_routes(conn, task_list, min_routes=1)
                descriptions = {t["task_id"]: t["description"] for t in taskdb.list_tasks(conn, task_list)}
            finally:
                conn.close()
            all_reactions = sorted({r for rs in routes.values() for reactions in rs.values() for r in reactions})
            ctx.list_routes[task_list] = routes
            ctx.list_complexes[task_list] = build_complex_cache(ctx.model, all_reactions)
            ctx.list_descriptions[task_list] = descriptions
        return ctx.list_routes[task_list], ctx.list_complexes[task_list], ctx.list_descriptions[task_list]


def _scoring_signal(method_id: str, params: dict, expr: pd.DataFrame):
    """(gene signal, aggregation, or_func) that a method scores routes with.

    TAS projects the expression itself; CellFie first turns it into gene
    activity levels (thresholds are taken over the whole dataset, so this is
    computed once per run, never per sample)."""
    if method_id == "TAS":
        return expr, params["aggregation"], params["or_func"]
    if method_id == "CellFie":
        gal = calculate_GAL(
            expr, thresh_type=params["threshold_type"], local_thresh_type=params["local_threshold_type"],
            minmaxmean_thresh_type=params["minmaxmean_threshold_type"], upper_bound=params["upper_bound"],
            lower_bound=params["lower_bound"], global_thresh_type=params["global_threshold_type"],
            global_value=params["global_value"], log_transformed=params["log_transformed"],
        )
        # A gene that is zero everywhere has threshold 0 -> 0/0; it simply has no activity.
        return gal.fillna(0.0), "mean", "max"
    raise ValueError(f"no scoring defined for method {method_id!r}")


def _execute_run(run, ctx: ModelContext, expr: pd.DataFrame):
    cfg = run.config
    routes, complex_cache, descriptions = _task_list_structures(ctx, cfg["task_list"])
    check_cancelled(run)
    run.total = expr.shape[1]
    signal, aggregation, or_func = _scoring_signal(cfg["method"], cfg["params"], expr)
    check_cancelled(run)

    def progress(done, total):
        run.done, run.total = done, total
        check_cancelled(run)

    frames = score_tasks_report(routes, ctx.model, signal, aggregation=aggregation, or_func=or_func,
                                progress=progress, complex_cache=complex_cache)
    samples = list(expr.columns)
    tasks = [{
        "task_id": task_id, "description": descriptions.get(task_id, ""), "n_routes": len(routes[task_id]),
        "scores": [float(x) for x in frames["scores"].loc[task_id]],
        "complete": [bool(x) for x in frames["complete"].loc[task_id]],
        "tied": [bool(x) for x in frames["tied"].loc[task_id]],
    } for task_id in frames["scores"].index]
    return {"samples": samples, "tasks": tasks, "frames": frames, "signal": signal,
            "aggregation": aggregation, "or_func": or_func}


@app.post("/api/runs", status_code=202)
async def api_start_run(body: dict = Body(...)):
    """Start scoring `task_list` of `model` against every sample of `dataset_id`.

    Body: ``{model, task_list, dataset_id, method, params}``. Returns the
    run's info at once; poll ``GET /api/runs/{run_id}`` for progress.
    """
    missing = [k for k in ("model", "task_list", "dataset_id", "method") if not body.get(k)]
    if missing:
        raise HTTPException(400, f"missing field(s): {missing}")
    try:
        method = methods.get_method(body["method"])
        params = methods.validate_params(method, body.get("params"))
    except (KeyError, ValueError) as exc:
        raise HTTPException(400, exc.args[0]) from None
    ctx = get_context(body["model"])
    _check_task_list(ctx, body["task_list"])
    expr = get_dataset(body["dataset_id"])

    model_genes = {g.id for g in ctx.model.genes}
    matched = len(model_genes & set(expr.index))
    if matched == 0:
        raise HTTPException(400, f"None of the dataset's {len(expr)} genes is a gene of this model "
                                 f"(the model uses {ctx.entry.gene_id_type or 'its own'} gene ids, e.g. "
                                 f"{sorted(model_genes)[:3]}); check the id column")
    config = {"model": body["model"], "task_list": body["task_list"], "dataset_id": body["dataset_id"],
              "method": method.id, "params": params}
    run = RUNS.submit(config, lambda r: _execute_run(r, ctx, expr), genes_matched=matched)
    return {**run.info(), "model_genes": len(model_genes)}


def _get_run(run_id: str):
    run = RUNS.get(run_id)
    if run is None:
        raise HTTPException(404, f"Unknown run {run_id!r}")
    return run


@app.get("/api/runs/{run_id}")
async def api_run_status(run_id: str):
    return _get_run(run_id).info()


@app.delete("/api/runs/{run_id}")
async def api_cancel_run(run_id: str):
    run = RUNS.cancel(run_id)
    if run is None:
        raise HTTPException(404, f"Unknown run {run_id!r}")
    return run.info()


@app.get("/api/runs/{run_id}/results")
async def api_run_results(run_id: str):
    run = _get_run(run_id)
    if run.status != "done":
        raise HTTPException(409, f"Run {run_id} is {run.status}" + (f": {run.error}" if run.error else ""))
    return {**run.info(), "samples": run.result["samples"], "tasks": run.result["tasks"]}


RESULTS_FORMAT = "mteapy-results"
RESULTS_FORMAT_VERSION = 1


def _signal_kind(method_id: str) -> str:
    return "gene_activity_levels" if method_id == "CellFie" else "expression"


def _relevant_genes(ctx: ModelContext, task_list: str, signal: pd.DataFrame) -> list[str]:
    """Genes of `signal` that any route of `task_list` can reach through a GPR --
    all the scoring (and the network view) ever looks at, a few thousand of ~20k."""
    routes, _, _ = _task_list_structures(ctx, task_list)
    reactions = {r for task_routes in routes.values() for rs in task_routes.values() for r in rs}
    genes = {g.id for rid in reactions for g in ctx.model.reactions.get_by_id(rid).genes}
    return [g for g in signal.index if g in genes]


@app.get("/api/runs/{run_id}/export")
async def api_export_run(run_id: str, format: str = "json"):
    """A finished run as a file: ``json`` is the self-contained bundle that
    ``POST /api/runs/import`` reads back (scores, parameters, the model it was
    scored on, and the scoring signal restricted to the model-relevant genes);
    ``tsv`` is just the task x sample score matrix."""
    run = _get_run(run_id)
    if run.status != "done":
        raise HTTPException(409, f"Run {run_id} is {run.status}")
    cfg, result = run.config, run.result
    stem = f"mteapy_{cfg['method']}_{cfg['model']}_{cfg['task_list']}"
    if format == "tsv":
        lines = ["\t".join(["task_id", "description", "n_routes", *result["samples"]])]
        for t in result["tasks"]:
            desc = (t["description"] or "").replace("\t", " ")
            lines.append("\t".join([t["task_id"], desc, str(t["n_routes"]), *(repr(x) for x in t["scores"])]))
        return Response("\n".join(lines) + "\n", media_type="text/tab-separated-values",
                        headers={"Content-Disposition": f'attachment; filename="{stem}.tsv"'})
    if format != "json":
        raise HTTPException(400, "format must be 'json' or 'tsv'")

    ctx = get_context(cfg["model"])
    signal = result["signal"]
    genes = _relevant_genes(ctx, cfg["task_list"], signal)
    sub = signal.loc[genes, result["samples"]]
    try:
        version = importlib.metadata.version("mteapy")
    except importlib.metadata.PackageNotFoundError:
        version = None
    bundle = {
        "format": RESULTS_FORMAT, "format_version": RESULTS_FORMAT_VERSION,
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"), "mteapy_version": version,
        "model": {"key": ctx.entry.key, "name": ctx.entry.name, "version": ctx.entry.version,
                  "sha256": ctx.entry.sha256},
        "task_list": {"name": cfg["task_list"],
                      "routes_fingerprint": taskdb.task_list_routes_fingerprint(ctx.db, cfg["task_list"])},
        "method": {"id": cfg["method"], "params": cfg["params"]},
        "samples": result["samples"],
        "tasks": result["tasks"],
        "signal": {"kind": _signal_kind(cfg["method"]), "genes": genes, "samples": result["samples"],
                   "values": [[float(f"{x:.6g}") for x in row] for row in sub.to_numpy()]},
    }
    return Response(json.dumps(bundle, allow_nan=False), media_type="application/json",
                    headers={"Content-Disposition": f'attachment; filename="{stem}.json"'})


@app.post("/api/runs/import", status_code=201)
async def api_import_run(bundle: dict = Body(...)):
    """Load a bundle written by ``/export`` as a finished run. It is only
    accepted against the very same model file and route sets it was scored on;
    anything else is refused with the reason."""
    if bundle.get("format") != RESULTS_FORMAT:
        raise HTTPException(400, f"Not an mteapy results file (format {bundle.get('format')!r})")
    if bundle.get("format_version") != RESULTS_FORMAT_VERSION:
        raise HTTPException(400, f"Unsupported results format version {bundle.get('format_version')!r} "
                                 f"(this server reads version {RESULTS_FORMAT_VERSION})")
    try:
        model_info, tl_info, method_info = bundle["model"], bundle["task_list"], bundle["method"]
        samples, tasks, sig = bundle["samples"], bundle["tasks"], bundle["signal"]
        task_list = tl_info["name"]
        method = methods.get_method(method_info["id"])
        params = methods.validate_params(method, method_info.get("params"))
    except (KeyError, TypeError, ValueError) as exc:
        raise HTTPException(400, f"Malformed results file: {exc.args[0] if exc.args else exc!r}") from None

    if model_info.get("key") not in MODELS:
        raise HTTPException(404, f"The results were scored on model {model_info.get('key')!r}, which is not "
                                 f"installed here (available: {list(MODELS)})")
    ctx = get_context(model_info["key"])
    if model_info.get("sha256") != ctx.entry.sha256:
        raise HTTPException(409, f"Model {ctx.entry.key} differs from the one the results were scored on "
                                 f"(file sha256 {str(model_info.get('sha256'))[:12]}..., installed "
                                 f"{ctx.entry.sha256[:12]}...)")
    _check_task_list(ctx, task_list)
    if tl_info.get("routes_fingerprint") != taskdb.task_list_routes_fingerprint(ctx.db, task_list):
        raise HTTPException(409, f"The routes stored for task list {task_list!r} changed since these results "
                                 f"were saved, so they no longer describe this database")
    routes, _, descriptions = _task_list_structures(ctx, task_list)
    try:
        for t in tasks:
            if t["task_id"] not in routes or len(routes[t["task_id"]]) != t["n_routes"]:
                raise ValueError(f"task {t['task_id']!r} does not match the database")
            if not (len(t["scores"]) == len(t["complete"]) == len(t["tied"]) == len(samples)):
                raise ValueError(f"task {t['task_id']!r} does not have one value per sample")
        if list(sig["samples"]) != list(samples):
            raise ValueError("signal samples differ from the result samples")
        signal = pd.DataFrame(sig["values"], index=sig["genes"], columns=sig["samples"], dtype=float)
        if signal.shape != (len(sig["genes"]), len(samples)):
            raise ValueError("signal matrix has the wrong shape")
    except (KeyError, TypeError, ValueError) as exc:
        raise HTTPException(400, f"Malformed results file: {exc.args[0] if exc.args else exc!r}") from None

    aggregation, or_func = (("mean", "max") if method.id == "CellFie"
                            else (params["aggregation"], params["or_func"]))
    dataset_id = None
    if sig.get("kind") == "expression":
        # The embedded subset is a valid expression matrix for the model's genes, so it can be re-run.
        dataset_id = str(uuid.uuid4())
        _DATASETS[dataset_id] = signal
    config = {"model": ctx.entry.key, "task_list": task_list, "dataset_id": dataset_id, "method": method.id,
              "params": params, "loaded": True, "saved_at": bundle.get("created")}
    result = {"samples": list(samples), "tasks": tasks, "signal": signal,
              "aggregation": aggregation, "or_func": or_func}
    run = RUNS.add_finished(config, result, genes_matched=len(signal))
    return {**run.info(), "dataset_id": dataset_id, "n_tasks": len(tasks), "samples": list(samples)}


MAX_TOPOLOGY_PANELS = 15


@app.get("/api/models/{model}/task_lists/{task_list}/tasks/{task_id}/network")
async def api_network(model: str, task_list: str, task_id: str, run_id: str | None = None,
                      sample: str | None = None):
    """A task's route network for one sample of a finished run (scored with that
    run's own parameters, so it matches the results table). With no
    `run_id`/`sample`, scores every
    route against an empty signal (gene_dict={}) -- every reaction reports
    `no_evidence`/`no_gpr` and every route ties at score 0, which means
    `score_task` returns *all* routes as "winning": exactly the plain
    topology view (every enumerated route variant, no data-driven
    winner) this mode is for. The frontend distinguishes this from a real
    all-tied result by simply not having asked for a sample.

    That all-tied set is capped at `MAX_TOPOLOGY_PANELS`: a task can have
    up to ~100 enumerated routes, and building each one's network graph
    is real work (annotation, GPR scoring, layout payload) -- rendering all
    of them eagerly inside one request would be a poor "browse the routes"
    UX (100 tabs to click through) and would hold this single-process
    server's event loop for the whole duration. A real scored result never
    hits this cap in practice (ties over genuine evidence are rare -- see
    docs/GTEX_ROUTE_SCORING_FINDINGS.md), so this only ever bites the
    topology-only, no-sample case."""
    if (run_id is None) != (sample is None):
        raise HTTPException(400, "run_id and sample must be given together, or both omitted")

    aggregation, or_func = "min", "max"
    if run_id is not None:
        run = _get_run(run_id)
        if run.status != "done":
            raise HTTPException(409, f"Run {run_id} is {run.status}")
        cfg = run.config
        if (cfg["model"], cfg["task_list"]) != (model, task_list):
            raise HTTPException(400, f"Run {run_id} scored {cfg['model']}/{cfg['task_list']}, "
                                     f"not {model}/{task_list}")
        df = run.result["signal"]
        if sample not in df.columns:
            raise HTTPException(400, f"Unknown sample {sample!r}; available: {list(df.columns)}")
        gene_dict = df[sample].to_dict()
        aggregation, or_func = run.result["aggregation"], run.result["or_func"]
    else:
        gene_dict = {}

    ctx = get_context(model)
    _check_task_list(ctx, task_list)
    task = get_task(ctx, task_list, task_id)
    routes = get_routes(ctx, task_list, task_id)
    if not routes:
        raise HTTPException(404, f"No enumerated routes for task {task_id!r} (task list {task_list!r})")
    complex_cache = get_complex_cache(ctx, task_list, task_id, routes)
    report = score_task(routes, complex_cache, gene_dict, aggregation=aggregation, task_id=task_id, or_func=or_func)
    input_ids, output_ids = get_boundary_ids(ctx, task_list, task_id)

    shown_route_ids = report.winning_route_ids
    truncated = False
    if run_id is None and len(shown_route_ids) > MAX_TOPOLOGY_PANELS:
        shown_route_ids = shown_route_ids[:MAX_TOPOLOGY_PANELS]
        truncated = True

    panels = []
    for route_id in shown_route_ids:
        reactions = routes[route_id]
        graph = build_route_graph(ctx.model, set(reactions), gene_dict, get_flux(ctx, route_id),
                                   ctx.met_bigg, ctx.rxn_bigg, input_ids, output_ids, or_func=or_func)
        panels.append({"route_id": route_id, **graph})

    return {
        "model": model,
        "task_list": task_list,
        "task_id": task_id,
        "task_description": task.description,
        "tissue": sample,
        "has_data": run_id is not None,
        "score": report.score,
        "is_complete": report.is_complete,
        "n_routes_total": len(report.winning_route_ids),
        "truncated": truncated,
        "panels": panels,
    }


@app.get("/")
def index():
    return FileResponse(str(Path(__file__).parent / "static" / "index.html"))


class _RevalidatingStatic(StaticFiles):
    """Static files that browsers must revalidate (cheap 304s), so an upgraded
    app never runs against a cached copy of its older scripts."""

    async def get_response(self, path, scope):
        response = await super().get_response(path, scope)
        response.headers["Cache-Control"] = "no-cache"
        return response


app.mount("/static", _RevalidatingStatic(directory=str(Path(__file__).parent / "static")), name="static")
