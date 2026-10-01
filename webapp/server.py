"""Local GUI for context-aware metabolic task scoring (mteapy).

This is the generalization of the earlier GTEx-only route-visualizer
prototype into an actual CellFie/TIDE-style tool: a user uploads their own
gene expression data (any CSV/TSV of gene_id + one or more sample columns),
the server scores every task in the routes database against whichever
sample they pick (mteapy.context_scoring, live -- no GTEx dependency), and
the existing D3/dagre network view renders whichever route(s) won for a
task the user clicks on.

Deliberately a local app, not a hosted Claude Artifact: scoring a *route*
against a sample is cheap pure arithmetic, but building a route's flux
vector the first time (mteapy.network.compute_route_fluxes) is an LP solve
over a genome-scale model copy, which needs COBRApy and a real solver --
not something that runs in a browser. Run with:

    uvicorn server:app --reload --port 8765   (from this directory)

then open http://127.0.0.1:8765/.
"""

from __future__ import annotations

import io
import uuid
from pathlib import Path

import pandas as pd
from cobra.io import read_sbml_model
from fastapi import FastAPI, File, HTTPException, UploadFile
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from mteapy.context_scoring import build_complex_cache, score_task
from mteapy.network import build_route_graph, compute_route_fluxes, load_bigg_ids, resolve_boundary_ids
from mteapy.routes import (
    connect, get_task_source, latest_model_id, list_tasks, load_route_fluxes, load_task_routes, save_route_fluxes,
)
from mteapy.task_model import build_metabolite_lookup
from mteapy.tasks import parse_task_file

MTEAPY_ROOT = Path(__file__).resolve().parent.parent
MTEAPY_DATA = MTEAPY_ROOT / "src" / "mteapy" / "data"
# Human-GEM and the GTEx example dataset are sibling checkouts of the whole
# workspace, not bundled into mteapy -- see PROVENANCE.md for why the
# model/routes DB moved into mteapy's own package data but these didn't.
WORKSPACE = MTEAPY_ROOT.parent

MODEL_PATH = MTEAPY_DATA / "HumanGEM_v201.xml"
DB_PATH = MTEAPY_DATA / "routes_human2.db"
METABOLITES_TSV = WORKSPACE / "Human-GEM" / "model" / "metabolites.tsv"
REACTIONS_TSV = WORKSPACE / "Human-GEM" / "model" / "reactions.tsv"
GTEX_PATH = WORKSPACE / "metabolic-variability" / "data" / "raw" / "gtex" / "GTEx_v10_gene_median_tpm.gct.gz"

app = FastAPI(title="mteapy task visualizer")

print("Loading model (this takes a few seconds)...")
MODEL = read_sbml_model(str(MODEL_PATH))
MET_LOOKUP = build_metabolite_lookup(MODEL)
MET_BIGG, RXN_BIGG = load_bigg_ids(str(METABOLITES_TSV), str(REACTIONS_TSV))
DB = connect(str(DB_PATH))
MODEL_ID = latest_model_id(DB)

# Different sources (e.g. "full" vs "cellfie_consensus_gurobi") read from
# different task-list files and reuse the same numeric task_id for
# unrelated tasks -- task "1" is a different task under each source. Every
# task lookup below is keyed by (source, task_id), never task_id alone;
# each source's task file is parsed lazily, on first use, via whatever
# path `task_sources` recorded for it at enumeration time (see
# `mteapy.cmds.enumerate_routes`/`register_task_source`).
_TASKS_BY_SOURCE: dict[str, dict[str, object]] = {}


def _tasks_for_source(source: str) -> dict[str, object]:
    if source not in _TASKS_BY_SOURCE:
        row = get_task_source(DB, source)
        if row is None:
            raise HTTPException(404, f"Unknown source {source!r} (no task_sources entry)")
        task_file = row["file_path"]
        path = Path(task_file)
        if not path.is_absolute():
            path = WORKSPACE / path
        _TASKS_BY_SOURCE[source] = {t.id: t for t in parse_task_file(str(path))}
    return _TASKS_BY_SOURCE[source]


print(f"Ready: model_id={MODEL_ID}.")

# Per-(source, task)-model structural caches -- independent of any sample,
# computed once and reused for every dataset/sample scored against that task.
_ROUTES_CACHE: dict[tuple[str, str], dict[int, frozenset[str]]] = {}
_COMPLEX_CACHE: dict[tuple[str, str], dict] = {}
_FLUX_CACHE: dict[tuple[str, str, int], dict[str, float]] = {}
_BOUNDARY_CACHE: dict[tuple[str, str], tuple[set[str], set[str]]] = {}

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


def get_task(source: str, task_id: str):
    task = _tasks_for_source(source).get(task_id)
    if task is None:
        raise HTTPException(404, f"Task {task_id!r} not found in source {source!r}'s task file")
    return task


def get_routes(source: str, task_id: str) -> dict[int, frozenset[str]]:
    key = (source, task_id)
    if key not in _ROUTES_CACHE:
        _ROUTES_CACHE[key] = load_task_routes(DB, source, task_id, MODEL_ID)
    return _ROUTES_CACHE[key]


def get_complex_cache(source: str, task_id: str, routes: dict[int, frozenset[str]]) -> dict:
    key = (source, task_id)
    if key not in _COMPLEX_CACHE:
        all_reactions = set().union(*routes.values()) if routes else set()
        _COMPLEX_CACHE[key] = build_complex_cache(MODEL, all_reactions)
    return _COMPLEX_CACHE[key]


def get_flux(source: str, task_id: str, route_id: int, reactions: frozenset[str]) -> dict[str, float]:
    """A route's flux never depends on any sample, so it's looked up from
    the routes DB first (persisted at enumeration time, or by an earlier
    visualizer request -- see `mteapy.routes.load_route_fluxes`/
    `save_route_fluxes`). Only on a genuine miss (an older route enumerated
    before flux persistence existed, and never viewed since) does this fall
    back to solving the LP fresh, immediately saving the result so that
    solve never has to happen again for this route, even across server
    restarts."""
    key = (source, task_id, route_id)
    if key in _FLUX_CACHE:
        return _FLUX_CACHE[key]

    stored = load_route_fluxes(DB, route_id)
    if len(stored) == len(reactions):
        _FLUX_CACHE[key] = stored
        return stored

    task = get_task(source, task_id)
    fluxes = compute_route_fluxes(MODEL, task, set(reactions))
    save_route_fluxes(DB, route_id, fluxes)
    _FLUX_CACHE[key] = fluxes
    return fluxes


def get_boundary_ids(source: str, task_id: str) -> tuple[set[str], set[str]]:
    key = (source, task_id)
    if key not in _BOUNDARY_CACHE:
        task = get_task(source, task_id)
        _BOUNDARY_CACHE[key] = resolve_boundary_ids(task, MET_LOOKUP)
    return _BOUNDARY_CACHE[key]


def get_dataset(dataset_id: str) -> pd.DataFrame:
    df = _DATASETS.get(dataset_id)
    if df is None:
        raise HTTPException(404, f"Unknown dataset {dataset_id!r} (upload it first)")
    return df


@app.get("/api/tasks")
async def api_tasks():
    tasks = list_tasks(DB, MODEL_ID)
    return {"tasks": tasks, "gtex_dataset_id": GTEX_DATASET_ID}


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


@app.get("/api/datasets/{dataset_id}/scores")
async def api_scores(dataset_id: str, sample: str):
    df = get_dataset(dataset_id)
    if sample not in df.columns:
        raise HTTPException(400, f"Unknown sample {sample!r}; available: {list(df.columns)}")
    gene_dict = df[sample].to_dict()

    results = []
    for row in list_tasks(DB, MODEL_ID):
        source, task_id = row["source"], row["task_id"]
        if task_id not in _tasks_for_source(source):
            continue  # a route exists but this task id isn't in that source's parsed task file
        routes = get_routes(source, task_id)
        complex_cache = get_complex_cache(source, task_id, routes)
        report = score_task(routes, complex_cache, gene_dict, task_id=task_id)
        results.append({
            "source": source,
            "task_id": task_id,
            "description": row["description"],
            "n_routes": row["n_routes"],
            "score": report.score,
            "is_tied": report.is_tied,
            "is_complete": report.is_complete,
            "winning_route_ids": list(report.winning_route_ids),
        })
    return {"sample": sample, "results": results}


MAX_TOPOLOGY_PANELS = 15


@app.get("/api/tasks/{source}/{task_id}/network")
async def api_network(source: str, task_id: str, dataset_id: str | None = None, sample: str | None = None):
    """A task's route network. With no `dataset_id`/`sample`, scores every
    route against an empty signal (gene_dict={}) -- every reaction reports
    `no_evidence`/`no_gpr` and every route ties at score 0, which means
    `score_task` returns *all* routes as "winning": exactly the plain
    topology view (every enumerated route variant, no data-driven
    winner) this mode is for. The frontend distinguishes this from a real
    all-tied result by simply not having asked for a sample.

    That all-tied set is capped at `MAX_TOPOLOGY_PANELS`: a task can have
    up to ~100 enumerated routes, and building each one's network graph
    means a flux solve on first view (`get_flux`'s cache miss path) --
    rendering all of them eagerly, synchronously, inside one request would
    both be a poor "browse the routes" UX (100 tabs to click through) and,
    worse, block this single-process server's event loop for the whole
    fallback-solve duration, freezing every other request (including
    trivial ones) until it finishes. A real scored result never hits this
    cap in practice (ties over genuine evidence are rare -- see
    docs/GTEX_ROUTE_SCORING_FINDINGS.md), so this only ever bites the
    topology-only, no-sample case."""
    if (dataset_id is None) != (sample is None):
        raise HTTPException(400, "dataset_id and sample must be given together, or both omitted")

    if dataset_id is not None:
        df = get_dataset(dataset_id)
        if sample not in df.columns:
            raise HTTPException(400, f"Unknown sample {sample!r}; available: {list(df.columns)}")
        gene_dict = df[sample].to_dict()
    else:
        gene_dict = {}

    task = get_task(source, task_id)
    routes = get_routes(source, task_id)
    if not routes:
        raise HTTPException(404, f"No enumerated routes for task {task_id!r} (source {source!r})")
    complex_cache = get_complex_cache(source, task_id, routes)
    report = score_task(routes, complex_cache, gene_dict, task_id=task_id)
    input_ids, output_ids = get_boundary_ids(source, task_id)

    shown_route_ids = report.winning_route_ids
    truncated = False
    if dataset_id is None and len(shown_route_ids) > MAX_TOPOLOGY_PANELS:
        shown_route_ids = shown_route_ids[:MAX_TOPOLOGY_PANELS]
        truncated = True

    panels = []
    for route_id in shown_route_ids:
        reactions = routes[route_id]
        fluxes = get_flux(source, task_id, route_id, reactions)
        graph = build_route_graph(MODEL, set(reactions), gene_dict, fluxes,
                                   MET_BIGG, RXN_BIGG, input_ids, output_ids)
        panels.append({"route_id": route_id, **graph})

    return {
        "source": source,
        "task_id": task_id,
        "task_description": task.description,
        "tissue": sample,
        "has_data": dataset_id is not None,
        "score": report.score,
        "is_complete": report.is_complete,
        "n_routes_total": len(report.winning_route_ids),
        "truncated": truncated,
        "panels": panels,
    }


@app.get("/")
def index():
    return FileResponse(str(Path(__file__).parent / "static" / "index.html"))


app.mount("/static", StaticFiles(directory=str(Path(__file__).parent / "static")), name="static")
