"""Scoring a task list with a method, and the self-contained results file (the "bundle") that records it.

One place for what the CLI (``run-mtea analyze ... --bundle``) and the web
visualizer both do, so a bundle written by either is read by the other:

* `score_run` scores every task of a routes database against every sample
  (TAS, CellFie or TIDE, always the context-aware way) and returns the table the
  visualizer shows;
* `build_bundle` / `write_bundle` / `read_bundle` turn a finished run into a JSON
  file (optionally gzip-compressed, picked by the ``.gz`` suffix) holding the
  scores, the method and its parameters, which model and routes they were scored
  on, and the scoring signal restricted to the genes the routes can reach -- all
  the visualizer needs to redraw any task's network without recomputing.

The file deliberately holds no per-task network or evidence: those are rebuilt
from the signal with the installed model, so the scoring rules live in one place
and a file can never disagree with the code that reads it.
"""

from __future__ import annotations

import gzip
import importlib.metadata
import json
import os
from datetime import datetime, timezone

import pandas as pd

from mteapy import taskdb
from mteapy.cellfie import calculate_GAL
from mteapy.context_scoring import score_tasks_report
from mteapy.tide import compute_TIDE

RESULTS_FORMAT = "mteapy-results"
RESULTS_FORMAT_VERSION = 1


def scoring_conventions(method_id: str, params: dict) -> tuple[str, str]:
    """(aggregation, or_func) a method scores routes with."""
    if method_id == "TAS":
        return params["aggregation"], params["or_func"]
    if method_id == "CellFie":
        return "mean", "max"
    if method_id == "TIDE":
        return "mean", params["or_func"]
    raise ValueError(f"no scoring defined for method {method_id!r}")


def signal_kind(method_id: str) -> str:
    return {"CellFie": "gene_activity_levels", "TIDE": "log_fold_change"}.get(method_id, "expression")


def scoring_signal(method_id: str, params: dict, data: pd.DataFrame):
    """(gene signal, aggregation, or_func) that a method scores routes with.

    TAS projects the expression itself; CellFie first turns it into gene
    activity levels (thresholds are taken over the whole dataset, so this is
    computed once per run, never per sample); TIDE's `data` already is the
    (masked) log-fold-change column to score."""
    aggregation, or_func = scoring_conventions(method_id, params)
    if method_id == "CellFie":
        gal = calculate_GAL(
            data, thresh_type=params["threshold_type"], local_thresh_type=params["local_threshold_type"],
            minmaxmean_thresh_type=params["minmaxmean_threshold_type"], upper_bound=params["upper_bound"],
            lower_bound=params["lower_bound"], global_thresh_type=params["global_threshold_type"],
            global_value=params["global_value"], log_transformed=params["log_transformed"],
        )
        # A gene that is zero everywhere has threshold 0 -> 0/0; it simply has no activity.
        # (A gene the data has no value for stays missing, which is a different thing.)
        return gal.fillna(0.0).where(data.notna()), aggregation, or_func
    return data, aggregation, or_func


def tide_signal(table: pd.DataFrame, lfc_col: str, missing_genes: str) -> pd.DataFrame:
    """The one-column frame TIDE scores: `lfc_col` of an (already masked) table, with a gene
    the table has no value for kept missing, or read as no change when ``missing_genes="zero"``."""
    values = table[lfc_col].astype(float)
    if missing_genes == "zero":
        values = values.fillna(0.0)
    return pd.DataFrame({lfc_col: values})


def score_run(model, routes, complex_cache, descriptions, method_id: str, params: dict, data: pd.DataFrame,
              progress=None, on_permutations=None, pvalues: dict | None = None) -> dict:
    """Score every task of `routes` against every column of `data` with a method.

    `data` is what the method reads: an expression matrix (TAS, CellFie) or a
    one-column log-fold-change frame (TIDE). `progress(done, total)` is called
    per sample (and per batch of permutations for TIDE); `on_permutations(total)`
    just before TIDE starts them. `pvalues` ({task_id: p}) skips TIDE's permutations
    when they were already computed.

    Returns ``{samples, tasks, frames, signal, aggregation, or_func}``; each task
    is ``{task_id, description, n_routes, scores, complete, tied[, pvalues]}``."""
    signal, aggregation, or_func = scoring_signal(method_id, params, data)
    missing_genes = params["missing_genes"]
    if method_id == "TIDE":
        # The real scores are one cheap pass; the permutations are the run.
        frames = score_tasks_report(routes, model, signal, aggregation=aggregation, or_func=or_func,
                                    complex_cache=complex_cache, missing_genes=missing_genes)
        if pvalues is None:
            if on_permutations:
                on_permutations(params["n_permutations"])
            found = compute_TIDE(
                signal, params["lfc_col"], None, model, or_func, params["n_permutations"], params["n_cpus"],
                False, params["random_seed"], mapping_strategy="context-aware", tasks_routes=routes,
                complex_cache=complex_cache, permutation_strategy=params["permutation_strategy"],
                progress=progress, missing_genes=missing_genes,
            )
            pvalues = dict(zip(found["task_id"], (float(p) for p in found["pvalue"])))
    else:
        frames = score_tasks_report(routes, model, signal, aggregation=aggregation, or_func=or_func,
                                    progress=progress, complex_cache=complex_cache, missing_genes=missing_genes)
    tasks = []
    for task_id in frames["scores"].index:
        task = {
            "task_id": task_id, "description": descriptions.get(task_id, ""), "n_routes": len(routes[task_id]),
            "scores": [float(x) for x in frames["scores"].loc[task_id]],
            "complete": [bool(x) for x in frames["complete"].loc[task_id]],
            "tied": [bool(x) for x in frames["tied"].loc[task_id]],
        }
        if pvalues is not None:
            task["pvalues"] = [pvalues[task_id]]
        tasks.append(task)
    return {"samples": list(signal.columns), "tasks": tasks, "frames": frames, "signal": signal,
            "aggregation": aggregation, "or_func": or_func}


def relevant_genes(model, routes: dict, signal: pd.DataFrame) -> list[str]:
    """Genes of `signal` that any route can reach through a GPR -- all the scoring
    (and the network view) ever looks at, a few thousand of ~20k."""
    reactions = {r for task_routes in routes.values() for rs in task_routes.values() for r in rs}
    genes = {g.id for rid in reactions for g in model.reactions.get_by_id(rid).genes}
    return [g for g in signal.index if g in genes]


def build_bundle(*, model, model_info: dict, task_list: str, routes_fingerprint: str, routes: dict,
                 method_id: str, params: dict, result: dict, warnings=(), source: dict | None = None) -> dict:
    """The results file of a finished `score_run`.

    `model_info` is ``{key, name, version, sha256}`` of the model file the routes were
    enumerated on (what a reader checks before trusting the file); `source` is optional
    free provenance (e.g. the command line that produced it)."""
    signal = result["signal"]
    genes = relevant_genes(model, routes, signal)
    sub = signal.loc[genes, result["samples"]]
    try:
        version = importlib.metadata.version("mteapy")
    except importlib.metadata.PackageNotFoundError:
        version = None
    bundle = {
        "format": RESULTS_FORMAT, "format_version": RESULTS_FORMAT_VERSION,
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"), "mteapy_version": version,
        "model": {k: model_info.get(k) for k in ("key", "name", "version", "sha256")},
        "task_list": {"name": task_list, "routes_fingerprint": routes_fingerprint},
        "method": {"id": method_id, "params": params},
        "warnings": list(warnings),
        "samples": result["samples"],
        "tasks": result["tasks"],
        "signal": {"kind": signal_kind(method_id), "genes": genes, "samples": result["samples"],
                   "values": [[None if x != x else float(f"{x:.6g}") for x in row] for row in sub.to_numpy()]},
    }
    if source:
        bundle["source"] = source
    return bundle


def dumps_bundle(bundle: dict) -> bytes:
    return json.dumps(bundle, allow_nan=False).encode("utf-8")


def write_bundle(path: str, bundle: dict) -> str:
    """Write `bundle` to `path`: gzip-compressed when it ends in ``.gz``. Returns the path."""
    data = dumps_bundle(bundle)
    if str(path).endswith(".gz"):
        data = gzip.compress(data, compresslevel=6, mtime=0)
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    with open(path, "wb") as fh:
        fh.write(data)
    return path


def read_bundle(path: str) -> dict:
    """Read a results file, compressed or not (by content, not by name)."""
    with open(path, "rb") as fh:
        data = fh.read()
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    return json.loads(data.decode("utf-8"))


def model_identity(model_file: str, spec: str | None = None) -> dict:
    """``{key, name, version, sha256}`` of a model file: the registered model with the same
    hash when there is one (so the web app, which looks models up by key, recognises it),
    otherwise just the hash and the file's name."""
    from mteapy import registry     # local: registry imports taskdb, which imports nothing of ours

    sha = taskdb.model_file_sha256(model_file)
    for entry in registry.discover_models():
        if entry.sha256 == sha:
            return {"key": entry.key, "name": entry.name, "version": entry.version, "sha256": sha}
    return {"key": spec, "name": os.path.basename(model_file), "version": None, "sha256": sha}
