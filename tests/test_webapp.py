"""End-to-end tests of the webapp's API over a toy model folder."""

import importlib
import io
import sys
import time
from pathlib import Path

import pytest

pytest.importorskip("fastapi")
pytest.importorskip("httpx")

from fastapi.testclient import TestClient  # noqa: E402

from mteapy import taskdb  # noqa: E402
from mteapy.enumeration import enumerate_alternate_routes  # noqa: E402
from test_registry import _make_model_folder  # noqa: E402

WEBAPP = Path(__file__).resolve().parent.parent / "webapp"
EXPRESSION = "geneID\ts1\ts2\ng1\t10\t1\ng2\t10\t1\ng3\t1\t8\ng4\t5\t0\ng_tr\t7\t7\n"


@pytest.fixture
def client(toy_model, tmp_path, monkeypatch):
    def tweak(manifest, directory):
        manifest.pop("annotations")
    directory = _make_model_folder(tmp_path, toy_model, tweak=tweak)
    conn = taskdb.connect(str(directory / "routes.db"))
    pk = taskdb.get_task_pk(conn, "ToyTasks", "1")
    found = enumerate_alternate_routes(toy_model, taskdb.load_task(conn, pk), max_routes=5)
    taskdb.record_enumeration_result(conn, pk, [r.reactions for r in found], [r.fluxes for r in found],
                                     status="optimal", max_routes=5, hit_cap=False, truncated=False,
                                     elapsed_seconds=0.0)
    conn.close()
    monkeypatch.setenv("MTEAPY_MODELS", str(tmp_path))
    monkeypatch.syspath_prepend(str(WEBAPP))
    for name in ("server", "runs"):
        sys.modules.pop(name, None)
    server = importlib.import_module("server")
    with TestClient(server.app) as c:
        c.server = server
        yield c


def _upload(client, text=EXPRESSION):
    r = client.post("/api/datasets", files={"file": ("x.tsv", io.BytesIO(text.encode()), "text/plain")})
    assert r.status_code == 200, r.text
    return r.json()["dataset_id"]


def _run(client, dataset_id, **overrides):
    body = {"model": "toy-1", "task_list": "ToyTasks", "dataset_id": dataset_id, "method": "TAS", "params": {}}
    body.update(overrides)
    return client.post("/api/runs", json=body)


def _wait(client, run_id, timeout=30):
    deadline = time.time() + timeout
    while time.time() < deadline:
        info = client.get(f"/api/runs/{run_id}").json()
        if info["status"] not in ("queued", "running"):
            return info
        time.sleep(0.05)
    raise AssertionError("run did not finish")


def test_models_and_methods_are_listed(client):
    [model] = client.get("/api/models").json()["models"]
    assert model["key"] == "toy-1" and model["available"] and model["problems"] == []
    assert [t["name"] for t in model["task_lists"]] == ["ToyTasks"]
    described = {m["id"]: m for m in client.get("/api/methods").json()["methods"]}
    assert list(described) == ["TAS", "CellFie"]
    assert {p["name"] for p in described["TAS"]["params"]} == {"aggregation", "or_func"}
    upper = next(p for p in described["CellFie"]["params"] if p["name"] == "global_value")
    assert upper["when"] == [["threshold_type", "global"]]


def test_tasks_of_a_list(client):
    r = client.get("/api/models/toy-1/task_lists/ToyTasks/tasks").json()
    assert [(t["task_id"], t["n_routes"]) for t in r["tasks"]] == [("1", 2)]
    assert client.get("/api/models/toy-1/task_lists/Nope/tasks").status_code == 404
    assert client.get("/api/models/nope/task_lists/ToyTasks/tasks").status_code == 404


def test_full_run_produces_the_score_matrix(client):
    run = _run(client, _upload(client))
    assert run.status_code == 202
    info = _wait(client, run.json()["run_id"])
    assert info["status"] == "done" and info["done"] == info["total"] == 2
    assert info["genes_matched"] == 5
    results = client.get(f"/api/runs/{info['run_id']}/results").json()
    assert results["samples"] == ["s1", "s2"]
    [task] = results["tasks"]
    assert task["task_id"] == "1" and task["scores"] == [10.0, 8.0] and task["n_routes"] == 2


def test_parameters_change_the_result_and_are_recorded(client):
    # One sample where the aggregation matters: route R1+R2 has reaction scores (10, 2) -> min 2, mean 6;
    # route R3 scores 4. min picks R3 (4 > 2); mean picks R1+R2 (6 > 4).
    ds = _upload(client, "geneID\ts1\ng1\t10\ng2\t2\ng3\t4\ng4\t1\ng_tr\t1\n")
    r_min = _wait(client, _run(client, ds).json()["run_id"])
    r_mean = _wait(client, _run(client, ds, params={"aggregation": "mean"}).json()["run_id"])
    assert r_mean["config"]["params"] == {"aggregation": "mean", "or_func": "max"}
    scores = lambda info: client.get(f"/api/runs/{info['run_id']}/results").json()["tasks"][0]["scores"]
    assert scores(r_min) == [4.0]
    assert scores(r_mean) == [6.0]
    # the network of each run is scored with that run's own aggregation
    base = "/api/models/toy-1/task_lists/ToyTasks/tasks/1/network"
    assert client.get(base, params={"run_id": r_min["run_id"], "sample": "s1"}).json()["score"] == 4.0
    assert client.get(base, params={"run_id": r_mean["run_id"], "sample": "s1"}).json()["score"] == 6.0


def test_network_uses_the_runs_sample_and_parameters(client):
    run = _wait(client, _run(client, _upload(client)).json()["run_id"])
    base = "/api/models/toy-1/task_lists/ToyTasks/tasks/1/network"
    topo = client.get(base).json()
    assert topo["has_data"] is False and len(topo["panels"]) == 2
    scored = client.get(base, params={"run_id": run["run_id"], "sample": "s1"}).json()
    assert scored["has_data"] is True and scored["score"] == 10.0 and len(scored["panels"]) == 1
    other = client.get(base, params={"run_id": run["run_id"], "sample": "s2"}).json()
    assert other["score"] == 8.0
    assert client.get(base, params={"run_id": run["run_id"]}).status_code == 400
    assert client.get(base, params={"run_id": run["run_id"], "sample": "zz"}).status_code == 400
    assert client.get(base, params={"run_id": "nope", "sample": "s1"}).status_code == 404


def test_run_validation_errors(client):
    ds = _upload(client)
    assert _run(client, ds, method="Nope").status_code == 400
    assert _run(client, ds, params={"aggregation": "sum"}).status_code == 400
    assert _run(client, ds, params={"bogus": 1}).status_code == 400
    assert _run(client, ds, task_list="Nope").status_code == 404
    assert _run(client, ds, model="nope").status_code == 404
    assert _run(client, "no-such-dataset").status_code == 404
    assert client.post("/api/runs", json={"model": "toy-1"}).status_code == 400
    wrong_ids = _upload(client, "geneID\ts1\nENSG1\t1\nENSG2\t2\n")
    r = _run(client, wrong_ids)
    assert r.status_code == 400 and "None of the dataset's" in r.json()["detail"]


def test_results_before_completion_and_unknown_runs(client):
    assert client.get("/api/runs/nope").status_code == 404
    assert client.get("/api/runs/nope/results").status_code == 404
    assert client.delete("/api/runs/nope").status_code == 404


def test_a_run_can_be_cancelled(client):
    runs = client.server.RUNS
    started = __import__("threading").Event()

    def slow(run):
        from runs import check_cancelled
        started.set()
        for _ in range(200):
            check_cancelled(run)
            time.sleep(0.02)
        return {}
    run = runs.submit({"model": "x"}, slow)
    assert started.wait(5)
    assert client.delete(f"/api/runs/{run.run_id}").status_code == 200
    assert _wait(client, run.run_id)["status"] == "cancelled"


def test_a_failing_run_reports_its_error(client):
    def boom(run):
        raise RuntimeError("kaput")
    run = client.server.RUNS.submit({"model": "x"}, boom)
    info = _wait(client, run.run_id)
    assert info["status"] == "error" and "kaput" in info["error"]
    assert client.get(f"/api/runs/{run.run_id}/results").status_code == 409


def test_cellfie_run_matches_the_library_and_scores_networks_on_gals(client, toy_model):
    import pandas as pd
    from mteapy.cellfie import calculate_CellFie_scores_context_aware, calculate_GAL

    expr = pd.read_csv(io.StringIO(EXPRESSION), sep="\t", index_col=0)
    conn = taskdb.connect(str(client.server.MODELS["toy-1"].db_path))
    routes = taskdb.load_task_list_routes(conn, "ToyTasks")
    expected, _ = calculate_CellFie_scores_context_aware(calculate_GAL(expr), routes, toy_model)

    run = _wait(client, _run(client, _upload(client), method="CellFie").json()["run_id"])
    assert run["status"] == "done"
    assert run["config"]["params"]["threshold_type"] == "local"
    [task] = client.get(f"/api/runs/{run['run_id']}/results").json()["tasks"]
    assert task["scores"] == pytest.approx(list(expected.loc["1"]))

    base = "/api/models/toy-1/task_lists/ToyTasks/tasks/1/network"
    scored = client.get(base, params={"run_id": run["run_id"], "sample": "s1"}).json()
    assert scored["score"] == pytest.approx(expected.loc["1", "s1"])


def test_cellfie_global_threshold_changes_scores_and_bad_values_are_rejected(client):
    ds = _upload(client)
    local = _wait(client, _run(client, ds, method="CellFie").json()["run_id"])
    glob = _wait(client, _run(client, ds, method="CellFie", params={"threshold_type": "global"}).json()["run_id"])
    scores = lambda info: client.get(f"/api/runs/{info['run_id']}/results").json()["tasks"][0]["scores"]
    assert scores(local) != scores(glob)
    assert _run(client, ds, method="CellFie", params={"threshold_type": "sideways"}).status_code == 400
    assert _run(client, ds, method="CellFie", params={"aggregation": "min"}).status_code == 400


# ---------------------------------------------------------------------------
# saving and loading results
# ---------------------------------------------------------------------------

def _finished_run(client, **overrides):
    return _wait(client, _run(client, _upload(client), **overrides).json()["run_id"])


def test_export_import_round_trip_tas(client):
    run = _finished_run(client)
    scores = client.get(f"/api/runs/{run['run_id']}/results").json()["tasks"]
    exported = client.get(f"/api/runs/{run['run_id']}/export")
    assert exported.status_code == 200 and "attachment" in exported.headers["content-disposition"]
    bundle = exported.json()
    assert bundle["format"] == "mteapy-results" and bundle["signal"]["kind"] == "expression"
    assert bundle["model"]["key"] == "toy-1" and bundle["method"] == {"id": "TAS", "params": run["config"]["params"]}
    assert "g_tr" not in bundle["signal"]["genes"]       # a gene no route reaches is not embedded

    loaded = client.post("/api/runs/import", json=bundle)
    assert loaded.status_code == 201, loaded.text
    info = loaded.json()
    assert info["status"] == "done" and info["config"]["loaded"] and info["samples"] == ["s1", "s2"]
    assert client.get(f"/api/runs/{info['run_id']}/results").json()["tasks"] == scores
    base = "/api/models/toy-1/task_lists/ToyTasks/tasks/1/network"
    for sample, expected in (("s1", 10.0), ("s2", 8.0)):
        assert client.get(base, params={"run_id": info["run_id"], "sample": sample}).json()["score"] == expected
    # the embedded expression is a usable dataset: re-running reproduces the saved scores
    again = _wait(client, _run(client, info["dataset_id"]).json()["run_id"])
    assert client.get(f"/api/runs/{again['run_id']}/results").json()["tasks"] == scores


def test_export_import_round_trip_cellfie_keeps_gals(client):
    run = _finished_run(client, method="CellFie")
    bundle = client.get(f"/api/runs/{run['run_id']}/export").json()
    assert bundle["signal"]["kind"] == "gene_activity_levels"
    info = client.post("/api/runs/import", json=bundle).json()
    assert info["dataset_id"] is None            # GALs depend on the whole dataset: not offered for re-running
    base = "/api/models/toy-1/task_lists/ToyTasks/tasks/1/network"
    original = client.get(base, params={"run_id": run["run_id"], "sample": "s1"}).json()["score"]
    assert client.get(base, params={"run_id": info["run_id"], "sample": "s1"}).json()["score"] == pytest.approx(original)


def test_tsv_export(client):
    run = _finished_run(client)
    r = client.get(f"/api/runs/{run['run_id']}/export", params={"format": "tsv"})
    rows = [line.split("\t") for line in r.text.strip().split("\n")]
    assert rows[0] == ["task_id", "description", "n_routes", "s1", "s2"]
    assert rows[1][0] == "1" and rows[1][2] == "2" and [float(x) for x in rows[1][3:]] == [10.0, 8.0]
    assert client.get(f"/api/runs/{run['run_id']}/export", params={"format": "xml"}).status_code == 400
    assert client.get("/api/runs/nope/export").status_code == 404


def test_import_refuses_files_that_do_not_match_this_installation(client):
    run = _finished_run(client)
    good = client.get(f"/api/runs/{run['run_id']}/export").json()

    def attempt(mutate):
        import copy
        bundle = copy.deepcopy(good)
        mutate(bundle)
        return client.post("/api/runs/import", json=bundle)

    assert attempt(lambda b: b.update(format="other")).status_code == 400
    assert attempt(lambda b: b.update(format_version=99)).status_code == 400
    assert attempt(lambda b: b["model"].update(key="elsewhere")).status_code == 404
    r = attempt(lambda b: b["model"].update(sha256="0" * 64))
    assert r.status_code == 409 and "differs" in r.json()["detail"]
    r = attempt(lambda b: b["task_list"].update(routes_fingerprint="0" * 64))
    assert r.status_code == 409 and "routes" in r.json()["detail"]
    assert attempt(lambda b: b["method"]["params"].update(aggregation="sum")).status_code == 400
    assert attempt(lambda b: b["tasks"][0].update(task_id="zz")).status_code == 400
    assert attempt(lambda b: b["tasks"][0].update(n_routes=7)).status_code == 400
    assert attempt(lambda b: b["tasks"][0].update(scores=[1.0])).status_code == 400
    assert attempt(lambda b: b.pop("signal")).status_code == 400
