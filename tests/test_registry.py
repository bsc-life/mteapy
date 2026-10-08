import shutil

import pytest
from cobra.io import write_sbml_model

from mteapy import registry, taskdb


def _make_db(directory, model, name="toy", version="1", db_name="routes.db", write_model=True):
    path = str(directory / "model.xml")
    write_sbml_model(model, path)
    conn = taskdb.create_db(str(directory / db_name))
    taskdb.register_model(conn, model, name, taskdb.model_file_sha256(path), version)
    conn.close()
    if not write_model:
        (directory / "model.xml").unlink()
    return str(directory / db_name)


def test_discovers_database_with_matching_model_file(toy_model, tmp_path):
    db = _make_db(tmp_path, toy_model)
    [entry] = registry.discover_models([str(tmp_path)])
    assert (entry.key, entry.name, entry.version) == ("toy-1", "toy", "1")
    assert entry.db_path == db and entry.model_path == str(tmp_path / "model.xml") and entry.available


def test_model_without_its_file_is_listed_but_unavailable(toy_model, tmp_path):
    _make_db(tmp_path, toy_model, write_model=False)
    [entry] = registry.discover_models([str(tmp_path)])
    assert entry.model_path is None and not entry.available


def test_a_different_xml_next_to_the_database_does_not_match(toy_model, tmp_path):
    _make_db(tmp_path, toy_model)
    (tmp_path / "model.xml").write_text("<sbml/>")
    [entry] = registry.discover_models([str(tmp_path)])
    assert entry.model_path is None


def test_scans_one_level_of_subdirectories_and_skips_legacy_databases(toy_model, tmp_path):
    sub = tmp_path / "ToyModel"
    sub.mkdir()
    _make_db(sub, toy_model)
    import sqlite3
    sqlite3.connect(tmp_path / "legacy.db").execute("CREATE TABLE routes (x)").connection.commit()
    keys = [e.key for e in registry.discover_models([str(tmp_path)])]
    assert keys == ["toy-1"]


def test_colliding_keys_are_made_unique(toy_model, tmp_path):
    for d in ("a", "b"):
        (tmp_path / d).mkdir()
        _make_db(tmp_path / d, toy_model)
    keys = [e.key for e in registry.discover_models([str(tmp_path)])]
    assert len(set(keys)) == 2


# ---------------------------------------------------------------------------
# Manifests
# ---------------------------------------------------------------------------

import json

TASKS_TXT = ("\tID\tDESCRIPTION\tSHOULD FAIL\tIN\tIN LB\tIN UB\tOUT\tOUT LB\tOUT UB\tEQU\tEQU LB\tEQU UB\tCHANGED RXN\t"
             "CHANGED LB\tCHANGED UB\tPRINT FLUX\tCOMMENTS\n\t1\tA to C\t\tA[c]\t1\t1\tC[c]\t1\t1\t\t\t\t\t\t\t\t\n")


def _make_model_folder(parent, toy_model, name="toy", version="1", folder="ToyModel", tweak=None):
    """A complete, valid model folder: model, schema-v2 db with one imported task list, manifest."""
    directory = parent / folder
    (directory / "tasks").mkdir(parents=True)
    model_path = directory / "model.xml"
    write_sbml_model(toy_model, str(model_path))
    sha = taskdb.model_file_sha256(str(model_path))
    (directory / "tasks" / "t.txt").write_text(TASKS_TXT)
    conn = taskdb.create_db(str(directory / "routes.db"))
    model_id = taskdb.register_model(conn, toy_model, name, sha, version)
    taskdb.import_task_list(conn, toy_model, model_id, "ToyTasks", str(directory / "tasks" / "t.txt"),
                            origin_sha256=taskdb.task_source_sha256(str(directory / "tasks" / "t.txt")))
    conn.close()
    manifest = {
        "manifest_version": 1, "name": name, "version": version, "description": "a toy", "license": "MIT",
        "gene_id_type": "symbol",
        "model": {"file": "model.xml", "format": "sbml", "sha256": sha},
        "database": {"file": "routes.db"},
        "task_lists": [{"name": "ToyTasks", "file": "tasks/t.txt",
                        "sha256": taskdb.task_source_sha256(str(directory / "tasks" / "t.txt"))}],
        "annotations": {"metabolites": "m.tsv"},
    }
    if tweak:
        tweak(manifest, directory)
    (directory / "manifest.json").write_text(json.dumps(manifest))
    return directory


def test_manifest_model_is_discovered_with_its_metadata(toy_model, tmp_path):
    directory = _make_model_folder(tmp_path, toy_model)
    [entry] = registry.discover_models([str(tmp_path)])
    assert entry.available and entry.problems == ()
    assert (entry.key, entry.gene_id_type, entry.license) == ("toy-1", "symbol", "MIT")
    assert entry.model_path == str(directory / "model.xml") and entry.db_path == str(directory / "routes.db")
    assert entry.annotations == {"metabolites": str(directory / "m.tsv")}
    assert [t["name"] for t in entry.task_lists] == ["ToyTasks"]


def test_a_model_folder_can_be_given_directly(toy_model, tmp_path):
    directory = _make_model_folder(tmp_path, toy_model)
    assert [e.key for e in registry.discover_models([str(directory)])] == ["toy-1"]


def test_modified_model_file_is_flagged_not_hidden(toy_model, tmp_path):
    directory = _make_model_folder(tmp_path, toy_model)
    (directory / "model.xml").write_text("tampered")
    [entry] = registry.discover_models([str(tmp_path)])
    assert not entry.available and entry.model_path is None
    assert any("sha256" in p for p in entry.problems)


def test_database_built_against_another_model_is_flagged(toy_model, tmp_path):
    def tweak(manifest, directory):
        manifest["model"]["sha256"] = "0" * 64
    _make_model_folder(tmp_path, toy_model, tweak=tweak)
    [entry] = registry.discover_models([str(tmp_path)])
    assert not entry.available
    assert any("different model" in p for p in entry.problems)


def test_task_list_missing_from_database_or_with_different_hash_is_flagged(toy_model, tmp_path):
    def tweak(manifest, directory):
        manifest["task_lists"].append({"name": "Ghost", "file": "tasks/t.txt"})
        manifest["task_lists"][0]["sha256"] = "f" * 64
    _make_model_folder(tmp_path, toy_model, tweak=tweak)
    [entry] = registry.discover_models([str(tmp_path)])
    assert any("Ghost" in p and "not in the database" in p for p in entry.problems)
    assert any("ToyTasks" in p and "hash differs" in p for p in entry.problems)


def test_missing_database_is_flagged(toy_model, tmp_path):
    directory = _make_model_folder(tmp_path, toy_model)
    (directory / "routes.db").unlink()
    [entry] = registry.discover_models([str(tmp_path)])
    assert any("database not found" in p for p in entry.problems)


@pytest.mark.parametrize("mutate, expected", [
    (lambda m: m.pop("name"), "'name'"),
    (lambda m: m.update(manifest_version=2), "manifest_version"),
    (lambda m: m["model"].update(file="/etc/passwd"), "model.file"),
    (lambda m: m["database"].update(file="../outside.db"), "database.file"),
    (lambda m: m["model"].update(sha256="xyz"), "model.sha256"),
    (lambda m: m["task_lists"].append(dict(m["task_lists"][0])), "duplicate task list"),
    (lambda m: m.update(annotations={"x": "/abs.tsv"}), "annotations"),
])
def test_validate_manifest_rejects_malformed_or_escaping_paths(mutate, expected):
    manifest = {"manifest_version": 1, "name": "n", "version": "1", "model": {"file": "m.xml", "sha256": "a" * 64},
                "database": {"file": "r.db"}, "task_lists": [{"name": "T", "file": "t.txt"}]}
    assert registry.validate_manifest(manifest) == []
    mutate(manifest)
    assert any(expected in p for p in registry.validate_manifest(manifest))


def test_unreadable_manifest_is_listed_with_the_reason(tmp_path):
    (tmp_path / "Broken").mkdir()
    (tmp_path / "Broken" / "manifest.json").write_text("{not json")
    [entry] = registry.discover_models([str(tmp_path)])
    assert not entry.available and "cannot read manifest" in entry.problems[0]


def test_resolve_model_by_key_folder_and_default(toy_model, tmp_path):
    directory = _make_model_folder(tmp_path, toy_model)
    assert registry.resolve_model("toy-1", [str(tmp_path)]).directory == str(directory)
    assert registry.resolve_model(str(directory)).key == "toy-1"
    with pytest.raises(KeyError, match="available: \\['toy-1'\\]"):
        registry.resolve_model("nope", [str(tmp_path)])
    assert registry.resolve_model().key == registry.DEFAULT_MODEL_KEY      # the bundled model


def test_resolve_model_refuses_an_unusable_model_with_reasons(toy_model, tmp_path):
    directory = _make_model_folder(tmp_path, toy_model)
    (directory / "model.xml").write_text("tampered")
    with pytest.raises(ValueError, match="sha256"):
        registry.resolve_model("toy-1", [str(tmp_path)])


def test_search_paths_come_from_the_environment_when_set(monkeypatch, tmp_path):
    monkeypatch.setenv(registry.SEARCH_PATHS_ENV, f"{tmp_path}{__import__('os').pathsep}/elsewhere")
    assert registry.default_search_paths() == [str(tmp_path), "/elsewhere"]
    monkeypatch.delenv(registry.SEARCH_PATHS_ENV)
    assert registry.default_search_paths()[-1] == registry.BUNDLED_MODELS_DIR


def test_bundled_model_is_complete_and_consistent():
    entry = registry.resolve_model()
    assert entry.available and entry.gene_id_type == "ensembl"
    assert {t["name"] for t in entry.task_lists} == {"HumanGEM-Full", "CellFie"}
    assert set(entry.annotations) == {"metabolites", "reactions"}
