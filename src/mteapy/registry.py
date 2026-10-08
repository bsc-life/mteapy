"""The model registry: which models mteapy can work with, and where their files are.

A model lives in its own folder:

    HumanGEM/
      manifest.json          what the model is and where its files are
      HumanGEM_v201.xml      the model (SBML)
      routes.db              schema-v2 task/route database (mteapy.taskdb)
      tasks/*.txt            the task-list files the database was imported from
      annotations/*.tsv      optional BiGG/EC annotation tables

``manifest.json`` (all paths relative to the folder, never absolute or
escaping it)::

    {
      "manifest_version": 1,
      "name": "Human-GEM", "version": "2.0.1", "description": "...",
      "license": "CC-BY-4.0",
      "origin": {"repo": "...", "note": "..."},          # optional
      "gene_id_type": "ensembl",
      "model": {"file": "HumanGEM_v201.xml", "format": "sbml", "sha256": "..."},
      "database": {"file": "routes.db"},
      "task_lists": [{"name": "...", "file": "tasks/x.txt", "sha256": "...", "description": "..."}],
      "annotations": {"metabolites": "annotations/m.tsv", "reactions": "annotations/r.tsv"}   # optional
    }

Nothing model-specific is hard-coded in the code: names, file locations and
annotation tables all come from here. `discover_models` also still accepts a
bare schema-v2 database next to a model file whose hash it recorded (no
manifest) -- handy for scratch work -- but a manifest is the supported form.

A model that cannot be used (a missing or modified file, a database built
against a different model, ...) is still *listed*, with the reasons in
`ModelEntry.problems`, so callers can say why instead of silently hiding it.
"""

from __future__ import annotations

import json
import os
import re
import sqlite3
from dataclasses import dataclass, field

from mteapy import taskdb

MANIFEST_NAME = "manifest.json"
MANIFEST_VERSION = 1
BUNDLED_MODELS_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data", "models")
DEFAULT_MODEL_KEY = "Human-GEM-2.0.1"
SEARCH_PATHS_ENV = "MTEAPY_MODELS"

_hash_cache: dict[tuple[str, int, int], str] = {}


@dataclass(frozen=True)
class ModelEntry:
    key: str                              # URL/CLI-safe unique handle, e.g. "Human-GEM-2.0.1"
    name: str
    version: str | None
    sha256: str                           # the model file's expected hash
    db_path: str
    model_path: str | None                # None: the model file was not found
    directory: str | None = None          # the model folder (None for a bare database)
    description: str = ""
    license: str | None = None
    gene_id_type: str | None = None
    annotations: dict = field(default_factory=dict)      # {"metabolites": path, "reactions": path}
    task_lists: tuple = ()                # manifest's task-list entries (name, file, sha256, description)
    problems: tuple = ()                  # why this entry cannot be used; empty when it can

    @property
    def available(self) -> bool:
        return self.model_path is not None and not self.problems


def _cached_sha256(path: str) -> str:
    st = os.stat(path)
    cache_key = (path, st.st_mtime_ns, st.st_size)
    if cache_key not in _hash_cache:
        _hash_cache[cache_key] = taskdb.model_file_sha256(path)
    return _hash_cache[cache_key]


def _slug(name: str, version: str | None) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "-", f"{name}-{version}" if version else name).strip("-")


def default_search_paths() -> list[str]:
    """``MTEAPY_MODELS`` (os.pathsep-separated) if set, else the user's
    ``~/.mteapy/models`` (if it exists) followed by the bundled models."""
    env = [p for p in os.environ.get(SEARCH_PATHS_ENV, "").split(os.pathsep) if p]
    if env:
        return env
    user = os.path.join(os.path.expanduser("~"), ".mteapy", "models")
    return ([user] if os.path.isdir(user) else []) + [BUNDLED_MODELS_DIR]


# ---------------------------------------------------------------------------
# Manifests
# ---------------------------------------------------------------------------

def _safe_relative(path) -> bool:
    return (isinstance(path, str) and path != "" and not os.path.isabs(path)
            and ".." not in os.path.normpath(path).split(os.sep))


def validate_manifest(manifest) -> list[str]:
    """Structural problems with a parsed manifest (empty list: well-formed).
    Checks shape only -- not that the files exist or hash correctly."""
    problems: list[str] = []
    if not isinstance(manifest, dict):
        return ["manifest is not a JSON object"]
    if manifest.get("manifest_version") != MANIFEST_VERSION:
        problems.append(f"manifest_version must be {MANIFEST_VERSION}, got {manifest.get('manifest_version')!r}")
    if not isinstance(manifest.get("name"), str) or not manifest.get("name"):
        problems.append("'name' is required")
    version = manifest.get("version")
    if version is not None and not isinstance(version, str):
        problems.append("'version' must be a string")
    for section in ("model", "database"):
        entry = manifest.get(section)
        if not isinstance(entry, dict) or not _safe_relative(entry.get("file")):
            problems.append(f"'{section}.file' must be a relative path inside the model folder")
    model = manifest.get("model")
    if isinstance(model, dict) and not re.fullmatch(r"[0-9a-f]{64}", str(model.get("sha256", ""))):
        problems.append("'model.sha256' must be a 64-character hex digest")
    task_lists = manifest.get("task_lists", [])
    if not isinstance(task_lists, list):
        problems.append("'task_lists' must be a list")
    else:
        names = set()
        for tl in task_lists:
            if not isinstance(tl, dict) or not isinstance(tl.get("name"), str) or not _safe_relative(tl.get("file")):
                problems.append(f"each task list needs a 'name' and a relative 'file': {tl!r}")
                continue
            if tl["name"] in names:
                problems.append(f"duplicate task list name {tl['name']!r}")
            names.add(tl["name"])
    annotations = manifest.get("annotations", {})
    if not isinstance(annotations, dict) or not all(_safe_relative(v) for v in annotations.values()):
        problems.append("'annotations' must map names to relative paths")
    return problems


def _entry_from_manifest(directory: str) -> ModelEntry:
    path = os.path.join(directory, MANIFEST_NAME)
    base = os.path.basename(os.path.normpath(directory))
    try:
        manifest = json.load(open(path))
    except (OSError, ValueError) as exc:
        return ModelEntry(base, base, None, "", "", None, directory, problems=(f"cannot read manifest: {exc}",))
    problems = validate_manifest(manifest)
    if problems:
        return ModelEntry(_slug(str(manifest.get("name", base)), None), str(manifest.get("name", base)), None, "",
                          "", None, directory, problems=tuple(problems))

    def inside(rel):
        return os.path.join(directory, rel)

    name, version = manifest["name"], manifest.get("version")
    model_path, db_path = inside(manifest["model"]["file"]), inside(manifest["database"]["file"])
    expected = manifest["model"]["sha256"]

    if not os.path.isfile(model_path):
        problems.append(f"model file not found: {manifest['model']['file']}")
        model_path = None
    elif _cached_sha256(model_path) != expected:
        problems.append(f"model file {manifest['model']['file']} does not match the manifest's sha256 "
                        f"(modified or wrong version)")
        model_path = None

    if not os.path.isfile(db_path):
        problems.append(f"database not found: {manifest['database']['file']}")
    else:
        try:
            conn = taskdb.connect(db_path)
        except (ValueError, sqlite3.DatabaseError) as exc:
            problems.append(f"database unusable: {exc}")
        else:
            try:
                db_models = {sha for (sha,) in conn.execute("SELECT sha256 FROM models")}
                if expected not in db_models:
                    problems.append("the database was built against a different model than the manifest names")
                db_sha = dict(conn.execute("SELECT name, origin_sha256 FROM task_lists"))
            finally:
                conn.close()
            for tl in manifest.get("task_lists", []):
                if tl["name"] not in db_sha:
                    problems.append(f"task list {tl['name']!r} is in the manifest but not in the database")
                elif tl.get("sha256") and db_sha[tl["name"]] and tl["sha256"] != db_sha[tl["name"]]:
                    problems.append(f"task list {tl['name']!r}: the manifest's file hash differs from the one the "
                                    f"database was imported from")

    return ModelEntry(
        key=_slug(name, version), name=name, version=version, sha256=expected, db_path=db_path,
        model_path=model_path, directory=directory, description=manifest.get("description", ""),
        license=manifest.get("license"), gene_id_type=manifest.get("gene_id_type"),
        annotations={k: inside(v) for k, v in manifest.get("annotations", {}).items()},
        task_lists=tuple(manifest.get("task_lists", [])), problems=tuple(problems),
    )


# ---------------------------------------------------------------------------
# Discovery and resolution
# ---------------------------------------------------------------------------

def _entries_from_bare_database(db_path: str) -> list[ModelEntry]:
    try:
        conn = taskdb.connect(db_path)
    except (ValueError, sqlite3.DatabaseError):
        return []
    try:
        rows = conn.execute("SELECT name, version, sha256 FROM models ORDER BY model_id").fetchall()
    finally:
        conn.close()
    directory = os.path.dirname(db_path)
    candidates = [os.path.join(directory, f) for f in sorted(os.listdir(directory)) if f.endswith((".xml", ".sbml"))]
    by_sha = {_cached_sha256(p): p for p in candidates}
    return [ModelEntry(_slug(n, v), n, v, sha, db_path, by_sha.get(sha)) for n, v, sha in rows]


def _candidate_directories(path: str) -> list[str]:
    """`path` itself if it holds a manifest, else its immediate subdirectories that do."""
    if os.path.isfile(os.path.join(path, MANIFEST_NAME)):
        return [path]
    if not os.path.isdir(path):
        return []
    return [os.path.join(path, d) for d in sorted(os.listdir(path))
            if os.path.isfile(os.path.join(path, d, MANIFEST_NAME))]


def discover_models(search_paths=None) -> list[ModelEntry]:
    """Every model found under `search_paths` (default: `default_search_paths()`).

    Each path may be a model folder (holds a manifest), a directory of model
    folders, or a bare schema-v2 database file / directory of them. Keys are
    made unique if two models would share one.
    """
    search_paths = default_search_paths() if search_paths is None else search_paths
    found: list[ModelEntry] = []
    for path in search_paths:
        manifest_dirs = _candidate_directories(path)
        for directory in manifest_dirs:
            found.append(_entry_from_manifest(directory))
        if os.path.isfile(path):
            found += _entries_from_bare_database(os.path.abspath(path))
        elif os.path.isdir(path) and not manifest_dirs:
            # bare databases: directly in the directory, or one level down
            dbs = [os.path.join(path, f) for f in sorted(os.listdir(path)) if f.endswith(".db")]
            for sub in sorted(os.listdir(path)):
                full = os.path.join(path, sub)
                if os.path.isdir(full) and not os.path.isfile(os.path.join(full, MANIFEST_NAME)):
                    dbs += [os.path.join(full, f) for f in sorted(os.listdir(full)) if f.endswith(".db")]
            for db in dbs:
                found += _entries_from_bare_database(os.path.abspath(db))

    unique: list[ModelEntry] = []
    seen: set[str] = set()
    for entry in found:
        key = entry.key
        if key in seen:
            stem = os.path.basename(os.path.normpath(entry.directory or os.path.dirname(entry.db_path) or "model"))
            key = _slug(f"{key}-{stem}", None)
            n = 2
            while key in seen:
                key, n = f"{key}-{n}", n + 1
            entry = ModelEntry(**{**entry.__dict__, "key": key})
        seen.add(key)
        unique.append(entry)
    return unique


def resolve_model(spec: str | None = None, search_paths=None) -> ModelEntry:
    """The usable `ModelEntry` for `spec`: a key (``Human-GEM-2.0.1``), a model
    folder, or a database file; None means `DEFAULT_MODEL_KEY`.

    Raises `KeyError` (listing what exists) for an unknown key, and
    `ValueError` (with the reasons) for a model that was found but cannot be used.
    """
    spec = spec or DEFAULT_MODEL_KEY
    if os.path.exists(spec):
        entries = discover_models([spec])
        if not entries:
            raise KeyError(f"{spec} is neither a model folder (manifest.json) nor a schema-v2 database")
        if len(entries) > 1:
            raise KeyError(f"{spec} holds several models ({[e.key for e in entries]}); pass one of their keys")
        entry = entries[0]
    else:
        entries = discover_models(search_paths)
        entry = next((e for e in entries if e.key == spec), None)
        if entry is None:
            raise KeyError(f"unknown model {spec!r}; available: {[e.key for e in entries]}")
    if not entry.available:
        reasons = "; ".join(entry.problems) or "its model file was not found next to the database"
        raise ValueError(f"model {entry.key!r} cannot be used: {reasons}")
    return entry
