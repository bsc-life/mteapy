import sqlite3

import pytest

from mteapy import taskdb
from mteapy.tasks import task_definition_hash

HEADER = "\tID\tDESCRIPTION\tSHOULD FAIL\tIN\tIN LB\tIN UB\tOUT\tOUT LB\tOUT UB\tEQU\tEQU LB\tEQU UB\tCHANGED RXN\tCHANGED LB\tCHANGED UB\tPRINT FLUX\tCOMMENTS\n"


def _row(task_id, description="", inp="", in_lb="", in_ub="", out="", out_lb="", out_ub="", equ="", equ_lb="",
         equ_ub="", ch="", ch_lb="", ch_ub="", comments=""):
    return "\t".join(["", task_id, description, "", inp, in_lb, in_ub, out, out_lb, out_ub, equ, equ_lb, equ_ub,
                      ch, ch_lb, ch_ub, "", comments]) + "\n"


@pytest.fixture
def db(toy_model, tmp_path):
    conn = taskdb.create_db(str(tmp_path / "t.db"))
    model_id = taskdb.register_model(conn, toy_model, "toy", "sha", "1")
    return conn, toy_model, model_id, tmp_path


def _import(db, rows, name="list"):
    conn, model, model_id, tmp_path = db
    path = tmp_path / f"{name}.txt"
    path.write_text(HEADER + "".join(rows))
    return taskdb.import_task_list(conn, model, model_id, name, str(path))


def test_model_entities_stored_once(db):
    conn, model, model_id, _ = db
    assert conn.execute("SELECT COUNT(*) FROM reactions").fetchone()[0] == len(model.reactions)
    assert conn.execute("SELECT COUNT(*) FROM metabolites").fetchone()[0] == len(model.metabolites)


def test_task_round_trips_with_identical_hash(db):
    conn = db[0]
    ids = _import(db, [
        _row("1", "A to C", inp="A[c]", in_lb="1", in_ub="1", out="C[c]", out_lb="1", out_ub="1"),
        _row("2", "everything", inp="A[c];X[e]", in_lb="1;0", in_ub="1;1000", out="C[c]", out_lb="0", out_ub="1000",
             equ="A[c] => B[c]", equ_lb="0", equ_ub="5", ch="R3", ch_lb="0", ch_ub="0"),
    ])
    from mteapy.tasks import parse_task_file
    parsed = {t.id: t for t in parse_task_file(str(db[3] / "list.txt"))}
    for task_id, pk in ids.items():
        assert task_definition_hash(taskdb.load_task(conn, pk)) == task_definition_hash(parsed[task_id])
    assert conn.execute("SELECT COUNT(*) FROM task_equation_terms").fetchone()[0] == 2


def test_unresolvable_metabolite_is_flagged_not_dropped(db):
    conn = db[0]
    ids = _import(db, [_row("1", inp="NOPE[c]", in_lb="1", in_ub="1", out="C[c]")])
    valid, note = conn.execute("SELECT valid, import_note FROM tasks WHERE task_pk = ?", (ids["1"],)).fetchone()
    assert valid == 0 and "NOPE[c]" in note
    with pytest.raises(ValueError):
        taskdb.load_task(conn, ids["1"])


def test_unknown_changed_reaction_is_flagged(db):
    conn = db[0]
    ids = _import(db, [_row("1", inp="A[c]", out="C[c]", ch="R_MISSING", ch_lb="0", ch_ub="0")])
    assert conn.execute("SELECT valid FROM tasks WHERE task_pk = ?", (ids["1"],)).fetchone()[0] == 0


def test_invalid_marker_in_comments_flags_task(db):
    conn = db[0]
    ids = _import(db, [_row("1", inp="A[c]", out="C[c]", comments="INVALID: needs curation")])
    valid, note = conn.execute("SELECT valid, import_note FROM tasks WHERE task_pk = ?", (ids["1"],)).fetchone()
    assert valid == 0 and "needs curation" in note


def test_duplicate_input_metabolite_is_flagged(db):
    conn = db[0]
    ids = _import(db, [_row("1", inp="A[c];A[c]", out="C[c]")])
    assert conn.execute("SELECT valid FROM tasks WHERE task_pk = ?", (ids["1"],)).fetchone()[0] == 0


def test_task_ids_unique_within_a_list_but_reusable_across_lists(db):
    _import(db, [_row("1", inp="A[c]", out="C[c]")], name="one")
    _import(db, [_row("1", inp="A[c]", out="D[c]")], name="two")
    with pytest.raises(sqlite3.IntegrityError):
        _import(db, [_row("1", inp="A[c]", out="C[c]")], name="one")


def test_route_flux_is_not_null(db):
    conn = db[0]
    ids = _import(db, [_row("1", inp="A[c]", out="C[c]")])
    rid = conn.execute(
        "INSERT INTO routes (task_pk, reaction_set_hash, n_reactions, first_seen_at) VALUES (?, 'h', 1, 'now')",
        (ids["1"],)).lastrowid
    with pytest.raises(sqlite3.IntegrityError):
        conn.execute("INSERT INTO route_reactions (route_id, rxn_pk, flux) VALUES (?, 1, NULL)", (rid,))


def test_create_db_refuses_to_overwrite(db):
    with pytest.raises(FileExistsError):
        taskdb.create_db(str(db[3] / "t.db"))


# ---------------------------------------------------------------------------
# Querying and recording
# ---------------------------------------------------------------------------

def _two_task_list(db):
    return _import(db, [
        _row("10", "A to C", inp="A[c]", in_lb="1", in_ub="1", out="C[c]", out_lb="1", out_ub="1"),
        _row("2", "A to D", inp="A[c]", in_lb="1", in_ub="1", out="D[c]", out_lb="1", out_ub="1"),
        _row("x", "text id", inp="A[c]", out="B[c]"),
    ], name="L")


def test_list_tasks_sorts_numeric_ids_first_then_text(db):
    _two_task_list(db)
    assert [t["task_id"] for t in taskdb.list_tasks(db[0], "L")] == ["2", "10", "x"]


def test_get_task_list_names_the_alternatives_on_a_miss(db):
    _two_task_list(db)
    with pytest.raises(KeyError, match="available: \\['L'\\]"):
        taskdb.get_task_list(db[0], "nope")


def test_record_enumeration_result_stores_routes_and_fluxes(db):
    conn = db[0]
    ids = _two_task_list(db)
    pk = ids["10"]
    r1, r2 = frozenset({"R1", "R2"}), frozenset({"R3"})
    route_ids = taskdb.record_enumeration_result(
        conn, pk, [r1, r2], [{"R1": 1.0, "R2": 1.0}, {"R3": 2.0}], status="optimal", max_routes=10,
        hit_cap=False, truncated=False, elapsed_seconds=1.5, solver="glpk",
    )
    assert len(set(route_ids)) == 2
    assert set(taskdb.load_task_routes(conn, pk).values()) == {r1, r2}
    assert taskdb.load_route_fluxes(conn, route_ids[1]) == {"R3": 2.0}
    status = taskdb.get_enumeration_status(conn, pk)
    assert status["is_exhaustive"] and status["n_routes"] == 2 and status["solver"] == "glpk"
    assert taskdb.load_task_list_routes(conn, "L") == {"10": taskdb.load_task_routes(conn, pk)}


def test_record_enumeration_result_dedups_and_accumulates_time(db):
    conn = db[0]
    pk = _two_task_list(db)["10"]
    kw = dict(status="optimal", max_routes=10, hit_cap=False, truncated=False, solver="glpk")
    first = taskdb.record_enumeration_result(conn, pk, [frozenset({"R3"})], [{"R3": 1.0}], elapsed_seconds=1.0, **kw)
    again = taskdb.record_enumeration_result(conn, pk, [frozenset({"R3"})], [{"R3": 1.0}], elapsed_seconds=2.0, **kw)
    assert first == again
    assert conn.execute("SELECT COUNT(*) FROM routes").fetchone()[0] == 1
    assert taskdb.get_enumeration_status(conn, pk)["cumulative_time_seconds"] == 3.0


def test_record_enumeration_result_flags_repeated_supports_as_truncated(db):
    conn = db[0]
    pk = _two_task_list(db)["10"]
    r = frozenset({"R3"})
    with pytest.warns(RuntimeWarning, match="repeated route support"):
        taskdb.record_enumeration_result(
            conn, pk, [r, r], [{"R3": 1.0}, {"R3": 1.0}], status="optimal", max_routes=2, hit_cap=True,
            truncated=False, elapsed_seconds=0.0,
        )
    status = taskdb.get_enumeration_status(conn, pk)
    assert status["n_routes"] == 1 and status["truncated"] and not status["hit_cap"]
    assert not status["is_exhaustive"]


def test_record_enumeration_result_requires_flux_for_every_reaction(db):
    conn = db[0]
    pk = _two_task_list(db)["10"]
    kw = dict(status="optimal", max_routes=10, hit_cap=False, truncated=False, elapsed_seconds=0.0)
    with pytest.raises(ValueError, match="no flux"):
        taskdb.record_enumeration_result(conn, pk, [frozenset({"R1", "R2"})], [{"R1": 1.0}], **kw)
    with pytest.raises(ValueError, match="not in the model"):
        taskdb.record_enumeration_result(conn, pk, [frozenset({"R_NOPE"})], [{"R_NOPE": 1.0}], **kw)
    with pytest.raises(ValueError, match="same length"):
        taskdb.record_enumeration_result(conn, pk, [frozenset({"R3"})], [], **kw)
    assert conn.execute("SELECT COUNT(*) FROM routes").fetchone()[0] == 0


def test_capped_enumeration_is_not_exhaustive(db):
    conn = db[0]
    pk = _two_task_list(db)["10"]
    taskdb.record_enumeration_result(conn, pk, [frozenset({"R3"})], [{"R3": 1.0}], status="optimal", max_routes=1,
                                     hit_cap=True, truncated=False, elapsed_seconds=0.0)
    assert not taskdb.get_enumeration_status(conn, pk)["is_exhaustive"]


def test_reset_task_list_routes_keeps_tasks(db):
    conn = db[0]
    ids = _two_task_list(db)
    taskdb.record_enumeration_result(conn, ids["10"], [frozenset({"R3"})], [{"R3": 1.0}], status="optimal",
                                     max_routes=10, hit_cap=False, truncated=False, elapsed_seconds=0.0)
    taskdb.reset_task_list_routes(conn, "L")
    assert conn.execute("SELECT COUNT(*) FROM routes").fetchone()[0] == 0
    assert conn.execute("SELECT COUNT(*) FROM route_reactions").fetchone()[0] == 0
    assert conn.execute("SELECT COUNT(*) FROM tasks").fetchone()[0] == 3
    assert taskdb.get_enumeration_status(conn, ids["10"]) is None


def test_connect_rejects_legacy_and_missing_databases(tmp_path):
    legacy = tmp_path / "legacy.db"
    sqlite3.connect(legacy).execute("CREATE TABLE routes (route_id INTEGER)").connection.commit()
    with pytest.raises(ValueError, match="migrate_db"):
        taskdb.connect(str(legacy))
    with pytest.raises(FileNotFoundError):
        taskdb.connect(str(tmp_path / "missing.db"))


def test_verify_model_file_detects_a_different_model(db, tmp_path):
    conn, _, model_id, _ = db
    other = tmp_path / "other.xml"
    other.write_text("not the model")
    with pytest.raises(ValueError, match="does not match"):
        taskdb.verify_model_file(conn, model_id, str(other))


def test_enumerated_toy_routes_round_trip_through_the_database(db):
    """End to end: enumerate on the toy model, record, and read back the same routes with fluxes."""
    from mteapy.enumeration import enumerate_alternate_routes

    conn, model, _, _ = db
    ids = _import(db, [_row("1", "A to C", inp="A[c]", in_lb="1", in_ub="1", out="C[c]", out_lb="1", out_ub="1")])
    task = taskdb.load_task(conn, ids["1"])
    found = enumerate_alternate_routes(model, task, max_routes=5)
    assert len(found) == 2  # R3 alone, or R1+R2
    taskdb.record_enumeration_result(
        conn, ids["1"], [r.reactions for r in found], [r.fluxes for r in found], status="optimal", max_routes=5,
        hit_cap=False, truncated=False, elapsed_seconds=0.0,
    )
    stored = taskdb.load_task_routes(conn, ids["1"])
    assert sorted(map(sorted, stored.values())) == sorted(map(sorted, (r.reactions for r in found)))
    for route_id, reactions in stored.items():
        assert set(taskdb.load_route_fluxes(conn, route_id)) == set(reactions)
