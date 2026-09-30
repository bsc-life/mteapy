from mteapy.routes import (
    connect,
    get_enumeration_status,
    get_sources_for_task_list,
    get_task_definition_hash,
    get_task_source,
    list_tasks,
    load_multiroute_tasks,
    load_route_fluxes,
    load_task_routes,
    reaction_set_hash,
    record_enumeration_result,
    register_model,
    register_task,
    register_task_source,
    reset_source,
    save_route_fluxes,
)


def _fresh_db():
    return connect(":memory:")


def test_register_model_is_idempotent_by_content_hash():
    conn = _fresh_db()
    id1 = register_model(conn, "model.xml", sha256="abc123", n_reactions=10, n_genes=5)
    id2 = register_model(conn, "model.xml", sha256="abc123", n_reactions=10, n_genes=5)
    assert id1 == id2

    other = register_model(conn, "model.xml", sha256="def456", n_reactions=10, n_genes=5)
    assert other != id1


def test_register_task_source_records_provenance_and_reports_no_change_on_first_call():
    conn = _fresh_db()
    changed = register_task_source(
        conn, "cellfie_consensus", "Human-GEM/data/metabolicTasks/metabolicTasks_CellfieConsensus.txt",
        sha256="abc123", origin_repo="bsc-life/Human-GEM", origin_ref="4f25a6a",
    )
    assert changed is False

    record = get_task_source(conn, "cellfie_consensus")
    assert record["source"] == "cellfie_consensus"
    assert record["sha256"] == "abc123"
    assert record["origin_repo"] == "bsc-life/Human-GEM"
    assert record["origin_ref"] == "4f25a6a"


def test_register_task_source_reports_change_when_hash_differs():
    conn = _fresh_db()
    register_task_source(conn, "full", "tasks.txt", sha256="original_hash")
    changed = register_task_source(conn, "full", "tasks.txt", sha256="different_hash")
    assert changed is True
    assert get_task_source(conn, "full")["sha256"] == "different_hash"


def test_register_task_source_reregistering_same_hash_reports_no_change():
    conn = _fresh_db()
    register_task_source(conn, "full", "tasks.txt", sha256="same_hash")
    changed = register_task_source(conn, "full", "tasks.txt", sha256="same_hash")
    assert changed is False


def test_get_task_source_returns_none_when_never_registered():
    conn = _fresh_db()
    assert get_task_source(conn, "unknown_source") is None


def test_register_task_source_groups_multiple_sources_under_one_task_list():
    conn = _fresh_db()
    register_task_source(conn, "cellfie_consensus", "tasks.txt", sha256="abc", task_list="cellfie")
    register_task_source(conn, "cellfie_consensus_gurobi", "tasks.txt", sha256="abc", task_list="cellfie")
    register_task_source(conn, "full", "full_tasks.txt", sha256="def", task_list="full")

    assert get_sources_for_task_list(conn, "cellfie") == ["cellfie_consensus", "cellfie_consensus_gurobi"]
    assert get_sources_for_task_list(conn, "full") == ["full"]
    assert get_sources_for_task_list(conn, "unknown") == []


def test_register_task_source_omitting_task_list_does_not_clear_a_previous_value():
    conn = _fresh_db()
    register_task_source(conn, "cellfie_consensus", "tasks.txt", sha256="abc", task_list="cellfie")
    # A later call (e.g. a re-enumeration run) that doesn't know/pass
    # task_list must not silently erase the earlier classification.
    register_task_source(conn, "cellfie_consensus", "tasks.txt", sha256="abc")
    assert get_task_source(conn, "cellfie_consensus")["task_list"] == "cellfie"


def test_register_task_source_omitting_origin_does_not_clear_a_previous_value():
    conn = _fresh_db()
    register_task_source(conn, "full", "tasks.txt", sha256="abc",
                          origin_repo="bsc-life/Human-GEM", origin_ref="4f25a6a")
    # A caller that couldn't determine origin_repo/origin_ref for this call
    # (e.g. _git_provenance's own best-effort failure mode) must not wipe
    # out an already-recorded origin.
    register_task_source(conn, "full", "tasks.txt", sha256="abc")
    record = get_task_source(conn, "full")
    assert record["origin_repo"] == "bsc-life/Human-GEM"
    assert record["origin_ref"] == "4f25a6a"


def test_record_and_load_single_task_routes():
    conn = _fresh_db()
    model_id = register_model(conn, "model.xml", sha256="abc123")
    register_task(conn, "full", "1", "test task")

    routes = [frozenset({"R1", "R2"}), frozenset({"R3"})]
    route_ids = record_enumeration_result(
        conn, "full", "1", model_id, routes,
        status="optimal", max_routes=10, hit_cap=False, truncated=False, elapsed_seconds=1.5,
    )
    assert len(route_ids) == 2
    assert len(set(route_ids)) == 2  # distinct routes get distinct ids

    loaded = load_task_routes(conn, "full", "1", model_id)
    assert loaded[route_ids[0]] == frozenset({"R1", "R2"})
    assert loaded[route_ids[1]] == frozenset({"R3"})


def test_record_deduplicates_identical_route_sets():
    conn = _fresh_db()
    model_id = register_model(conn, "model.xml", sha256="abc123")
    register_task(conn, "full", "1", "test task")

    first_ids = record_enumeration_result(
        conn, "full", "1", model_id, [frozenset({"R1", "R2"})],
        status="optimal", max_routes=10, hit_cap=False, truncated=False, elapsed_seconds=1.0,
    )
    # Re-running enumeration with the same result (e.g. a rerun after a
    # model/task update with no actual change) must not create a duplicate
    # route row -- it should resolve back to the same route_id.
    second_ids = record_enumeration_result(
        conn, "full", "1", model_id, [frozenset({"R2", "R1"})],  # same set, different order
        status="optimal", max_routes=10, hit_cap=False, truncated=False, elapsed_seconds=1.0,
    )
    assert first_ids == second_ids

    routes = load_task_routes(conn, "full", "1", model_id)
    assert len(routes) == 1


def test_load_multiroute_tasks_filters_single_route_tasks():
    conn = _fresh_db()
    model_id = register_model(conn, "model.xml", sha256="abc123")
    register_task(conn, "full", "single", "one route only")
    register_task(conn, "full", "multi", "two routes")

    record_enumeration_result(
        conn, "full", "single", model_id, [frozenset({"R1"})],
        status="optimal", max_routes=10, hit_cap=False, truncated=False, elapsed_seconds=1.0,
    )
    record_enumeration_result(
        conn, "full", "multi", model_id, [frozenset({"R1"}), frozenset({"R2"})],
        status="optimal", max_routes=10, hit_cap=False, truncated=False, elapsed_seconds=1.0,
    )

    result = load_multiroute_tasks(conn, model_id)
    assert ("full", "multi") in result
    assert ("full", "single") not in result
    assert len(result[("full", "multi")]) == 2


def test_reaction_set_hash_is_order_independent():
    assert reaction_set_hash(["R1", "R2"]) == reaction_set_hash(["R2", "R1"])
    assert reaction_set_hash(["R1", "R2"]) != reaction_set_hash(["R1", "R3"])


def test_list_tasks_returns_description_and_route_count_ordered_numerically():
    conn = _fresh_db()
    model_id = register_model(conn, "model.xml", sha256="abc123")
    # Registered out of numeric order and with task "10" before "2" to check
    # the listing sorts by numeric task id, not lexicographic string order.
    register_task(conn, "full", "10", "tenth task")
    register_task(conn, "full", "2", "second task")
    record_enumeration_result(
        conn, "full", "10", model_id, [frozenset({"R1"})],
        status="optimal", max_routes=10, hit_cap=False, truncated=False, elapsed_seconds=1.0,
    )
    record_enumeration_result(
        conn, "full", "2", model_id, [frozenset({"R1"}), frozenset({"R2"})],
        status="optimal", max_routes=10, hit_cap=False, truncated=False, elapsed_seconds=1.0,
    )

    result = list_tasks(conn, model_id)
    assert [r["task_id"] for r in result] == ["2", "10"]
    assert result[0]["description"] == "second task"
    assert result[0]["n_routes"] == 2
    assert result[1]["n_routes"] == 1


def test_list_tasks_sorts_non_numeric_ids_after_numeric_ones():
    """A plain CAST(task_id AS INTEGER) would collapse every non-numeric id
    to 0, clustering "ER" together with whichever numeric ids happen to
    also collapse to 0 -- non-numeric ids must sort after all numeric ones,
    not get silently interleaved at the front."""
    conn = _fresh_db()
    model_id = register_model(conn, "model.xml", sha256="abc123")
    for task_id in ["ER", "10", "2", "AA"]:
        record_enumeration_result(
            conn, "full", task_id, model_id, [frozenset({"R1"})],
            status="optimal", max_routes=10, hit_cap=False, truncated=False, elapsed_seconds=1.0,
        )

    result = list_tasks(conn, model_id)
    assert [r["task_id"] for r in result] == ["2", "10", "AA", "ER"]


def test_list_tasks_blank_description_for_unregistered_task():
    conn = _fresh_db()
    model_id = register_model(conn, "model.xml", sha256="abc123")
    # A route recorded without ever calling register_task for it -- the
    # LEFT JOIN must not drop the row, just leave description blank.
    record_enumeration_result(
        conn, "full", "1", model_id, [frozenset({"R1"})],
        status="optimal", max_routes=10, hit_cap=False, truncated=False, elapsed_seconds=1.0,
    )
    result = list_tasks(conn, model_id)
    assert result == [{"source": "full", "task_id": "1", "description": "", "n_routes": 1}]


def test_record_enumeration_result_n_routes_reflects_distinct_count_not_input_length():
    """A caller-supplied `routes` list containing a duplicate reaction-set
    (both entries resolve to the same route_id) must not inflate n_routes
    to the raw input length -- it should report the true distinct count."""
    conn = _fresh_db()
    model_id = register_model(conn, "model.xml", sha256="abc123")
    route_ids = record_enumeration_result(
        conn, "full", "1", model_id,
        [frozenset({"R1"}), frozenset({"R1"}), frozenset({"R2"})],  # R1 listed twice
        status="optimal", max_routes=10, hit_cap=False, truncated=False, elapsed_seconds=1.0,
    )
    assert len(set(route_ids)) == 2  # R1 and R2, R1's duplicate resolves to the same route_id
    status = get_enumeration_status(conn, "full", "1", model_id)
    assert status["n_routes"] == 2


def test_record_enumeration_result_stores_flux_for_new_routes():
    conn = _fresh_db()
    model_id = register_model(conn, "model.xml", sha256="abc123")
    route_ids = record_enumeration_result(
        conn, "full", "1", model_id, [frozenset({"R1", "R2"})],
        status="optimal", max_routes=10, hit_cap=False, truncated=False, elapsed_seconds=1.0,
        fluxes=[{"R1": 1.0, "R2": -0.5}],
    )
    assert load_route_fluxes(conn, route_ids[0]) == {"R1": 1.0, "R2": -0.5}


def test_record_enumeration_result_without_fluxes_leaves_flux_null():
    conn = _fresh_db()
    model_id = register_model(conn, "model.xml", sha256="abc123")
    route_ids = record_enumeration_result(
        conn, "full", "1", model_id, [frozenset({"R1"})],
        status="optimal", max_routes=10, hit_cap=False, truncated=False, elapsed_seconds=1.0,
    )
    assert load_route_fluxes(conn, route_ids[0]) == {}


def test_record_enumeration_result_backfills_flux_for_already_existing_route():
    conn = _fresh_db()
    model_id = register_model(conn, "model.xml", sha256="abc123")
    # First recorded with no flux (as an older run, before this feature, would have).
    first_ids = record_enumeration_result(
        conn, "full", "1", model_id, [frozenset({"R1"})],
        status="optimal", max_routes=10, hit_cap=False, truncated=False, elapsed_seconds=1.0,
    )
    assert load_route_fluxes(conn, first_ids[0]) == {}

    # Re-running enumeration now supplies flux for that same (already
    # existing, same reaction-set-hash) route -- it should get backfilled,
    # not silently dropped because the route row itself wasn't new.
    second_ids = record_enumeration_result(
        conn, "full", "1", model_id, [frozenset({"R1"})],
        status="optimal", max_routes=10, hit_cap=False, truncated=False, elapsed_seconds=1.0,
        fluxes=[{"R1": 2.5}],
    )
    assert second_ids == first_ids
    assert load_route_fluxes(conn, first_ids[0]) == {"R1": 2.5}


def test_record_enumeration_result_partial_flux_backfill_does_not_null_out_existing_values():
    """A caller-supplied partial fluxes dict for an already-existing route
    must only touch the reactions it actually mentions -- never null out a
    previously-stored value for a reaction it simply didn't include."""
    conn = _fresh_db()
    model_id = register_model(conn, "model.xml", sha256="abc123")
    route_ids = record_enumeration_result(
        conn, "full", "1", model_id, [frozenset({"R1", "R2"})],
        status="optimal", max_routes=10, hit_cap=False, truncated=False, elapsed_seconds=1.0,
        fluxes=[{"R1": 1.0, "R2": 2.0}],
    )
    assert load_route_fluxes(conn, route_ids[0]) == {"R1": 1.0, "R2": 2.0}

    # Re-recording the same route with a partial flux dict (only R1) must
    # leave R2's already-stored flux alone, not null it out.
    record_enumeration_result(
        conn, "full", "1", model_id, [frozenset({"R1", "R2"})],
        status="optimal", max_routes=10, hit_cap=False, truncated=False, elapsed_seconds=1.0,
        fluxes=[{"R1": 9.0}],
    )
    assert load_route_fluxes(conn, route_ids[0]) == {"R1": 9.0, "R2": 2.0}


def test_save_route_fluxes_persists_for_later_load():
    conn = _fresh_db()
    model_id = register_model(conn, "model.xml", sha256="abc123")
    route_ids = record_enumeration_result(
        conn, "full", "1", model_id, [frozenset({"R1", "R2"})],
        status="optimal", max_routes=10, hit_cap=False, truncated=False, elapsed_seconds=1.0,
    )
    assert load_route_fluxes(conn, route_ids[0]) == {}

    save_route_fluxes(conn, route_ids[0], {"R1": 3.0, "R2": -1.0})
    assert load_route_fluxes(conn, route_ids[0]) == {"R1": 3.0, "R2": -1.0}


def test_flux_column_migration_on_pre_existing_database(tmp_path):
    db_path = str(tmp_path / "legacy.db")
    # Simulate a database created before the flux column existed: build it
    # from the schema string with the ALTER TABLE-based migration skipped.
    import sqlite3

    from mteapy.routes import SCHEMA

    legacy_conn = sqlite3.connect(db_path)
    legacy_conn.executescript(SCHEMA)
    legacy_conn.close()

    # connect() must transparently add the missing column, not error.
    conn = connect(db_path)
    model_id = register_model(conn, "model.xml", sha256="abc123")
    route_ids = record_enumeration_result(
        conn, "full", "1", model_id, [frozenset({"R1"})],
        status="optimal", max_routes=10, hit_cap=False, truncated=False, elapsed_seconds=1.0,
        fluxes=[{"R1": 1.0}],
    )
    assert load_route_fluxes(conn, route_ids[0]) == {"R1": 1.0}


def test_register_task_stores_and_preserves_definition_hash():
    conn = _fresh_db()
    register_task(conn, "full", "1", "a task", definition_hash="abc123")
    assert get_task_definition_hash(conn, "full", "1") == "abc123"

    # A later call that doesn't pass a hash (e.g. from code that doesn't
    # know/care about resume-safety) must not erase the one already stored.
    register_task(conn, "full", "1", "a task, renamed")
    assert get_task_definition_hash(conn, "full", "1") == "abc123"

    # An explicit new hash still overwrites, e.g. after the task genuinely changed.
    register_task(conn, "full", "1", "a task, renamed", definition_hash="def456")
    assert get_task_definition_hash(conn, "full", "1") == "def456"


def test_get_task_definition_hash_none_when_never_set():
    conn = _fresh_db()
    register_task(conn, "full", "1", "a task")
    assert get_task_definition_hash(conn, "full", "1") is None
    assert get_task_definition_hash(conn, "full", "nonexistent") is None


def test_get_enumeration_status_none_before_any_run():
    conn = _fresh_db()
    model_id = register_model(conn, "model.xml", sha256="abc123")
    assert get_enumeration_status(conn, "full", "1", model_id) is None


def test_get_enumeration_status_is_exhaustive_only_when_not_capped_or_truncated():
    conn = _fresh_db()
    model_id = register_model(conn, "model.xml", sha256="abc123")

    record_enumeration_result(
        conn, "full", "capped", model_id, [frozenset({"R1"})],
        status="optimal", max_routes=1, hit_cap=True, truncated=False, elapsed_seconds=1.0,
    )
    record_enumeration_result(
        conn, "full", "truncated", model_id, [frozenset({"R1"})],
        status="optimal", max_routes=10, hit_cap=False, truncated=True, elapsed_seconds=1.0,
    )
    record_enumeration_result(
        conn, "full", "done", model_id, [frozenset({"R1"})],
        status="optimal", max_routes=10, hit_cap=False, truncated=False, elapsed_seconds=1.0,
    )

    assert get_enumeration_status(conn, "full", "capped", model_id)["is_exhaustive"] is False
    assert get_enumeration_status(conn, "full", "truncated", model_id)["is_exhaustive"] is False
    assert get_enumeration_status(conn, "full", "done", model_id)["is_exhaustive"] is True


def test_reset_source_clears_routes_but_keeps_task_metadata():
    conn = _fresh_db()
    model_id = register_model(conn, "model.xml", sha256="abc123")
    register_task(conn, "full", "1", "a task", definition_hash="abc123")
    record_enumeration_result(
        conn, "full", "1", model_id, [frozenset({"R1"}), frozenset({"R2"})],
        status="optimal", max_routes=10, hit_cap=True, truncated=False, elapsed_seconds=1.0,
        fluxes=[{"R1": 1.0}, {"R2": 1.0}],
    )
    assert len(load_task_routes(conn, "full", "1", model_id)) == 2

    reset_source(conn, "full", model_id)

    assert load_task_routes(conn, "full", "1", model_id) == {}
    assert get_enumeration_status(conn, "full", "1", model_id) is None
    # Task metadata (description, definition_hash) survives the reset.
    assert get_task_definition_hash(conn, "full", "1") == "abc123"


def test_reset_source_does_not_touch_other_sources_or_models():
    conn = _fresh_db()
    model_id = register_model(conn, "model.xml", sha256="abc123")
    other_model_id = register_model(conn, "other.xml", sha256="def456")

    record_enumeration_result(
        conn, "full", "1", model_id, [frozenset({"R1"})],
        status="optimal", max_routes=10, hit_cap=False, truncated=False, elapsed_seconds=1.0,
    )
    record_enumeration_result(
        conn, "cellfie_consensus", "1", model_id, [frozenset({"R2"})],
        status="optimal", max_routes=10, hit_cap=False, truncated=False, elapsed_seconds=1.0,
    )
    record_enumeration_result(
        conn, "full", "1", other_model_id, [frozenset({"R3"})],
        status="optimal", max_routes=10, hit_cap=False, truncated=False, elapsed_seconds=1.0,
    )

    reset_source(conn, "full", model_id)

    assert load_task_routes(conn, "full", "1", model_id) == {}
    assert len(load_task_routes(conn, "cellfie_consensus", "1", model_id)) == 1
    assert len(load_task_routes(conn, "full", "1", other_model_id)) == 1
