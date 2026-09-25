from mteapy.routes import (
    connect,
    load_multiroute_tasks,
    load_task_routes,
    reaction_set_hash,
    record_enumeration_result,
    register_model,
    register_task,
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
