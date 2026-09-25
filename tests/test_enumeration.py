from mteapy.enumeration import compute_task_alternate_routes, enumerate_alternate_routes
from mteapy.tasks import BoundedMetabolite, MetabolicTask

# Same LINEAR task as in test_task_model.py: A[c] -> C[c] has two redundant
# routes (R1+R2, or the direct R3 shortcut). Any convex combination of the
# two is also optimal for total real-reaction flux only at the t=0 extreme
# (R3 alone costs 1, R1+R2 costs 2t+... -- minimizing total flux strictly
# prefers R3 alone as the reference route; see conftest.py docstring).
LINEAR = MetabolicTask(
    id="LINEAR", description="A to C (two redundant routes)",
    inputs=[BoundedMetabolite("A[c]", 1, 1)],
    outputs=[BoundedMetabolite("C[c]", 1, 1)],
)

# A to D has a single possible route (R4) -- no alternate optima.
NO_ALT = MetabolicTask(
    id="NO_ALT", description="A to D (single route via R4)",
    inputs=[BoundedMetabolite("A[c]", 1, 1)],
    outputs=[BoundedMetabolite("D[c]", 1, 1)],
)

UNREACHABLE = MetabolicTask(
    id="UNREACHABLE", description="D to A (wrong direction, infeasible)",
    inputs=[BoundedMetabolite("D[c]", 1, 1)],
    outputs=[BoundedMetabolite("A[c]", 1, 1)],
)


def test_enumerate_finds_both_routes_sparsest_first(toy_model):
    routes = enumerate_alternate_routes(toy_model, LINEAR, max_routes=10)

    # Reference route minimizes total real flux -> R3 alone (cost 1) beats
    # R1+R2 (cost 2). Once R3 is cut out, R1+R2 is the only way left to
    # satisfy the task, then nothing remains.
    assert routes == [frozenset({"R3"}), frozenset({"R1", "R2"})]


def test_enumerate_respects_max_routes(toy_model):
    routes = enumerate_alternate_routes(toy_model, LINEAR, max_routes=1)
    assert routes == [frozenset({"R3"})]


def test_enumerate_no_alt_single_route(toy_model):
    routes = enumerate_alternate_routes(toy_model, NO_ALT, max_routes=10)
    assert routes == [frozenset({"R4"})]


def test_enumerate_infeasible_task_returns_empty(toy_model):
    assert enumerate_alternate_routes(toy_model, UNREACHABLE, max_routes=10) == []


def test_compute_task_alternate_routes_batch(toy_model):
    routes_by_task, summary_df = compute_task_alternate_routes(
        toy_model, [LINEAR, NO_ALT, UNREACHABLE], verbose=False
    )
    assert summary_df.loc["LINEAR", "n_routes"] == 2
    assert summary_df.loc["NO_ALT", "n_routes"] == 1
    assert not summary_df.loc["UNREACHABLE", "included"]
    assert "UNREACHABLE" not in routes_by_task
    assert routes_by_task["LINEAR"][0] == frozenset({"R3"})
