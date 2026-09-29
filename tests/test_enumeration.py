import pytest

from mteapy.enumeration import compute_task_alternate_routes, enumerate_alternate_routes
from mteapy.tasks import BoundedMetabolite, MetabolicTask


def _reaction_sets(routes):
    return [r.reactions for r in routes]

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
    assert _reaction_sets(routes) == [frozenset({"R3"}), frozenset({"R1", "R2"})]


def test_enumerate_captures_each_route_own_flux(toy_model):
    routes = enumerate_alternate_routes(toy_model, LINEAR, max_routes=10)

    # R3 alone: A[c] -> C[c] direct, task fixes 1 unit in and out.
    r3_route = next(r for r in routes if r.reactions == frozenset({"R3"}))
    assert r3_route.fluxes == {"R3": pytest.approx(1.0)}

    # R1+R2: same 1 unit has to move through both steps in series.
    r1r2_route = next(r for r in routes if r.reactions == frozenset({"R1", "R2"}))
    assert r1r2_route.fluxes["R1"] == pytest.approx(1.0)
    assert r1r2_route.fluxes["R2"] == pytest.approx(1.0)


def test_enumerate_respects_max_routes(toy_model):
    routes = enumerate_alternate_routes(toy_model, LINEAR, max_routes=1)
    assert _reaction_sets(routes) == [frozenset({"R3"})]


@pytest.mark.parametrize("max_routes", [0, -1])
def test_enumerate_rejects_non_positive_max_routes(toy_model, max_routes):
    # The reference-route solve happens unconditionally before the
    # "additional routes" loop, so without this guard max_routes=0 would
    # still silently yield 1 route instead of 0.
    with pytest.raises(ValueError, match="max_routes"):
        enumerate_alternate_routes(toy_model, LINEAR, max_routes=max_routes)


def test_enumerate_no_alt_single_route(toy_model):
    routes = enumerate_alternate_routes(toy_model, NO_ALT, max_routes=10)
    assert _reaction_sets(routes) == [frozenset({"R4"})]


def test_enumerate_infeasible_task_returns_empty(toy_model):
    assert enumerate_alternate_routes(toy_model, UNREACHABLE, max_routes=10) == []


def test_seed_routes_finds_only_the_remaining_alternative(toy_model):
    # Seed with the reference route already "known" -- enumeration should
    # skip straight past it and find only the one genuinely new route.
    routes = enumerate_alternate_routes(toy_model, LINEAR, max_routes=10, seed_routes=[frozenset({"R3"})])
    assert _reaction_sets(routes) == [frozenset({"R1", "R2"})]


def test_seed_routes_with_everything_already_known_finds_nothing_new(toy_model):
    # Both of LINEAR's routes seeded -- a resumed enumeration correctly
    # reports "nothing further exists" rather than rediscovering either.
    seeds = [frozenset({"R3"}), frozenset({"R1", "R2"})]
    routes = enumerate_alternate_routes(toy_model, LINEAR, max_routes=10, seed_routes=seeds)
    assert routes == []


def test_seed_routes_empty_list_behaves_like_no_seeding(toy_model):
    routes = enumerate_alternate_routes(toy_model, LINEAR, max_routes=10, seed_routes=[])
    assert _reaction_sets(routes) == [frozenset({"R3"}), frozenset({"R1", "R2"})]


def test_seed_routes_sharing_a_reaction_does_not_duplicate_constraints(toy_model):
    # R1+R2 and R3 share no reactions in this fixture, so exercise the
    # "reaction already has a y-var from an earlier seed" path directly by
    # seeding the same route reaction twice via two overlapping (fabricated)
    # seed sets -- this must not raise (e.g. from re-adding a variable/
    # constraint with a name that already exists).
    seeds = [frozenset({"R3"}), frozenset({"R3"})]  # deliberately duplicated
    routes = enumerate_alternate_routes(toy_model, LINEAR, max_routes=10, seed_routes=seeds)
    assert _reaction_sets(routes) == [frozenset({"R1", "R2"})]


def test_stop_info_reports_max_routes_reached(toy_model):
    stop_info = {}
    enumerate_alternate_routes(toy_model, LINEAR, max_routes=1, stop_info=stop_info)
    assert stop_info["reason"] == "max_routes_reached"


def test_stop_info_reports_infeasible_when_genuinely_exhausted(toy_model):
    stop_info = {}
    enumerate_alternate_routes(toy_model, LINEAR, max_routes=10, stop_info=stop_info)
    # LINEAR has exactly 2 routes; the 3rd cut makes the problem infeasible.
    assert stop_info["reason"] == "infeasible"


def test_stop_info_reports_infeasible_for_an_infeasible_task(toy_model):
    stop_info = {}
    enumerate_alternate_routes(toy_model, UNREACHABLE, max_routes=10, stop_info=stop_info)
    assert stop_info["reason"] == "infeasible"


def test_degenerate_duplicate_support_stops_without_reyielding(toy_model, monkeypatch):
    """Regression test for a real bug found on large Human-GEM tasks: the
    cut constraints are enforced over the solver's own binary y variables,
    but the *reported* support is independently re-derived from each
    solve's raw flux magnitudes -- for a large enough/degenerate problem,
    numerical noise can let two solves the cuts consider distinct still
    threshold out to the exact same reaction set, which must never be
    silently recorded as a second genuine route. Reproduced here by forcing
    `_support` to report the same set twice in a row, independent of any
    real solver's numerical behavior.
    """
    import mteapy.enumeration as enum_mod

    call_count = {"n": 0}

    def fake_support(solution, flux_threshold):
        call_count["n"] += 1
        return frozenset({"R3"})  # always "finds" the same support

    monkeypatch.setattr(enum_mod, "_support", fake_support)

    stop_info = {}
    routes = enumerate_alternate_routes(toy_model, LINEAR, max_routes=10, stop_info=stop_info)

    assert _reaction_sets(routes) == [frozenset({"R3"})]  # only yielded once
    assert stop_info["reason"] == "degenerate_duplicate"
    assert call_count["n"] == 2  # solved once more to detect the repeat, then stopped


def test_seed_matching_the_first_new_solve_is_also_caught_as_degenerate(toy_model, monkeypatch):
    """The same duplicate-detection must also apply against `seed_routes`,
    not just routes found earlier in the current call -- a safety net in
    case the numerical degeneracy strikes on the very first post-seed
    solve during a resume."""
    import mteapy.enumeration as enum_mod

    monkeypatch.setattr(enum_mod, "_support", lambda solution, flux_threshold: frozenset({"R1", "R2"}))

    stop_info = {}
    routes = enumerate_alternate_routes(
        toy_model, LINEAR, max_routes=10, seed_routes=[frozenset({"R1", "R2"})], stop_info=stop_info,
    )
    assert routes == []
    assert stop_info["reason"] == "degenerate_duplicate"


def test_compute_task_alternate_routes_batch(toy_model):
    routes_by_task, summary_df = compute_task_alternate_routes(
        toy_model, [LINEAR, NO_ALT, UNREACHABLE], verbose=False
    )
    assert summary_df.loc["LINEAR", "n_routes"] == 2
    assert summary_df.loc["NO_ALT", "n_routes"] == 1
    assert not summary_df.loc["UNREACHABLE", "included"]
    assert "UNREACHABLE" not in routes_by_task
    assert routes_by_task["LINEAR"][0].reactions == frozenset({"R3"})
