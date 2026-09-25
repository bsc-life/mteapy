import pandas as pd
import pytest

from mteapy.context_scoring import (
    RouteScore,
    build_complex_cache,
    score_reaction,
    score_task,
    score_tasks_matrix,
    tied_argmax,
)


@pytest.mark.parametrize("scores,expected", [
    ({}, ()),
    ({"A": 1.0}, ("A",)),
    ({"A": 1.0, "B": 2.0}, ("B",)),
    ({"A": 2.0, "B": 2.0, "C": 1.0}, ("A", "B")),
    ({"A": 1.0, "B": 1.0 + 1e-12}, ("A", "B")),  # within floating tolerance -> tied
])
def test_tied_argmax(scores, expected):
    assert tied_argmax(scores) == expected


def test_score_reaction_and_within_complex_or_across_complexes():
    # (A and B) or (C) -- AND takes the min within a complex, OR takes the
    # max across complexes; C alone should win here since it's unopposed.
    complexes = (("A", "B"), ("C",))
    gene_dict = {"A": 5.0, "B": 1.0, "C": 3.0}
    score, winner = score_reaction(complexes, gene_dict)
    assert score == 3.0
    assert winner == ("C",)


def test_score_reaction_missing_gene_defaults_to_zero():
    complexes = (("A", "B"),)
    score, winner = score_reaction(complexes, {"A": 5.0})  # B absent from gene_dict
    assert score == 0.0
    assert winner == ("A", "B")


def test_score_reaction_no_complexes_is_zero_with_no_winner():
    assert score_reaction((), {"A": 5.0}) == (0.0, None)


def test_score_reaction_ambiguous_winner_is_none():
    complexes = (("A",), ("B",))
    score, winner = score_reaction(complexes, {"A": 2.0, "B": 2.0})
    assert score == 2.0
    assert winner is None  # the score is well-defined, which complex earned it is not


def test_score_task_single_route_is_never_tied():
    complex_cache = {"R1": (("g1",),), "R2": (("g2",),)}
    task_routes = {1: frozenset({"R1", "R2"})}
    gene_dict = {"g1": 3.0, "g2": 1.0}

    result = score_task(task_routes, complex_cache, gene_dict, task_id="T")
    assert not result.is_tied
    assert result.winning_route_ids == (1,)
    assert result.score == 1.0  # min(3.0, 1.0)


def test_score_task_detects_genuine_route_tie():
    complex_cache = {"R1": (("g1",),), "R2": (("g2",),)}
    # Two single-reaction routes that happen to score identically.
    task_routes = {1: frozenset({"R1"}), 2: frozenset({"R2"})}
    gene_dict = {"g1": 5.0, "g2": 5.0}

    result = score_task(task_routes, complex_cache, gene_dict, task_id="T")
    assert result.is_tied
    assert result.winning_route_ids == (1, 2)
    assert result.score == 5.0


def test_score_task_min_vs_median_can_change_the_winner():
    complex_cache = {
        "R1": (("g1",),), "R2": (("g2",),), "R3": (("g3",),),
    }
    # Route A: two reactions, one very low (bottleneck) one very high.
    # Route B: two reactions, both middling.
    task_routes = {
        "A": frozenset({"R1", "R2"}),
        "B": frozenset({"R3"}),
    }
    gene_dict = {"g1": 0.1, "g2": 100.0, "g3": 5.0}

    min_result = score_task(task_routes, complex_cache, gene_dict, aggregation="min")
    assert min_result.winning_route_ids == ("B",)  # A's min (0.1) loses to B's 5.0

    median_result = score_task(task_routes, complex_cache, gene_dict, aggregation="median")
    # median of a 1-reaction route is just that reaction's score (5.0);
    # median of A's two reactions is the mean of 0.1 and 100 -> 50.05
    assert median_result.winning_route_ids == ("A",)


def test_score_task_unsupported_aggregation_raises():
    complex_cache = {"R1": (("g1",),)}
    with pytest.raises(ValueError):
        score_task({1: frozenset({"R1"})}, complex_cache, {"g1": 1.0}, aggregation="mean")


def test_build_complex_cache_and_score_tasks_matrix(toy_model):
    reaction_ids = ["R1", "R2", "R3", "R4"]
    cache = build_complex_cache(toy_model, reaction_ids)
    assert cache["R1"] == (("g1",),)
    assert cache["R3"] == (("g3",),)

    tasks_routes = {
        "LINEAR": {1: frozenset({"R1", "R2"}), 2: frozenset({"R3"})},
    }
    expr_df = pd.DataFrame({
        "sample_A": {"g1": 10.0, "g2": 10.0, "g3": 1.0},   # R1+R2 route wins
        "sample_B": {"g1": 1.0, "g2": 1.0, "g3": 10.0},    # R3 route wins
    })

    matrix = score_tasks_matrix(tasks_routes, toy_model, expr_df)
    assert matrix.shape == (1, 2)
    assert matrix.loc["LINEAR", "sample_A"] == 10.0  # min(10, 10) beats R3's 1.0
    assert matrix.loc["LINEAR", "sample_B"] == 10.0  # R3 alone scores 10.0
