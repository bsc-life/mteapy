import pytest
import pandas as pd
import numpy as np

from cobra.core.gene import GPR
from mteapy.context_scoring import build_complex_cache
from mteapy.tide import (
    calculate_TIDE_scores, calculate_TIDE_scores_context_aware, calculate_random_TIDE_scores_context_aware,
    compute_TIDE,
)


@pytest.mark.parametrize("gene_dict,expected", [
    (dict(zip(["G" + str(i+1) for i in range(10)], [0,0,0,0,0,0,0,0,0,0])), [0, 0]),
    (dict(zip(["G" + str(i+1) for i in range(10)], [-5,-5,5,5,5,5,5,-5,0,0])), [0, -2.5]),
    (dict(zip(["G" + str(i+1) for i in range(10)], [0,-1,2,-1.2,-6,0.4,-0.9,-3.2,0,2.3])), [-3.5, -1.6])
])
def test_TIDE_scores(gene_dict, expected):
    task_structure = pd.DataFrame({"X1": [1,1,0,0], "X2": [0,0,1,1]}, index=["R1","R2","R3","R4"]).astype(bool)
    gpr_dict = dict(zip(
        task_structure.index,
        [GPR.from_string("G1 and G2"), GPR.from_string("G3 or G4 or G5"), \
         GPR.from_string("(G6 or G7) and G8"), GPR.from_string("G9 and G10")]
    ))

    assert all(calculate_TIDE_scores(gene_dict, task_structure, gpr_dict, or_func="absmax") == expected)


def test_calculate_TIDE_scores_context_aware_picks_the_absmax_supported_route(toy_model):
    # R1+R2 (g1=-5, g2=-5) has mean score -5; R3 alone (g3=2) has mean score
    # 2. or_func="absmax" picks the route by largest magnitude, so {R1,R2}
    # (|-5| > |2|) should win, not the higher raw value.
    tasks_routes = {"AC": {1: frozenset({"R1", "R2"}), 2: frozenset({"R3"})}}
    complex_cache = build_complex_cache(toy_model, {"R1", "R2", "R3"})
    gene_dict = {"g1": -5, "g2": -5, "g3": 2}

    scores = calculate_TIDE_scores_context_aware(gene_dict, tasks_routes, complex_cache, or_func="absmax")

    assert list(scores) == [-5]


def test_compute_TIDE_context_aware_requires_routes_and_complex_cache(toy_model):
    expr_data = pd.DataFrame({"lfc": [-5, -5, 2]}, index=["g1", "g2", "g3"])
    with pytest.raises(ValueError):
        compute_TIDE(expr_data, "lfc", None, toy_model, or_func="absmax", mapping_strategy="context-aware")


def test_compute_TIDE_context_aware_scores_the_winning_route(toy_model):
    tasks_routes = {"AC": {1: frozenset({"R1", "R2"}), 2: frozenset({"R3"})}}
    complex_cache = build_complex_cache(toy_model, {"R1", "R2", "R3"})
    expr_data = pd.DataFrame({"lfc": [-5, -5, 2]}, index=["g1", "g2", "g3"])

    results = compute_TIDE(
        expr_data, "lfc", None, toy_model, or_func="absmax", n_permutations=5, n_cpus=1,
        mapping_strategy="context-aware", tasks_routes=tasks_routes, complex_cache=complex_cache,
    )

    assert list(results["task_id"]) == ["AC"]
    assert results.loc[0, "score"] == -5


def test_compute_TIDE_context_aware_rejects_unsupported_permutation_strategy(toy_model):
    tasks_routes = {"AC": {1: frozenset({"R1", "R2"}), 2: frozenset({"R3"})}}
    complex_cache = build_complex_cache(toy_model, {"R1", "R2", "R3"})
    expr_data = pd.DataFrame({"lfc": [-5, -5, 2]}, index=["g1", "g2", "g3"])

    with pytest.raises(ValueError):
        compute_TIDE(
            expr_data, "lfc", None, toy_model, or_func="absmax",
            mapping_strategy="context-aware", tasks_routes=tasks_routes, complex_cache=complex_cache,
            permutation_strategy="bogus",
        )


def test_compute_TIDE_fixed_route_permutation_only_uses_the_winning_routes_genes(toy_model):
    # Route {R1,R2} (g1=-5, g2=-5) wins under absmax over {R3} (g3=2).
    # "fixed-route" should build its null distribution only from {g1, g2}
    # -- i.e. identical to directly permuting just that single route.
    tasks_routes = {"AC": {1: frozenset({"R1", "R2"}), 2: frozenset({"R3"})}}
    complex_cache = build_complex_cache(toy_model, {"R1", "R2", "R3"})
    gene_dict = {"g1": -5, "g2": -5, "g3": 2}
    expr_data = pd.DataFrame({"lfc": [-5, -5, 2]}, index=["g1", "g2", "g3"])

    results = compute_TIDE(
        expr_data, "lfc", None, toy_model, or_func="absmax", n_permutations=5, n_cpus=1, random_seed=0,
        mapping_strategy="context-aware", tasks_routes=tasks_routes, complex_cache=complex_cache,
        permutation_strategy="fixed-route", random_scores_flag=True,
    )

    direct = calculate_random_TIDE_scores_context_aware(
        dict(gene_dict), {"AC": {1: frozenset({"R1", "R2"})}}, complex_cache, "absmax",
        n_permutations=5, n_cpus=1, random_seed=0,
    )

    got = [float(v) for v in results.loc[0, "random_score_array"].split(";")]
    assert got == list(direct["AC"])