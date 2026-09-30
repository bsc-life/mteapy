import pytest
import pandas as pd
import numpy as np

from mteapy.cellfie import calculate_GAL, calculate_CellFie_scores_context_aware, compute_CellFie


def test_calculate_GAL_rejects_unsupported_thresh_type():
    expr_data = pd.DataFrame({"S1": [1.0, 2.0]}, index=["G1", "G2"])

    with pytest.raises(ValueError):
        calculate_GAL(expr_data, thresh_type="Local")


def test_calculate_CellFie_scores_context_aware_picks_the_max_route(toy_model):
    # R1+R2 (g1=1, g2=1) has mean score 1; R3 alone (g3=3) has mean score 3.
    # or_func="max" (CellFie's convention, non-negative GALs) picks {R3}.
    tasks_routes = {"AC": {1: frozenset({"R1", "R2"}), 2: frozenset({"R3"})}}
    gal_df = pd.DataFrame({"S1": [1.0, 1.0, 3.0]}, index=["g1", "g2", "g3"])

    scores_df, binary_df = calculate_CellFie_scores_context_aware(gal_df, tasks_routes, toy_model)

    assert scores_df.loc["AC", "S1"] == 3.0
    assert binary_df.loc["AC", "S1"] == int(3.0 >= 5 * np.log(2))


def test_compute_CellFie_context_aware_requires_tasks_routes(toy_model):
    expr_data = pd.DataFrame({"S1": [1.0, 1.0, 3.0]}, index=["g1", "g2", "g3"])
    with pytest.raises(ValueError):
        compute_CellFie(expr_data, None, toy_model, mapping_strategy="context-aware")
