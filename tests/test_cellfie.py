import pytest
import pandas as pd
import numpy as np

from mteapy.cellfie import (
    calculate_GAL, calculate_RAL, calculate_CellFie_scores,
    calculate_CellFie_scores_context_aware, compute_CellFie,
)


def test_calculate_GAL_rejects_unsupported_thresh_type():
    expr_data = pd.DataFrame({"S1": [1.0, 2.0]}, index=["G1", "G2"])

    with pytest.raises(ValueError):
        calculate_GAL(expr_data, thresh_type="Local")


def test_calculate_GAL_percentile_matches_matlab_log_space_by_default():
    # The original MATLAB CellFie computes percentile thresholds by taking
    # the percentile in log10 space and exponentiating back
    # (10**quantile(log10(x), p)), not a percentile taken directly on
    # linear x -- these give different values for skewed data. Sanity-check
    # the two really do differ for this input, then confirm the default
    # (log_transformed=False) matches the log-space computation.
    expr_data = pd.DataFrame({"S1": [1.0, 10.0, 100.0, 1000.0]}, index=["G1", "G2", "G3", "G4"])
    linear_data = expr_data.to_numpy().flatten()

    log_space_threshold = 10 ** np.quantile(np.log10(linear_data), 0.5)
    linear_space_threshold = np.quantile(linear_data, 0.5)
    assert not np.isclose(log_space_threshold, linear_space_threshold)

    gal = calculate_GAL(expr_data, thresh_type="global", global_thresh_type="percentile", global_value=0.5)
    expected = 5 * np.log(1 + linear_data / log_space_threshold)
    assert np.allclose(gal.to_numpy().flatten(), expected)


def test_calculate_GAL_log_transformed_computes_percentile_directly():
    expr_data = pd.DataFrame({"S1": [1.0, 10.0, 100.0, 1000.0]}, index=["G1", "G2", "G3", "G4"])
    linear_data = expr_data.to_numpy().flatten()

    gal = calculate_GAL(expr_data, thresh_type="global", global_thresh_type="percentile", global_value=0.5,
                         log_transformed=True)
    threshold = np.quantile(linear_data, 0.5)
    expected = 5 * np.log(1 + linear_data / threshold)
    assert np.allclose(gal.to_numpy().flatten(), expected)


def test_calculate_RAL_reaction_with_unmeasured_gene_is_nan(toy_model):
    gpr_dict = {r.id: r.gpr for r in toy_model.reactions if r.gene_reaction_rule}
    gal_df = pd.DataFrame({"S1": [5.0, 3.0]}, index=["g1", "g2"])  # g3/g4/g_tr not measured

    ral_df = calculate_RAL(gal_df, gpr_dict)

    assert ral_df.loc["R1", "S1"] == 5.0
    assert ral_df.loc["R2", "S1"] == 3.0
    assert np.isnan(ral_df.loc["R3", "S1"])
    assert np.isnan(ral_df.loc["R4", "S1"])


def test_calculate_CellFie_scores_excludes_no_data_reactions_from_the_average(toy_model):
    # Old (buggy) behavior would zero-fill R3 and average 5.0 with 0.0 -> 2.5.
    # Original CellFie excludes a no-data reaction from the average entirely.
    gpr_dict = {r.id: r.gpr for r in toy_model.reactions if r.gene_reaction_rule}
    gal_df = pd.DataFrame({"S1": [5.0, 3.0]}, index=["g1", "g2"])
    ral_df = calculate_RAL(gal_df, gpr_dict)

    task_structure = pd.DataFrame(
        {"T1": [True, False, True, False, False]},
        index=["R1", "R2", "R3", "R4", "R_TR"],
    )
    scores_df, binary_df = calculate_CellFie_scores(ral_df, task_structure)

    assert scores_df.loc["T1", "S1"] == 5.0


def test_calculate_CellFie_scores_task_with_no_data_at_all_is_nan(toy_model):
    gpr_dict = {r.id: r.gpr for r in toy_model.reactions if r.gene_reaction_rule}
    gal_df = pd.DataFrame({"S1": [5.0]}, index=["g1"])  # only g1 measured
    ral_df = calculate_RAL(gal_df, gpr_dict)

    task_structure = pd.DataFrame(
        {"T1": [False, False, True, False, False]},  # only R3, which has no data
        index=["R1", "R2", "R3", "R4", "R_TR"],
    )
    scores_df, binary_df = calculate_CellFie_scores(ral_df, task_structure)

    assert np.isnan(scores_df.loc["T1", "S1"])
    assert np.isnan(binary_df.loc["T1", "S1"])


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
