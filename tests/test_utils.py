import os
import pytest
import pandas as pd
import numpy as np

from cobra.core.gene import GPR
from mteapy.utils import mask_lfc_values, absmax, map_gpr, map_gpr_w_names, calculate_pvalue, add_task_metadata


def test_mask_lfc_values():
    curdir = os.path.dirname(os.path.realpath(__file__))
    input_df = pd.read_csv(os.path.join(curdir, "./data/mask_lfc_test_input.tsv"), sep="\t")
    expected_df = pd.read_csv(os.path.join(curdir, "./data/mask_lfc_test_expected.tsv"), sep="\t")
    assert all(mask_lfc_values(input_df, "lfc", "pvalue", 0.05) == expected_df)


@pytest.mark.parametrize("input,expected", [
    (np.array([-0.9, -0.5, -0.1]), -0.9),
    (np.array([0.1, 1.3, -1]), 1.3),
    (np.array([0.1, 1.3, 1]), 1.3)
])
def test_absmax(input, expected):
    assert absmax(input) == expected


@pytest.mark.parametrize("input,expected", [
    (("A and B and C", {"A": 1, "B": 1.3, "C": -0.9}), -0.9),
    (("A or B or C", {"A": 1, "B": 1.3, "C": -0.9}), 1.3),
    (("A or B or C", {"A": -1, "B": -1.3, "C": -0.9}), -1.3),
    (("A or B or (C and D)", {"A": -1, "B": -1.3, "C": -0.9, "D": 1.3}), -1.3),
    (("(A or B) and (C and D)", {"A": 1, "B": 1.3, "C": -0.9, "D": -1.3}), -1.3),
])
def test_map_gpr(input, expected):
    expr, gene_dict = input
    assert map_gpr(GPR.from_string(expr), gene_dict, or_func="absmax") == expected


@pytest.mark.parametrize("input,expected", [
    (("A and B and C", {"A": 1, "B": 1.3, "C": -0.9}), (-0.9, "C")),
    (("A or B or C", {"A": 1, "B": 1.3, "C": -0.9}), (1.3, "B")),
    (("A or B or C", {"A": -1, "B": -1.3, "C": -0.9}), (-0.9, "C")),
    (("A or B or (C and D)", {"A": -1, "B": -1.3, "C": -0.9, "D": 1.3}), (-0.9, "C")),
    (("(A or B) and (C and D)", {"A": 1, "B": 1.3, "C": -0.9, "D": -1.3}), (-1.3, "D")),
])
def test_map_gpr_w_names(input, expected):
    expr, gene_dict = input
    assert map_gpr_w_names(GPR.from_string(expr), gene_dict) == expected


def test_map_gpr_w_names_and_clause_with_one_missing_gene_has_no_data():
    # A complex needs every subunit measured to be assessed at all --
    # matching the original MATLAB CellFie's -1-sentinel-propagates-
    # through-AND semantics (a partial complex assessment using only the
    # measured subunit would be more lenient than the published algorithm).
    score, gene = map_gpr_w_names(GPR.from_string("A and B"), {"A": 5.0})
    assert score is None
    assert gene == "B"


def test_map_gpr_w_names_or_clause_recovers_past_one_missing_gene():
    # An isoenzyme option missing data doesn't sink the whole OR as long as
    # another option has real data -- max naturally skips it.
    score, gene = map_gpr_w_names(GPR.from_string("A or B"), {"B": 5.0})
    assert (score, gene) == (5.0, "B")


def test_map_gpr_w_names_or_clause_all_missing_has_no_data():
    score, gene = map_gpr_w_names(GPR.from_string("A or B"), {})
    assert score is None


def test_map_gpr_w_names_single_missing_gene_has_no_data():
    score, gene = map_gpr_w_names(GPR.from_string("A"), {})
    assert score is None
    assert gene == "A"


def test_map_gpr_w_names_none_gpr_has_no_data():
    assert map_gpr_w_names(None, {"A": 5.0}) == (None, "0")


@pytest.mark.parametrize("input,expected", [
    ((0, np.array([0, 0, 0, 0, 0])), 1.0),
    ((0, np.array([0, 1, 0, 0, 0])), 0.8),
    ((1, np.array([0, 0, 0, 0, 0])), 0.0)
])
def test_calculate_pvalue(input, expected):
    score, array = input
    assert calculate_pvalue(score, array) == expected


def test_add_task_metadata_keeps_result_rows_without_matching_metadata():
    results_df = pd.DataFrame({"task_id": ["1", "2"], "score": [0.1, 0.2]})
    task_metadata_df = pd.DataFrame({
        "ID": ["1"],
        "SYSTEM": ["amino acid metabolism"],
        "DESCRIPTION": ["some task"],
        "SUBSYSTEM": ["some subsystem"],
    })

    annotated_df = add_task_metadata(results_df, task_metadata_df)

    assert list(annotated_df["task_id"]) == ["1", "2"]
    assert annotated_df.loc[annotated_df["task_id"] == "2", "task_description"].isna().all()