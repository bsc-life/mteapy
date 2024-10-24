import os
import pytest
import pandas as pd
import numpy as np

from cobra.core.gene import GPR
from mteapy.utils import mask_lfc_values, absmax, map_gpr, map_gpr_w_names, calculate_pvalue


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


@pytest.mark.parametrize("input,expected", [
    ((0, np.array([0, 0, 0, 0, 0])), 1.0),
    ((0, np.array([0, 1, 0, 0, 0])), 0.8),
    ((1, np.array([0, 0, 0, 0, 0])), 0.0)
])
def test_calculate_pvalue(input, expected):
    score, array = input
    assert calculate_pvalue(score, array) == expected