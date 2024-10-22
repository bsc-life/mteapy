import os
import pytest
import pandas as pd
import numpy as np

from cobra.core.gene import GPR
from mteapy.utils import mask_lfc_values, absmax, safe_eval_gpr, calculate_pvalue


def test_mask_lfc_values():
    curdir = os.path.dirname(os.path.realpath(__file__))
    input_df = pd.read_csv(os.path.join(curdir, "./data/mask_lfc_test_input.tsv"), sep="\t")
    expected_df = pd.read_csv(os.path.join(curdir, "./data/mask_lfc_test_expected.tsv"), sep="\t")
    assert all(mask_lfc_values(input_df, "lfc", "pvalue", 0.05) == expected_df)


@pytest.mark.parametrize("input,expected", [
    (np.array([-0.9, -0.5, -0.1]), -0.9),
    (np.array([0.1, 1.3, -1]), 1.3)
])
def test_absolute_minmax(input, expected):
    assert absmax(input) == expected


@pytest.mark.parametrize("input,expected", [
    (("A and B and C", {"A": 1, "B": 1.3, "C": -0.9}), -0.9),
    (("A or B or C", {"A": 1, "B": 1.3, "C": -0.9}), 1.3),
    (("A or B or C", {"A": -1, "B": -1.3, "C": -0.9}), -1.3)
])
def test_safe_eval_gpr(input, expected):
    expr, gene_dict = input
    assert safe_eval_gpr(GPR.from_string(expr), gene_dict, or_func="absmax") == expected


def test_calculate_pvalue():
    # TODO
    pass