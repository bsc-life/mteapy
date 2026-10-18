import numpy as np
import pandas as pd

from mteapy import inputcheck


def test_a_log_fold_change_column_passes():
    errors, warnings = inputcheck.check_lfc(pd.Series([2.0, -1.0, 0.3, -4.0]))
    assert errors == [] and warnings == []


def test_non_negative_values_look_like_expression_and_warn():
    errors, warnings = inputcheck.check_lfc(pd.Series([0.0, 3.0, 12.5, 400.0]), "baseMean")
    assert errors == [] and len(warnings) == 1
    assert "'baseMean' has no negative values" in warnings[0] and "expression" in warnings[0]


def test_a_single_negative_value_is_enough_to_not_warn_about_sign():
    assert inputcheck.check_lfc(pd.Series([0.0, 3.0, -0.1]))[1] == []


def test_huge_fold_changes_warn():
    errors, warnings = inputcheck.check_lfc(pd.Series([-500.0, 800.0, 300.0, -250.0] * 10))
    assert errors == [] and "very large" in warnings[0]


def test_unusable_columns_are_errors():
    assert "no numeric" in inputcheck.check_lfc(pd.Series([np.nan, np.inf]))[0][0]
    assert "nothing to score" in inputcheck.check_lfc(pd.Series([0.0, 0.0, 0.0]))[0][0]


def test_pvalue_range():
    assert inputcheck.check_pvalues(pd.Series([0.0, 0.5, 1.0, np.nan])) == ([], [])
    assert "outside [0, 1]" in inputcheck.check_pvalues(pd.Series([0.1, 12.0]))[0][0]


def test_expression_with_negative_values_warns_and_names_the_columns():
    frame = pd.DataFrame({"a": [1.0, 2.0], "b": [-1.0, 3.0], "c": [0.0, 5.0]})
    errors, warnings = inputcheck.check_expression(frame, "CellFie")
    assert errors == [] and "'b'" in warnings[0] and "CellFie" in warnings[0] and "'a'" not in warnings[0]
    assert inputcheck.check_expression(frame[["a", "c"]]) == ([], [])
