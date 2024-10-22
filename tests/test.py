import pytest

from cobra.core.gene import GPR
from mteapy.utils import safe_eval_gpr


@pytest.mark.parametrize("input,expected", [
    ((GPR().from_string("A and B and C"), {"A": 1, "B": 1.3, "C": -0.9}), -0.9),
    ((GPR().from_string("A or B or C"), {"A": 1, "B": 1.3, "C": -0.9}), 1.3)
])
def test_safe_eval_gpr(input, expected):
    expr, gene_dict = input
    assert safe_eval_gpr(expr, gene_dict, or_func="absmax") == expected