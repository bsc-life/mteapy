import pytest

from mteapy.complexes import get_enzymes


@pytest.mark.parametrize("input,expected", [
    ("", ()),
    ("A", (("A",),)),
    ("A and B", (("A", "B"),)),
    ("A or B", (("A",), ("B",))),
    ("A or A", (("A",),)),                                   # exact duplicate collapses
    ("A or (A and B)", (("A",),)),                           # absorption: (A,B) is not minimal
    ("(A and B) or (A and C)", (("A", "B"), ("A", "C"))),
    ("A and (B or C)", (("A", "B"), ("A", "C"))),             # requires distributive expansion
    ("A and (B or C) and D", (("A", "B", "D"), ("A", "C", "D"))),
    ("(A and B) or (C and (D or E))", (("A", "B"), ("C", "D"), ("C", "E"))),
    ("GENE-WITH-DASH or OTHER", (("GENE-WITH-DASH",), ("OTHER",))),
])
def test_get_enzymes(input, expected):
    assert get_enzymes(input) == expected
