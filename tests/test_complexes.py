"""mteapy.complexes is now a thin re-export of cobra_netgraph.gpr (see that
module's own test suite for the full GPR-decomposition coverage) -- this
just confirms the re-export itself is wired correctly."""

from mteapy.complexes import UnsupportedGPRExpression, get_enzymes
from cobra_netgraph.gpr import UnsupportedGPRExpression as _RealUnsupportedGPRExpression
from cobra_netgraph.gpr import get_enzymes as _real_get_enzymes


def test_get_enzymes_is_the_cobra_netgraph_implementation():
    assert get_enzymes is _real_get_enzymes
    assert get_enzymes("A and B") == (("A", "B"),)


def test_unsupported_gpr_expression_is_the_cobra_netgraph_implementation():
    assert UnsupportedGPRExpression is _RealUnsupportedGPRExpression
