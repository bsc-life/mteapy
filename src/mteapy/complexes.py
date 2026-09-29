"""GPR-to-enzyme-complex decomposition.

This is now a thin re-export: the actual implementation moved to
`cobra_netgraph.gpr`, a standalone package with no task-scoring dependency
(GPR decomposition is generic COBRApy-model logic, useful on its own --
e.g. to a reconstruction/gapfilling pipeline that has no notion of
"metabolic tasks" at all). Kept here, under the same names, so existing
code doing `from mteapy.complexes import get_enzymes` keeps working
unchanged.
"""
from cobra_netgraph.gpr import UnsupportedGPRExpression, get_enzymes

__all__ = ["get_enzymes", "UnsupportedGPRExpression"]
