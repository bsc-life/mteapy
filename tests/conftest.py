"""Shared test fixtures: a small synthetic COBRApy model.

Network (all reactions irreversible, lb=0 unless noted):

    X[e] --EX_X (boundary)--> (sink/source, bounds -1000..1000)
    X[e] --R_TR--> A[c]                      gene: g_tr
    A[c] --R1--> B[c]                        gene: g1
    B[c] --R2--> C[c]                        gene: g2
    A[c] --R3--> C[c]                        gene: g3   (alt. direct route A->C)
    A[c] --R4--> D[c]                        gene: g4   (only route to D)

So a task requiring A[c] -> C[c] has two redundant routes (R1+R2, or R3
alone) and no single essential gene, while a task requiring A[c] -> D[c]
has exactly one route (R4) and g4 is essential for it.
"""

import pytest
from cobra.core import Metabolite, Model, Reaction


@pytest.fixture
def toy_model() -> Model:
    model = Model("toy_model")

    mets = {}
    for name, compartment in [("X", "e"), ("A", "c"), ("B", "c"), ("C", "c"), ("D", "c")]:
        met = Metabolite(f"{name}_{compartment}", name=name, compartment=compartment)
        mets[f"{name}[{compartment}]"] = met

    ex_x = Reaction("EX_X")
    ex_x.add_metabolites({mets["X[e]"]: -1})
    ex_x.bounds = (-1000, 1000)

    r_tr = Reaction("R_TR")
    r_tr.add_metabolites({mets["X[e]"]: -1, mets["A[c]"]: 1})
    r_tr.bounds = (0, 1000)
    r_tr.gene_reaction_rule = "g_tr"

    r1 = Reaction("R1")
    r1.add_metabolites({mets["A[c]"]: -1, mets["B[c]"]: 1})
    r1.bounds = (0, 1000)
    r1.gene_reaction_rule = "g1"

    r2 = Reaction("R2")
    r2.add_metabolites({mets["B[c]"]: -1, mets["C[c]"]: 1})
    r2.bounds = (0, 1000)
    r2.gene_reaction_rule = "g2"

    r3 = Reaction("R3")
    r3.add_metabolites({mets["A[c]"]: -1, mets["C[c]"]: 1})
    r3.bounds = (0, 1000)
    r3.gene_reaction_rule = "g3"

    r4 = Reaction("R4")
    r4.add_metabolites({mets["A[c]"]: -1, mets["D[c]"]: 1})
    r4.bounds = (0, 1000)
    r4.gene_reaction_rule = "g4"

    model.add_reactions([ex_x, r_tr, r1, r2, r3, r4])
    model.objective = "R_TR"

    return model
