import pytest
import pandas as pd
import numpy as np

from cobra.core.gene import GPR
from mteapy.tide import calculate_TIDE_scores


@pytest.mark.parametrize("gene_dict,expected", [
    (dict(zip(["G" + str(i+1) for i in range(10)], [0,0,0,0,0,0,0,0,0,0])), [0, 0]),
    (dict(zip(["G" + str(i+1) for i in range(10)], [-5,-5,5,5,5,5,5,-5,0,0])), [0, -2.5]),
    (dict(zip(["G" + str(i+1) for i in range(10)], [0,-1,2,-1.2,-6,0.4,-0.9,-3.2,0,2.3])), [-3.5, -1.6])
])
def test_TIDE_scores(gene_dict, expected):
    task_structure = pd.DataFrame({"X1": [1,1,0,0], "X2": [0,0,1,1]}, index=["R1","R2","R3","R4"]).astype(bool)
    gpr_dict = dict(zip(
        task_structure.index,
        [GPR.from_string("G1 and G2"), GPR.from_string("G3 or G4 or G5"), \
         GPR.from_string("(G6 or G7) and G8"), GPR.from_string("G9 and G10")]
    ))

    assert all(calculate_TIDE_scores(gene_dict, task_structure, gpr_dict, or_func="absmax") == expected)