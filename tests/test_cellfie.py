import pytest
import pandas as pd

from mteapy.cellfie import calculate_GAL


def test_calculate_GAL_rejects_unsupported_thresh_type():
    expr_data = pd.DataFrame({"S1": [1.0, 2.0]}, index=["G1", "G2"])

    with pytest.raises(ValueError):
        calculate_GAL(expr_data, thresh_type="Local")
