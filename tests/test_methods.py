import json

import pandas as pd
import pytest

from mteapy import methods
from mteapy.context_scoring import RunCancelled, score_tasks_matrix, score_tasks_report
from mteapy.parser import mtea_parser


def test_methods_are_json_serializable_and_complete():
    described = methods.describe_methods()
    json.dumps(described)
    [tas] = [m for m in described if m["id"] == "TAS"]
    assert {p["name"] for p in tas["params"]} == {"aggregation", "or_func"}
    assert all(p["help"] and p["label"] for p in tas["params"])


def test_validate_params_fills_defaults_and_rejects_bad_input():
    tas = methods.get_method("TAS")
    assert methods.validate_params(tas, None) == {"aggregation": "min", "or_func": "max"}
    assert methods.validate_params(tas, {"aggregation": "mean"}) == {"aggregation": "mean", "or_func": "max"}
    with pytest.raises(ValueError, match="aggregation"):
        methods.validate_params(tas, {"aggregation": "sum"})
    with pytest.raises(ValueError, match="unknown parameter"):
        methods.validate_params(tas, {"threshold": 1})
    with pytest.raises(KeyError, match="available"):
        methods.get_method("nope")


def test_numeric_and_bool_params_are_coerced_and_range_checked():
    method = methods.Method("X", "X", "", (
        methods.Param("p", "P", "float", 0.5, "", min=0, max=1),
        methods.Param("n", "N", "int", 3, "", min=1),
        methods.Param("flag", "F", "bool", False, ""),
    ))
    assert methods.validate_params(method, {"p": "0.25", "n": 2.0, "flag": True}) == {"p": 0.25, "n": 2, "flag": True}
    for bad in ({"p": 2}, {"p": "abc"}, {"n": 1.5}, {"n": 0}, {"flag": "yes"}):
        with pytest.raises(ValueError):
            methods.validate_params(method, bad)


def test_cellfie_cli_options_match_the_spec():
    parser = mtea_parser()
    args = parser.parse_args(["analyze", "CellFie", "x.tsv"])
    for p in methods.CELLFIE.params:
        assert getattr(args, p.dest or p.name) == p.default
    args = parser.parse_args(["analyze", "CellFie", "x.tsv", "--threshold_type", "global", "--global_value", "0.9",
                              "--log_transformed"])
    assert (args.thresh_type, args.global_value, args.log_transformed) == ("global", 0.9, True)
    # every `when` condition refers to a real parameter and a legal value
    by_name = {p.name: p for p in methods.CELLFIE.params}
    for p in methods.CELLFIE.params:
        for name, value in p.when:
            assert value in by_name[name].choices


def test_cli_options_match_the_spec():
    """The CLI's TAS options must be exactly the spec's: same names, defaults and choices."""
    parser = mtea_parser()
    args = parser.parse_args(["analyze", "TAS", "x.tsv"])
    for p in methods.TAS.params:
        assert getattr(args, p.name) == p.default
    args = parser.parse_args(["analyze", "TAS", "x.tsv", "--aggregation", "mean", "--or_func", "absmax"])
    assert (args.aggregation, args.or_func) == ("mean", "absmax")
    with pytest.raises(SystemExit):
        parser.parse_args(["analyze", "TAS", "x.tsv", "--aggregation", "sum"])


# ---------------------------------------------------------------------------
# score_tasks_report
# ---------------------------------------------------------------------------

@pytest.fixture
def toy_inputs(toy_model):
    # Task "ac": A->C, two routes (R3, or R1+R2); task "ad": A->D, one route (R4).
    routes = {"ac": {1: frozenset({"R3"}), 2: frozenset({"R1", "R2"})}, "ad": {3: frozenset({"R4"})}}
    expr = pd.DataFrame({"s1": {"g1": 10, "g2": 10, "g3": 1, "g4": 5, "g_tr": 7},
                         "s2": {"g1": 1, "g2": 1, "g3": 8, "g4": 0, "g_tr": 7}})
    return routes, expr


def test_report_has_four_aligned_frames(toy_model, toy_inputs):
    routes, expr = toy_inputs
    report = score_tasks_report(routes, toy_model, expr)
    assert set(report) == {"scores", "complete", "tied", "winners"}
    for frame in report.values():
        assert list(frame.index) == ["ac", "ad"] and list(frame.columns) == ["s1", "s2"]
    assert report["winners"].loc["ac", "s1"] == (2,)      # R1+R2 (min 10) beats R3 (1)
    assert report["winners"].loc["ac", "s2"] == (1,)      # R3 (8) beats R1+R2 (min 1)
    assert report["scores"].loc["ac", "s1"] == 10 and report["scores"].loc["ac", "s2"] == 8
    assert not report["tied"].loc["ac", "s1"]


def test_matrix_wrapper_matches_the_report(toy_model, toy_inputs):
    routes, expr = toy_inputs
    scores, complete = score_tasks_matrix(routes, toy_model, expr)
    report = score_tasks_report(routes, toy_model, expr)
    pd.testing.assert_frame_equal(scores, report["scores"])
    pd.testing.assert_frame_equal(complete, report["complete"])


def test_progress_is_reported_per_sample_and_can_cancel(toy_model, toy_inputs):
    routes, expr = toy_inputs
    seen = []
    score_tasks_report(routes, toy_model, expr, progress=lambda done, total: seen.append((done, total)))
    assert seen == [(1, 2), (2, 2)]

    def stop(done, total):
        raise RunCancelled()
    with pytest.raises(RunCancelled):
        score_tasks_report(routes, toy_model, expr, progress=stop)
