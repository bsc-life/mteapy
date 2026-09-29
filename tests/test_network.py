"""Tests for mteapy.network: the task-scoring-specific decorators over
`cobra_netgraph.bipartite.build_bipartite_graph`.

The generic graph-shape behavior (currency-metabolite duplication, flux-
direction edges, BiGG/EC lookup) is delegated entirely to `cobra_netgraph`
and tested there -- this file only covers what mteapy's own wrapper adds:
evidence/score decoration on reaction nodes, io tagging on metabolite
nodes, `or_func` threading, and the task-specific
`compute_route_fluxes`/`resolve_boundary_ids` helpers. One integration
test (`test_build_route_graph_duplicates_currency_and_tags_boundary_io`)
exercises the generic delegation and the mteapy-specific decoration
together, to catch a regression in the wiring between them.
"""

import pytest
from cobra.core import Metabolite, Model, Reaction

from mteapy.context_scoring import ReactionEvidence
from mteapy.network import build_route_graph, compute_route_fluxes, resolve_boundary_ids
from mteapy.task_model import build_metabolite_lookup
from mteapy.tasks import BoundedMetabolite, MetabolicTask


def _task(inputs=(), outputs=()):
    return MetabolicTask(
        id="t1", description="test task", should_fail=False,
        inputs=list(inputs), outputs=list(outputs), equations=[], changed_bounds=[],
    )


def test_resolve_boundary_ids_skips_unresolvable_entries(toy_model):
    lookup = build_metabolite_lookup(toy_model)
    task = _task(
        inputs=[BoundedMetabolite(metabolite="X[e]", lower_bound=0, upper_bound=1000)],
        outputs=[
            BoundedMetabolite(metabolite="D[c]", lower_bound=0, upper_bound=1000),
            BoundedMetabolite(metabolite="nonexistent[z]", lower_bound=0, upper_bound=1000),
        ],
    )
    input_ids, output_ids = resolve_boundary_ids(task, lookup)
    assert input_ids == {"X_e"}
    assert output_ids == {"D_c"}  # the unresolvable OUT entry is silently skipped, not an error


def test_compute_route_fluxes_sign_and_magnitude(toy_model):
    task = _task(
        inputs=[BoundedMetabolite(metabolite="X[e]", lower_bound=0, upper_bound=1000)],
        outputs=[BoundedMetabolite(metabolite="D[c]", lower_bound=1, upper_bound=1)],
    )
    fluxes = compute_route_fluxes(toy_model, task, {"R_TR", "R4"})
    assert fluxes["R_TR"] == pytest.approx(1.0)
    assert fluxes["R4"] == pytest.approx(1.0)


def test_build_route_graph_duplicates_currency_and_tags_boundary_io():
    """Integration test: currency duplication (delegated to cobra_netgraph)
    and io tagging (mteapy's own metabolite_decorator) must both show up
    correctly on the same graph."""
    model = Model("m")
    a = Metabolite("a_c", name="A", compartment="c")
    b = Metabolite("b_c", name="B", compartment="c")
    c = Metabolite("c_c", name="C", compartment="c")
    atp = Metabolite("atp_c", name="ATP", compartment="c")
    adp = Metabolite("adp_c", name="ADP", compartment="c")

    r1 = Reaction("R1")  # A + ATP -> B + ADP  (currency ATP/ADP alongside real backbone A/B)
    r1.add_metabolites({a: -1, atp: -1, b: 1, adp: 1})
    r1.gene_reaction_rule = ""
    r2 = Reaction("R2")  # B + ATP -> C + ADP
    r2.add_metabolites({b: -1, atp: -1, c: 1, adp: 1})
    r2.gene_reaction_rule = ""
    model.add_reactions([r1, r2])

    fluxes = {"R1": 1.0, "R2": 1.0}
    graph = build_route_graph(
        model, {"R1", "R2"}, gene_dict={}, fluxes=fluxes,
        input_ids={"a_c"}, output_ids={"c_c"},
    )

    met_nodes = [n for n in graph["nodes"] if n["type"] == "metabolite"]
    atp_nodes = [n for n in met_nodes if n["full_id"] == "atp_c"]
    backbone_nodes = [n for n in met_nodes if n["full_id"] in ("a_c", "b_c", "c_c")]

    # ATP is duplicated once per reaction that uses it (delegated to cobra_netgraph).
    assert len(atp_nodes) == 2
    assert all(n["currency"] for n in atp_nodes)
    assert not any(n["currency"] for n in backbone_nodes)

    io_by_id = {n["full_id"]: n["io"] for n in met_nodes}
    assert io_by_id["a_c"] == "input"
    assert io_by_id["c_c"] == "output"
    assert io_by_id["b_c"] == "internal"


def test_build_route_graph_decorates_reaction_with_full_evidence_breakdown():
    """The reaction_decorator must attach the complete
    mteapy.context_scoring breakdown (not just a bare score), matching
    what a caller building a gene/complex detail panel needs."""
    model = Model("m")
    a = Metabolite("a_c", name="A", compartment="c")
    b = Metabolite("b_c", name="B", compartment="c")
    r = Reaction("R1")
    r.add_metabolites({a: -1, b: 1})
    r.gene_reaction_rule = "g1 and g2"
    model.add_reactions([r])

    graph = build_route_graph(model, {"R1"}, gene_dict={"g1": 5.0, "g2": 2.0}, fluxes={"R1": 1.0})
    rxn_node = next(n for n in graph["nodes"] if n["type"] == "reaction")

    assert rxn_node["complexes"] == [["g1", "g2"]]
    assert rxn_node["complex_scores"] == [2.0]  # AND: min(5.0, 2.0)
    assert rxn_node["gene_values"] == {"g1": 5.0, "g2": 2.0}
    assert rxn_node["score"] == 2.0
    assert rxn_node["evidence"] == ReactionEvidence.SUPPORTED.value
    assert rxn_node["winning_complexes"] == [["g1", "g2"]]
    assert rxn_node["flux"] == 1.0  # still present, set by the generic layer


def test_build_route_graph_or_func_reaches_assess_reaction():
    """Regression test: build_route_graph must pass or_func through to
    assess_reaction, not silently default to "max" -- a reaction node
    scored against signed (log-fold-change) data needs "absmax" or its
    evidence/winning-complex classification comes out wrong."""
    model = Model("m")
    a = Metabolite("a_c", name="A", compartment="c")
    b = Metabolite("b_c", name="B", compartment="c")
    r = Reaction("R1")
    r.add_metabolites({a: -1, b: 1})
    r.gene_reaction_rule = "g1 or g2"
    model.add_reactions([r])

    # g1 is mildly up (+0.2), g2 is strongly down (-3.0). Under "max" g1
    # wins; under "absmax" g2 must win instead.
    gene_dict = {"g1": 0.2, "g2": -3.0}

    default_graph = build_route_graph(model, {"R1"}, gene_dict, fluxes={"R1": 1.0})
    rxn_node = next(n for n in default_graph["nodes"] if n["type"] == "reaction")
    assert rxn_node["score"] == 0.2
    assert rxn_node["winning_complexes"] == [["g1"]]

    absmax_graph = build_route_graph(model, {"R1"}, gene_dict, fluxes={"R1": 1.0}, or_func="absmax")
    rxn_node = next(n for n in absmax_graph["nodes"] if n["type"] == "reaction")
    assert rxn_node["score"] == -3.0
    assert rxn_node["winning_complexes"] == [["g2"]]
