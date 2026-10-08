"""Build a directed metabolite-reaction network for one task route, scored
against one sample's expression.

This is the task-scoring-specific layer on top of the generic
`cobra_netgraph.bipartite.build_bipartite_graph`: it supplies the two
decorators that attach `mteapy.context_scoring`'s evidence/score fields to
each reaction node and task input/output tagging to each metabolite node.
The graph-shape logic itself (currency-metabolite duplication, BiGG/EC
lookup, flux-direction edges) lives in `cobra_netgraph` and has no
dependency on tasks or gene expression at all -- see that package if you
need a bipartite network for something other than a metabolic task (e.g. a
reconstruction/gapfilling pipeline's own reaction-status view).

Everything in this module that does NOT depend on sample expression --
which route uses which reactions, each reaction's GPR-derived candidate
complexes, its EC/BiGG annotations, and its solved flux for this specific
route -- is a property of the *model and task* alone. Scoring a route
against a sample's expression (`mteapy.context_scoring.assess_reaction`) is
comparatively cheap pure arithmetic once that structural information is
available. Callers that serve many samples against the same task/route
should reuse `compute_route_fluxes`'s result and a `build_complex_cache`
across requests rather than recomputing them per sample.
"""

from __future__ import annotations

import dataclasses

from cobra.core import Metabolite, Model
from cobra_netgraph.bipartite import (
    CURRENCY_DEGREE_THRESHOLD,
    CURRENCY_NAMES,
    build_bipartite_graph,
    classify_currency,
    ec_code_of,
    load_bigg_ids,
)
from cobra_netgraph.gpr import get_enzymes

from mteapy.context_scoring import assess_reaction
from mteapy.task_model import (
    _EQU_ARROW_PATTERN, _EQU_TERM_SPLIT_PATTERN, _parse_equation_term, build_metabolite_lookup,
    build_task_model, is_pseudo_reaction, set_min_total_flux_objective,
)
from mteapy.tasks import MetabolicTask

__all__ = [
    "build_route_graph", "compute_route_fluxes", "resolve_boundary_ids",
    "load_bigg_ids", "classify_currency", "ec_code_of",
    "CURRENCY_NAMES", "CURRENCY_DEGREE_THRESHOLD",
]


def resolve_boundary_ids(task: MetabolicTask, lookup: dict[str, str]) -> tuple[set[str], set[str]]:
    """The model metabolite ids for a task's declared IN/OUT boundary
    metabolites (e.g. "glucose[e]" -> "MAM01965e"), for input/output
    coloring. Deliberately not the EQU-constraint metabolites (e.g. the
    ATP/ADP/Pi/H+ bookkeeping equation most ATP-regeneration tasks use) --
    those are the same currency cofactors appearing throughout the route
    for unrelated reactions too, so tagging every occurrence would be
    misleading rather than clarifying.
    """
    input_ids = {lookup[bm.metabolite.strip()] for bm in task.inputs if bm.metabolite.strip() in lookup}
    output_ids = {lookup[bm.metabolite.strip()] for bm in task.outputs if bm.metabolite.strip() in lookup}
    return input_ids, output_ids


def compute_route_fluxes(model: Model, task: MetabolicTask, reactions: set[str]) -> dict[str, float]:
    """Solve for one signed flux value per reaction in `reactions`.

    Every other real (non-pseudo) reaction is closed to (0, 0) first, so
    the only degrees of freedom left are exactly this route's reactions
    plus the task's own pseudo source/demand/equation reactions -- the
    resulting flux vector is a genuine feasible solution for the task
    that uses only this route, not an arbitrary/independent FBA solve.
    The sign tells the true net direction of flow (which can run opposite
    the way the reaction's stoichiometry happens to be written), and the
    magnitude is what edge width is scaled by.

    This result depends only on the model, the task, and the route's
    reaction set -- never on any sample's expression -- so a caller
    serving many samples against the same route should compute it once
    and cache it, not re-solve per sample.
    """
    tmodel = build_task_model(model, task)
    for reaction in tmodel.reactions:
        if not is_pseudo_reaction(reaction.id) and reaction.id not in reactions:
            reaction.bounds = (0.0, 0.0)
    set_min_total_flux_objective(tmodel)
    solution = tmodel.optimize()
    return {rid: solution.fluxes[rid] for rid in reactions}


def compute_route_fluxes_submodel(model: Model, task: MetabolicTask, reactions: set[str],
                                  met_lookup: dict[str, str] | None = None) -> dict[str, float]:
    """Same result as `compute_route_fluxes`, but solved on a tiny model
    holding only this route's reactions (plus the task's pseudo reactions)
    instead of a copy of the whole genome-scale model.

    `compute_route_fluxes` closes every non-route reaction to (0, 0), so
    they cannot carry flux anyway; dropping them changes nothing about the
    feasible set, only how much model has to be copied and handed to the
    solver. Task metabolites the route's reactions don't touch are still
    added (as bare metabolites) so their pseudo source/demand reactions
    keep their bounds, exactly as in the full model. A `CHANGED RXN` entry
    for a reaction outside the route is dropped -- that reaction is closed
    in the full-model version regardless.
    """
    lookup = met_lookup if met_lookup is not None else build_metabolite_lookup(model)
    sub = Model(f"route_{task.id}")
    sub.add_reactions([model.reactions.get_by_id(rid).copy() for rid in reactions])

    tokens = [bm.metabolite for bm in task.inputs] + [bm.metabolite for bm in task.outputs]
    for eq in task.equations:
        arrow = _EQU_ARROW_PATTERN.search(eq.equation)
        for side in ((eq.equation[: arrow.start()], eq.equation[arrow.end():]) if arrow else ()):
            for term in _EQU_TERM_SPLIT_PATTERN.split(side):
                if term.strip():
                    tokens.append(_parse_equation_term(term.strip())[1])
    for token in tokens:
        met_id = lookup.get(token.strip())
        if met_id is not None and met_id not in sub.metabolites:
            src = model.metabolites.get_by_id(met_id)
            sub.add_metabolites([Metabolite(src.id, formula=src.formula, name=src.name, compartment=src.compartment)])

    sub_task = dataclasses.replace(
        task, changed_bounds=[cb for cb in task.changed_bounds if cb.reaction_id in sub.reactions]
    )
    tmodel = build_task_model(sub, sub_task, lookup)
    set_min_total_flux_objective(tmodel)
    solution = tmodel.optimize()
    return {rid: solution.fluxes[rid] for rid in reactions}


def build_route_graph(model: Model, reactions: set[str], gene_dict: dict[str, float], fluxes: dict[str, float],
                       met_bigg: dict[str, str] | None = None, rxn_bigg: dict[str, str] | None = None,
                       input_ids: set[str] | None = None, output_ids: set[str] | None = None,
                       or_func: str = "max") -> dict:
    """Directed bipartite metabolite-reaction graph for one route, scored
    against one sample's `gene_dict`.

    A thin task-scoring-specific wrapper over
    `cobra_netgraph.bipartite.build_bipartite_graph`: it decorates each
    reaction node with `mteapy.context_scoring.assess_reaction`'s
    evidence/score/complex breakdown, and each metabolite node with
    `io` ("input"/"output"/"internal", from `input_ids`/`output_ids` --
    see `resolve_boundary_ids`). Currency-metabolite duplication, BiGG/EC
    lookup, and flux-direction edges are all handled by the generic layer
    -- see that module's docstring for the full behavior.

    `met_bigg`/`rxn_bigg` (from `load_bigg_ids`) and `input_ids`/`output_ids`
    (from `resolve_boundary_ids`) are optional: omit them (e.g. for a model
    without a BiGG annotation table, or a task with no declared boundary)
    and every node just gets a blank `bigg_id` / `io="internal"`.

    `or_func` is passed straight through to `assess_reaction`: "max"
    (default) for a non-negative signal like expression, "absmax" for one
    that can be negative like a log-fold-change. Getting this wrong doesn't
    raise -- it silently mis-scores/mis-classifies evidence for a signed
    `gene_dict`, so a caller visualizing differential-expression data must
    pass `or_func="absmax"` explicitly (see `mteapy.context_scoring
    .assess_reaction`'s docstring for why the two need different rules).
    """
    input_ids = input_ids or set()
    output_ids = output_ids or set()

    def reaction_decorator(r, model):
        complexes = get_enzymes(r.gene_reaction_rule)
        assessment = assess_reaction(r.id, complexes, gene_dict, or_func=or_func)
        gene_values = {g: gene_dict.get(g, 0.0) for c in complexes for g in c}
        complex_scores = [assessment.complex_scores[c] for c in complexes]
        return dict(
            complexes=[list(c) for c in complexes], complex_scores=complex_scores,
            gene_values=gene_values, score=assessment.score, evidence=assessment.evidence.value,
            winning_complexes=[list(c) for c in assessment.winning_complexes],
        )

    def metabolite_decorator(m, model):
        if m.id in output_ids:
            io = "output"
        elif m.id in input_ids:
            io = "input"
        else:
            io = "internal"
        return {"io": io}

    return build_bipartite_graph(
        model, reactions, fluxes=fluxes, met_bigg=met_bigg, rxn_bigg=rxn_bigg,
        reaction_decorator=reaction_decorator, metabolite_decorator=metabolite_decorator,
    )
