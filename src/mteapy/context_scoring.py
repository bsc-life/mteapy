"""Two-dimensional ("context-aware") task-activity scoring.

CellFie/TIDE-style scoring collapses a task to a single reaction set (via
one pFBA solve) and each reaction's GPR to a single flattened number (via
`mteapy.utils.map_gpr`). This module scores both dimensions of degeneracy
that approach discards:

- **Topology**: a task can admit several, equally optimal *routes* (distinct
  reaction sets -- see `mteapy.enumeration`/`mteapy.routes`). Which one
  actually matches a sample's biology is itself a question, not a given.
- **Regulation**: within one reaction, its GPR can itself offer several
  candidate enzyme complexes (`mteapy.complexes.get_enzymes`) -- which one a
  sample's expression actually supports is a second, independent question.

Pipeline per (task, sample):

1. For each of the task's candidate routes, and each reaction in it,
   decompose the reaction's GPR into candidate complexes.
2. Score each complex as the min of its member genes' expression (an AND
   requirement), and take a reaction's score as the max across its own
   candidate complexes (an OR choice) -- together this reproduces exactly
   what `map_gpr` would compute for that reaction alone, while also
   recording *which* complex achieved the max, which `map_gpr` discards.
3. Aggregate a route's score from its reactions' scores: "min" (the
   strict weakest-link choice this project's own route-tie analysis
   validated) or "median" (offered as the more permissive alternative
   that analysis also tried, but found far too tie-insensitive to trust
   alone -- not the default).
4. A task's winning route(s) for that sample are whichever are tied for
   the top score. Ties are returned explicitly, never silently broken:
   this project found that a naive `max(dict, key=dict.get)` on tied
   scores systematically (and invisibly) favoured whichever route
   happened to be enumerated first -- see `tied_argmax` below.

`score_task` returns the full per-route, per-reaction, per-complex
breakdown for one (task, sample) -- the interpretability payoff of this
whole approach. `score_tasks_matrix` wraps it into a plain tasks x samples
score matrix, the same shape CellFie/TIDE-style scores use, for drop-in
comparison or downstream differential-activity testing.
"""

from __future__ import annotations

import math
import statistics
from dataclasses import dataclass, field

import pandas as pd
from cobra.core import Model

from mteapy.complexes import get_enzymes

_SUPPORTED_AGGREGATIONS = ("min", "median")


def tied_argmax(scores: dict, rel_tol: float = 1e-9, abs_tol: float = 1e-9) -> tuple:
    """All keys tied for the max value in `scores`, within floating tolerance.

    Returns an empty tuple for an empty `scores` (e.g. a route with no
    reactions), a single-element tuple for a clear winner, and a
    multi-element tuple when several keys are genuinely tied -- callers
    must not just take `[0]` and assume it means anything on its own; check
    `len(...)` (or use `TaskActivityResult.is_tied`) first.
    """
    if not scores:
        return ()
    max_val = max(scores.values())
    return tuple(sorted(k for k, v in scores.items() if math.isclose(v, max_val, rel_tol=rel_tol, abs_tol=abs_tol)))


def build_complex_cache(model: Model, reaction_ids) -> dict[str, tuple[tuple[str, ...], ...]]:
    """Precompute reaction_id -> candidate complexes for every id in `reaction_ids`.

    A reaction's GPR structure doesn't depend on the sample being scored,
    so this is computed once and reused across every sample scored against
    the same route set, rather than re-decomposing the same GPR per sample
    (the mistake that made an earlier version of this pipeline ~2 hours
    slower than it needed to be).
    """
    cache = {}
    for rid in reaction_ids:
        gpr = model.reactions.get_by_id(rid).gene_reaction_rule
        cache[rid] = get_enzymes(gpr)
    return cache


def score_reaction(
    complexes: tuple[tuple[str, ...], ...], gene_dict: dict[str, float]
) -> tuple[float, tuple[str, ...] | None]:
    """Score one reaction's candidate complexes against `gene_dict`.

    Returns `(reaction_score, winning_complex)`. `winning_complex` is the
    single tied_argmax winner, or `None` if the complex-level winner is
    itself ambiguous (two-plus complexes tied) -- the reaction's own score
    is still perfectly well-defined in that case, only the attribution of
    *which* complex earned it is ambiguous.
    """
    if not complexes:
        return 0.0, None
    scores = {c: min(gene_dict.get(g, 0.0) for g in c) for c in complexes}
    winners = tied_argmax(scores)
    winning_complex = winners[0] if len(winners) == 1 else None
    return scores[winners[0]], winning_complex


def _aggregate(values: list[float], aggregation: str) -> float:
    if not values:
        return 0.0
    if aggregation == "min":
        return min(values)
    if aggregation == "median":
        return statistics.median(values)
    raise ValueError(f"Unsupported aggregation {aggregation!r}; use one of {_SUPPORTED_AGGREGATIONS}")


@dataclass
class RouteScore:
    """One route's score for one sample, with its full reaction-level breakdown."""

    route_id: int
    score: float
    reaction_scores: dict[str, float] = field(default_factory=dict)
    # None where the complex-level winner was itself ambiguous for that reaction.
    reaction_complex: dict[str, tuple[str, ...] | None] = field(default_factory=dict)


@dataclass
class TaskActivityResult:
    """One task's context-aware activity result for one sample."""

    task_id: str
    winning_route_ids: tuple[int, ...]
    score: float
    route_scores: dict[int, RouteScore] = field(default_factory=dict)

    @property
    def is_tied(self) -> bool:
        """True when more than one route is tied for the top score.

        Route-level ties turned out to be the norm, not the exception (see
        module docstring) -- always check this before treating
        `winning_route_ids[0]` as *the* answer.
        """
        return len(self.winning_route_ids) > 1


def score_task(
    task_routes: dict[int, frozenset[str]],
    complex_cache: dict[str, tuple[tuple[str, ...], ...]],
    gene_dict: dict[str, float],
    aggregation: str = "min",
    task_id: str = "",
) -> TaskActivityResult:
    """Score every candidate route of one task against one sample's expression.

    Parameters
    ----------
    task_routes:
        `{route_id: reaction_id_set}` for one task, e.g. from
        `mteapy.routes.load_task_routes` (or `load_multiroute_tasks`, keyed
        down to one task). A single-route task (a plain dict of one entry)
        works fine here -- there is simply nothing for the topology
        dimension to distinguish, and `is_tied` will be False.
    complex_cache:
        `{reaction_id: candidate_complexes}`, from `build_complex_cache`.
        Must cover every reaction appearing in `task_routes`.
    gene_dict:
        One sample's gene expression, `{gene_id: value}`.
    aggregation:
        How a route's reaction scores combine into one route score:
        "min" (default, strict weakest-link) or "median" (permissive).
    task_id:
        Carried through to the result only for the caller's convenience
        (e.g. building a results table); not used in scoring.
    """
    route_scores: dict[int, RouteScore] = {}
    for route_id, reactions in task_routes.items():
        reaction_scores: dict[str, float] = {}
        reaction_complex: dict[str, tuple[str, ...] | None] = {}
        for rid in reactions:
            complexes = complex_cache.get(rid, ())
            score, winner = score_reaction(complexes, gene_dict)
            reaction_scores[rid] = score
            reaction_complex[rid] = winner
        agg_score = _aggregate(list(reaction_scores.values()), aggregation)
        route_scores[route_id] = RouteScore(route_id, agg_score, reaction_scores, reaction_complex)

    top = tied_argmax({rid: rs.score for rid, rs in route_scores.items()})
    final_score = route_scores[top[0]].score if top else 0.0

    return TaskActivityResult(
        task_id=task_id,
        winning_route_ids=top,
        score=final_score,
        route_scores=route_scores,
    )


def score_tasks_matrix(
    tasks_routes: dict[str, dict[int, frozenset[str]]],
    model: Model,
    expr_df: pd.DataFrame,
    aggregation: str = "min",
) -> pd.DataFrame:
    """Score every task in `tasks_routes` against every sample (column) in `expr_df`.

    Returns a plain tasks x samples score matrix -- the same shape a
    CellFie/TIDE score matrix has, so it can be dropped into the same
    downstream differential-activity workflow. This is the scalar-only
    convenience wrapper; call `score_task` directly (per task, per sample)
    when the route/complex breakdown itself is wanted, not just the final
    number.

    Parameters
    ----------
    tasks_routes:
        `{task_id: {route_id: reaction_id_set}}` for every task to score,
        e.g. built from one or more calls to `mteapy.routes.load_task_routes`.
    model:
        The COBRApy model the routes' reaction ids come from, used only to
        look up each reaction's GPR for `build_complex_cache`.
    expr_df:
        Gene expression, genes (index) x samples (columns).
    aggregation:
        Passed through to `score_task` ("min" or "median").
    """
    all_reactions = sorted({r for routes in tasks_routes.values() for reactions in routes.values() for r in reactions})
    complex_cache = build_complex_cache(model, all_reactions)

    rows: dict[str, dict[str, float]] = {}
    for task_id, routes in tasks_routes.items():
        row = {}
        for sample in expr_df.columns:
            gene_dict = expr_df[sample].to_dict()
            result = score_task(routes, complex_cache, gene_dict, aggregation=aggregation, task_id=task_id)
            row[sample] = result.score
        rows[task_id] = row

    return pd.DataFrame(rows).T
