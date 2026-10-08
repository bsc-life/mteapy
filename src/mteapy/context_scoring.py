"""Two-dimensional ("context-aware") task-activity scoring.

CellFie/TIDE-style scoring collapses a task to a single reaction set (via
one pFBA solve) and each reaction's GPR to a single flattened number (via
`mteapy.utils.map_gpr`). This module scores both dimensions of degeneracy
that approach discards:

- **Topology**: a task can admit several, equally optimal *routes* (distinct
  reaction sets -- see `mteapy.enumeration`/`mteapy.taskdb`). Which one
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
   recording *which* complex achieved the max (or why none did -- see
   `ReactionEvidence` below), which `map_gpr` discards entirely.
3. Aggregate a route's score from its reactions' scores: "min" (default,
   the strict weakest-link choice this project's own route-tie analysis
   validated for the CellFie-style expression case), "median" (a more
   permissive alternative that analysis also tried, but found far too
   tie-insensitive to trust alone), or "mean" (matches the published
   TIDE/TIDE-essential methodology this package also implements directly --
   use this one when comparing against or reproducing TIDE results).
4. A task's winning route(s) for that sample are whichever are tied for
   the top score. Ties are returned explicitly, never silently broken:
   this project found that a naive `max(dict, key=dict.get)` on tied
   scores systematically (and invisibly) favoured whichever route
   happened to be enumerated first -- see `tied_argmax` below.

Every call returns a full report, not just a number: `score_task` gives a
`TaskActivityReport` covering *every* candidate route (not only the
winner), and every reaction within each route keeps its full per-complex
score breakdown and an explicit `ReactionEvidence` classification -- a
score of exactly 0 for a reaction that has no GPR at all (a transporter by
diffusion, a spontaneous reaction) means something categorically different
from a score of 0 for a reaction with real candidate complexes that
happen to show no expression support anywhere, and this report keeps that
distinction visible rather than collapsing both to a bare "0.0".
"""

from __future__ import annotations

import math
import statistics
from dataclasses import dataclass, field
from enum import Enum

import pandas as pd
from cobra.core import Model

from mteapy.complexes import get_enzymes

_SUPPORTED_AGGREGATIONS = ("min", "median", "mean")

Complex = tuple[str, ...]


class ReactionEvidence(Enum):
    """How well-supported one reaction's score is, for one sample."""

    NO_GPR = "no_gpr"
    """The reaction has no gene association at all (diffusion, a spontaneous
    reaction, or simply un-annotated). There is nothing to evaluate -- this
    is not evidence of absence, it's absence of a gene-level control point."""

    NO_EVIDENCE = "no_evidence"
    """The reaction has one or more candidate complexes, but every one of
    them scored exactly 0 against this sample's expression: real genes are
    annotated, but none show any supporting expression."""

    AMBIGUOUS = "ambiguous"
    """More than one candidate complex tied for the (non-zero) top score.
    The reaction's score is well-defined; which specific complex earned it
    is not."""

    SUPPORTED = "supported"
    """Exactly one candidate complex uniquely achieved the (non-zero) top
    score: the cleanest possible case."""


_SUPPORTED_OR_FUNCS = ("max", "absmax")


def tied_argmax(scores: dict, rel_tol: float = 1e-9, abs_tol: float = 1e-9, or_func: str = "max") -> tuple:
    """All keys tied for the max value in `scores`, within floating tolerance.

    Returns an empty tuple for an empty `scores`, a single-element tuple for
    a clear winner, and a multi-element tuple when several keys are
    genuinely tied -- callers must not just take `[0]` and assume it means
    anything on its own; check `len(...)` first.

    `or_func` picks how "max" is judged: "max" (default) compares raw
    values, correct for non-negative signals like expression, where 0 is
    the floor and bigger is always more support. "absmax" compares by
    magnitude instead (`abs(v)`), needed for a signal that can be negative
    (e.g. a log-fold-change, TIDE-style) -- a strongly *down*-regulated
    complex has to be able to win an OR just as much as a strongly
    up-regulated one would, which plain "max" would blind by always
    preferring the least-negative value. This mirrors the `or_func`
    convention `mteapy.utils.map_gpr`/TIDE already use, just made
    tie-aware: `absmax` there silently returns one winner via
    `np.argmax`, which is exactly the kind of hidden tie-breaking this
    project's own `tied_argmax` exists to avoid.
    """
    if not scores:
        return ()
    if or_func not in _SUPPORTED_OR_FUNCS:
        raise ValueError(f"Unsupported or_func {or_func!r}; use one of {_SUPPORTED_OR_FUNCS}")
    key = abs if or_func == "absmax" else (lambda v: v)
    max_key_val = max(key(v) for v in scores.values())
    return tuple(sorted(
        (k for k, v in scores.items() if math.isclose(key(v), max_key_val, rel_tol=rel_tol, abs_tol=abs_tol)),
        key=lambda c: c,
    ))


def build_complex_cache(model: Model, reaction_ids) -> dict[str, tuple[Complex, ...]]:
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


@dataclass
class ReactionAssessment:
    """One reaction's full scoring breakdown for one sample."""

    reaction_id: str
    score: float
    complex_scores: dict[Complex, float] = field(default_factory=dict)
    """Every candidate complex considered, each with its own score -- not
    just the winner, so a near-miss alternative complex stays visible."""
    winning_complexes: tuple[Complex, ...] = ()
    """The tied_argmax winner set among `complex_scores`. Empty when
    `evidence` is NO_GPR or NO_EVIDENCE (there is no meaningful "winner"
    among complexes that are absent or all show zero support) -- even
    though, numerically, "all tied at 0" would otherwise also satisfy
    tied_argmax, that is not reported here as a winning set."""
    evidence: ReactionEvidence = ReactionEvidence.NO_GPR


def assess_reaction(reaction_id: str, complexes: tuple[Complex, ...], gene_dict: dict[str, float],
                     or_func: str = "max") -> ReactionAssessment:
    """Score one reaction's candidate complexes against `gene_dict` and classify the result.

    AND (a complex's own genes) is always `min` of the raw values,
    regardless of `or_func` -- that combination rule is meaningful whether
    `gene_dict` holds non-negative expression or signed log-fold-change
    values, so it never needs to change. `or_func` (see `tied_argmax`)
    only affects OR, i.e. which candidate complex counts as "the" winner
    among a reaction's alternatives: "max" (default) for expression-like
    signals, "absmax" for a signal that can be negative.
    """
    if or_func not in _SUPPORTED_OR_FUNCS:
        raise ValueError(f"Unsupported or_func {or_func!r}; use one of {_SUPPORTED_OR_FUNCS}")
    if not complexes:
        return ReactionAssessment(reaction_id, 0.0, {}, (), ReactionEvidence.NO_GPR)

    complex_scores = {c: min(gene_dict.get(g, 0.0) for g in c) for c in complexes}
    score = max(complex_scores.values(), key=abs) if or_func == "absmax" else max(complex_scores.values())

    if score == 0.0:
        return ReactionAssessment(reaction_id, 0.0, complex_scores, (), ReactionEvidence.NO_EVIDENCE)

    winners = tied_argmax(complex_scores, or_func=or_func)
    evidence = ReactionEvidence.SUPPORTED if len(winners) == 1 else ReactionEvidence.AMBIGUOUS
    return ReactionAssessment(reaction_id, score, complex_scores, winners, evidence)


def _aggregate(values: list[float], aggregation: str) -> float:
    if not values:
        return 0.0
    if aggregation == "min":
        return min(values)
    if aggregation == "median":
        return statistics.median(values)
    if aggregation == "mean":
        return statistics.mean(values)
    raise ValueError(f"Unsupported aggregation {aggregation!r}; use one of {_SUPPORTED_AGGREGATIONS}")


@dataclass
class RouteReport:
    """One route's full scoring breakdown for one sample."""

    route_id: int
    score: float
    aggregation: str
    reactions: dict[str, ReactionAssessment] = field(default_factory=dict)

    @property
    def is_complete(self) -> bool:
        """True if every gene-associated reaction along this route is
        cleanly SUPPORTED. NO_GPR reactions don't count against
        completeness -- there is nothing to evaluate for them, so their
        presence is not a gap in the evidence. AMBIGUOUS and NO_EVIDENCE
        reactions do count against it: they are real gaps in what this
        route's score can be said to be backed by.
        """
        return all(
            r.evidence == ReactionEvidence.SUPPORTED
            for r in self.reactions.values()
            if r.evidence != ReactionEvidence.NO_GPR
        )

    @property
    def incomplete_reactions(self) -> tuple[str, ...]:
        """reaction_ids responsible for `is_complete` being False."""
        return tuple(sorted(
            rid for rid, r in self.reactions.items()
            if r.evidence in (ReactionEvidence.NO_EVIDENCE, ReactionEvidence.AMBIGUOUS)
        ))


@dataclass
class TaskActivityReport:
    """One task's full context-aware activity report for one sample."""

    task_id: str
    score: float
    """The task's global activity score for this sample -- the winning
    route(s)' score (all tied winners share the same score, by definition)."""
    winning_route_ids: tuple[int, ...]
    routes: dict[int, RouteReport] = field(default_factory=dict)
    """Every candidate route's report, not only the winner(s) -- so a
    near-miss alternate route stays inspectable."""

    @property
    def is_tied(self) -> bool:
        """True when more than one route is tied for the top score.

        Route-level ties turned out to be the norm, not the exception in
        this project's own analysis -- always check this before treating
        `winning_route_ids[0]` as *the* answer.
        """
        return len(self.winning_route_ids) > 1

    @property
    def is_complete(self) -> bool:
        """True only if every winning route (all of them, in case of a tie)
        is itself complete. False (partial) if the winning route relies on
        any reaction with ambiguous or absent expression evidence, or if
        the task itself has no winning route at all (e.g. an empty route
        set)."""
        if not self.winning_route_ids:
            return False
        return all(self.routes[rid].is_complete for rid in self.winning_route_ids)


def score_task(
    task_routes: dict[int, frozenset[str]],
    complex_cache: dict[str, tuple[Complex, ...]],
    gene_dict: dict[str, float],
    aggregation: str = "min",
    task_id: str = "",
    or_func: str = "max",
    reaction_cache: dict[str, "ReactionAssessment"] | None = None,
) -> TaskActivityReport:
    """Score every candidate route of one task against one sample's expression.

    Parameters
    ----------
    task_routes:
        `{route_id: reaction_id_set}` for one task, e.g. from
        `mteapy.taskdb.load_task_routes` (or `load_task_list_routes`, keyed
        down to one task). A single-route task (a plain dict of one entry)
        works fine here -- there is simply nothing for the topology
        dimension to distinguish, and `is_tied` will be False.
    complex_cache:
        `{reaction_id: candidate_complexes}`, from `build_complex_cache`.
        Must cover every reaction appearing in `task_routes`.
    gene_dict:
        One sample's signal, `{gene_id: value}` -- non-negative expression
        (CellFie-style; use `or_func="max"`) or signed log-fold-change
        (TIDE-style; use `or_func="absmax"`).
    aggregation:
        How a route's reaction scores combine into one route score: "min"
        (default, strict weakest-link -- this project's own route-tie
        analysis validated this over "median" for the CellFie-style
        expression case, `or_func="max"`), "median" (permissive), or "mean"
        (matches the published TIDE/TIDE-essential methodology this
        package also implements directly -- see `mteapy.tide`'s
        `calculate_TIDE_scores`/`calculate_TIDEe_scores` -- so use "mean"
        when comparing against or reproducing TIDE results with
        `or_func="absmax"`, rather than carrying the expression-case "min"
        default over to the differential-expression case unexamined).
    task_id:
        Carried through to the result only for the caller's convenience
        (e.g. building a results table); not used in scoring.
    or_func:
        Passed through to `assess_reaction`/`tied_argmax`: "max" (default)
        for a non-negative signal, "absmax" for one that can be negative.
    reaction_cache:
        Optional `{reaction_id: ReactionAssessment}` cache, reused (and
        filled in) across calls. `assess_reaction`'s result for a given
        reaction depends only on that reaction's (fixed) candidate
        complexes and `gene_dict` -- never on which route or task is
        asking -- so a caller scoring many tasks/routes against the *same*
        `gene_dict` (e.g. every task for one sample, or one permutation's
        shuffled values across the whole task list) should build one cache
        and pass it to every call: routes of the same task typically share
        90%+ of their reactions, and reactions repeat across tasks too, so
        without this a task list can trigger orders of magnitude more
        `assess_reaction` calls than there are distinct reactions. Passing
        None (default) scopes the cache to just this one call, matching
        the previous (uncached-across-calls) behavior.
    """
    if aggregation not in _SUPPORTED_AGGREGATIONS:
        raise ValueError(f"Unsupported aggregation {aggregation!r}; use one of {_SUPPORTED_AGGREGATIONS}")

    cache = reaction_cache if reaction_cache is not None else {}

    route_reports: dict[int, RouteReport] = {}
    for route_id, reactions in task_routes.items():
        assessments = {}
        for rid in reactions:
            if rid not in cache:
                cache[rid] = assess_reaction(rid, complex_cache.get(rid, ()), gene_dict, or_func=or_func)
            assessments[rid] = cache[rid]
        agg_score = _aggregate([a.score for a in assessments.values()], aggregation)
        route_reports[route_id] = RouteReport(route_id, agg_score, aggregation, assessments)

    winning_route_ids = tied_argmax({rid: r.score for rid, r in route_reports.items()}, or_func=or_func)
    final_score = route_reports[winning_route_ids[0]].score if winning_route_ids else 0.0

    return TaskActivityReport(
        task_id=task_id,
        score=final_score,
        winning_route_ids=winning_route_ids,
        routes=route_reports,
    )


class RunCancelled(Exception):
    """Raised by a `progress` callback to stop `score_tasks_report` early."""


def score_tasks_report(
    tasks_routes: dict[str, dict[int, frozenset[str]]],
    model: Model,
    expr_df: pd.DataFrame,
    aggregation: str = "min",
    or_func: str = "max",
    progress=None,
    complex_cache: dict | None = None,
) -> dict[str, pd.DataFrame]:
    """Score every task in `tasks_routes` against every sample (column) in `expr_df`.

    Returns four tasks x samples DataFrames (task ids as the index, sample
    names as the columns):

    - ``scores``: each cell's activity score (`TaskActivityReport.score`)
    - ``complete``: bool, `TaskActivityReport.is_complete` -- whether the
      score rests on fully supported evidence rather than ambiguous or
      absent reactions
    - ``tied``: bool, `TaskActivityReport.is_tied` -- more than one route
      shares the top score
    - ``winners``: the tuple of winning route ids (`winning_route_ids`)

    Call `score_task` directly (per task, per sample) for the full
    per-route/per-reaction/per-complex report.

    Parameters
    ----------
    tasks_routes:
        `{task_id: {route_id: reaction_id_set}}` for every task to score,
        e.g. `mteapy.taskdb.load_task_list_routes`.
    model:
        The COBRApy model the routes' reaction ids come from, used only to
        look up each reaction's GPR for `build_complex_cache`.
    expr_df:
        Gene expression, genes (index) x samples (columns).
    aggregation:
        Passed through to `score_task` ("min", "median", or "mean").
    or_func:
        Passed through to `score_task`: "max" (default, for a non-negative
        signal like expression) or "absmax" (for a signal that can be
        negative, like a log-fold-change).
    progress:
        Optional callable ``progress(samples_done, samples_total)``, called
        after each sample; it may raise `RunCancelled` to stop the run.
    complex_cache:
        A prebuilt `build_complex_cache(model, ...)` covering every reaction
        in `tasks_routes`, to skip rebuilding it when scoring the same task
        list repeatedly.
    """
    if complex_cache is None:
        all_reactions = sorted(
            {r for routes in tasks_routes.values() for reactions in routes.values() for r in reactions})
        complex_cache = build_complex_cache(model, all_reactions)

    cells: dict[str, dict[str, dict]] = {"scores": {}, "complete": {}, "tied": {}, "winners": {}}
    for key in cells:
        cells[key] = {task_id: {} for task_id in tasks_routes}

    # Looped sample-first (not task-first) so one reaction_cache can be
    # shared across every task for a given sample -- assess_reaction's
    # result only depends on (reaction, gene_dict), and routes/tasks
    # overlap heavily in which reactions they use, so scoring task-first
    # would rebuild the same (reaction, sample) result redundantly for
    # every task that happens to share it.
    samples = list(expr_df.columns)
    for done, sample in enumerate(samples, 1):
        gene_dict = expr_df[sample].to_dict()
        reaction_cache: dict = {}
        for task_id, routes in tasks_routes.items():
            report = score_task(
                routes, complex_cache, gene_dict, aggregation=aggregation, task_id=task_id, or_func=or_func,
                reaction_cache=reaction_cache,
            )
            cells["scores"][task_id][sample] = report.score
            cells["complete"][task_id][sample] = report.is_complete
            cells["tied"][task_id][sample] = report.is_tied
            cells["winners"][task_id][sample] = tuple(report.winning_route_ids)
        if progress is not None:
            progress(done, len(samples))

    out = {}
    for key, rows in cells.items():
        frame = pd.DataFrame(rows).T
        out[key] = frame.reindex(index=list(tasks_routes), columns=samples)
    return out


def score_tasks_matrix(
    tasks_routes: dict[str, dict[int, frozenset[str]]],
    model: Model,
    expr_df: pd.DataFrame,
    aggregation: str = "min",
    or_func: str = "max",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """`(scores, complete)` of `score_tasks_report` -- the two tasks x samples
    DataFrames a CellFie/TIDE score matrix has, so a caller working
    matrix-first still sees which scores rest on ambiguous or absent
    evidence rather than that information being silently dropped."""
    report = score_tasks_report(tasks_routes, model, expr_df, aggregation=aggregation, or_func=or_func)
    return report["scores"], report["complete"]


def compute_TAS(
        expr_data: pd.DataFrame,
        task_structure: pd.DataFrame | None,
        model: Model,
        aggregation: str = "min",
        or_func: str = "max",
        mapping_strategy: str = "classic",
        tasks_routes: dict | None = None,
    ):
    """Task Activity Score (TAS): project expression straight through each
    task's GPR (AND = min, OR = `or_func`) and aggregate across its
    reaction(s) by `aggregation` -- no percentile-threshold gene-activity
    transform, no essentiality/permutation step. This is CellFie/TIDE's
    peer for the case where that transform is itself the thing being
    second-guessed: a reaction/task scoring zero after thresholding is
    ambiguous between "genuinely off" and "just under this threshold
    convention's floor" (see the `metabolic-variability` project's own
    `01_gtex_reaction_expression_qc.ipynb`/`model_curation.qmd` for a worked
    case), and TAS's score is exactly `TaskActivityReport.score` -- the
    plain, untransformed quantity every other method's own activity claim
    is ultimately built from.

    Parameters
    ----------
    expr_data: pandas.DataFrame
        Gene expression, genes (index) x samples (columns). Non-negative
        expression pairs with `or_func="max"`; a signed signal (e.g.
        log-fold-change) needs `or_func="absmax"`.

    task_structure: pandas.DataFrame
        A boolean matrix, reactions (index) x tasks (columns) -- the
        classic CellFie/TIDE single-fixed-reaction-set representation.
        Required when `mapping_strategy="classic"`; ignored otherwise.

    model: cobra.core.Model
        The COBRA model the task structure's/routes' reaction ids come
        from, used to look up each reaction's GPR.

    aggregation: str ["min" | "median" | "mean"]
        How a task's (or one route's) reaction scores combine into one
        task score (default: "min", the strict weakest-link reading this
        project's own route-tie analysis validated for the expression
        case -- see `score_task`'s own docstring for when "median"/"mean"
        are the better fit instead).

    or_func: str ["max" | "absmax"]
        Passed through to `assess_reaction`/`tied_argmax`: "max" (default)
        for a non-negative signal, "absmax" for one that can be negative.

    mapping_strategy: str ["classic" | "context-aware"]
        "classic" (default) scores each task's single, fixed reaction set
        from `task_structure`. "context-aware" instead scores every
        enumerated alternate route of each task and takes the
        best-supported one for each sample (`score_tasks_matrix`), which
        needs `tasks_routes` instead of `task_structure`.

    tasks_routes: dict
        `{task_id: {route_id: reaction_id_set}}`. Required when
        mapping_strategy="context-aware"; ignored otherwise.

    Returns
    -------
    scores_df: pandas.DataFrame
        Tasks (index, named "task_id") x samples: each cell's TAS.

    complete_df: pandas.DataFrame
        Same shape, boolean: `TaskActivityReport.is_complete` for that
        (task, sample) -- whether the score rests on fully SUPPORTED
        evidence rather than any AMBIGUOUS/NO_EVIDENCE reaction. Unlike
        CellFie's `binary_scores_df`, this is not an activity call against
        a threshold -- TAS makes no such call -- it is purely an evidence-
        completeness flag alongside the raw score.
    """
    if mapping_strategy == "context-aware":
        if tasks_routes is None:
            raise ValueError("mapping_strategy='context-aware' requires tasks_routes.")
    elif mapping_strategy == "classic":
        if task_structure is None:
            raise ValueError("mapping_strategy='classic' requires task_structure.")
        task_structure = task_structure.astype(bool)
        tasks_routes = {
            task_id: {0: frozenset(task_structure.index[task_structure[task_id]])}
            for task_id in task_structure.columns
        }
    else:
        raise ValueError(f"Unsupported mapping_strategy {mapping_strategy!r}. Please, use 'classic' or 'context-aware'.")

    scores_df, complete_df = score_tasks_matrix(tasks_routes, model, expr_data, aggregation=aggregation, or_func=or_func)
    scores_df.index.name = "task_id"
    complete_df.index.name = "task_id"
    return scores_df, complete_df
