"""Build task-constrained COBRApy models from `MetabolicTask` objects.

To evaluate a task, a copy of the model is made where:

1. All boundary reactions (exchange/sink/demand) are closed.
2. The original objective is cleared.
3. A pseudo "source" reaction is added for each task input, and a pseudo
   "demand" reaction for each task output, bounded as specified by the task.
4. A pseudo reaction is added for each `EQU` constraint, bounded by its
   net-flux range.
5. Existing reaction bounds are overridden for each `CHANGED RXN` entry.
6. The objective is set to (the negative of) the sum of all pseudo
   reactions added in steps 3-4, so that solving the model finds a
   feasible point that consumes/produces as little as possible from the
   task's boundary conditions. This is only meant for feasibility
   checking (see `check_task`/`check_task_list` below) and as a base for
   a pFBA reference-route step; a route enumerator should override this
   objective entirely with its own (minimizing total flux over the
   model's real reactions, not these pseudo ones -- see
   `set_min_total_flux_objective`).

Pseudo reactions are prefixed with ``mtask_`` so they can be filtered out
of downstream flux/essentiality results.

This is a per-task-copy engine: every call to `build_task_model` makes a
full copy of `model`. For gene-essentiality analysis specifically, which
needs to test every gene against every task, a shared-model, bound-toggling
architecture is far faster at genome scale -- this module is for
per-task/per-route analysis, not that.
"""

from __future__ import annotations

import re
from difflib import get_close_matches

from cobra.core import Metabolite, Model, Reaction

from mteapy.tasks import MetabolicTask

SOURCE_PREFIX = "mtask_source_"
DEMAND_PREFIX = "mtask_demand_"
EQUATION_PREFIX = "mtask_equation_"

_MET_PATTERN = re.compile(r"^(?P<name>.*)\[(?P<compartment>[A-Za-z])\]$")
_EQU_ARROW_PATTERN = re.compile(r"<=>|=>")
_EQU_TERM_PATTERN = re.compile(r"^\s*(?:(?P<coef>\d+(?:\.\d+)?)\s+)?(?P<met>\S.*\S|\S)\s*$")
_EQU_TERM_SPLIT_PATTERN = re.compile(r"(?<=\s)\+(?=\s)")


class TaskModelError(ValueError):
    """Raised when a task cannot be translated into a valid model."""


def is_pseudo_reaction(reaction_id: str) -> bool:
    return reaction_id.startswith((SOURCE_PREFIX, DEMAND_PREFIX, EQUATION_PREFIX))


def build_metabolite_lookup(model: Model) -> dict[str, str]:
    """Map ``"name[compartment]"`` -> metabolite id for every model metabolite."""
    return {f"{m.name}[{m.compartment}]": m.id for m in model.metabolites}


def resolve_metabolite(model: Model, lookup: dict[str, str], token: str) -> Metabolite:
    """Resolve a ``"name[compartment]"`` task reference to a model `Metabolite`.

    Raises `TaskModelError` (with the closest-matching names as a hint) if
    the token can't be parsed or doesn't match any model metabolite.
    """
    match = _MET_PATTERN.match(token.strip())
    if not match:
        raise TaskModelError(
            f"Could not parse metabolite reference {token!r} (expected 'name[compartment]')"
        )
    met_id = lookup.get(token.strip())
    if met_id is None:
        suggestions = get_close_matches(token.strip(), lookup.keys(), n=3)
        hint = f" Closest matches: {suggestions}" if suggestions else ""
        raise TaskModelError(f"Metabolite {token!r} not found in model.{hint}")
    return model.metabolites.get_by_id(met_id)


def _add_boundary_reaction(model: Model, metabolite: Metabolite, prefix: str, coefficient: float, lb: float, ub: float) -> Reaction:
    reaction_id = f"{prefix}{metabolite.id}"
    if reaction_id in model.reactions:
        # cobra silently drops (rather than errors on) a reaction whose id already
        # exists, which would otherwise leave this pseudo-reaction detached from
        # the model and crash later when its objective coefficient is set. This
        # happens when the same metabolite is listed twice as a task IN (or twice
        # as an OUT) -- a real data issue seen in curated task lists (e.g. a
        # metabolite repeated with two different bound windows).
        raise TaskModelError(
            f"Duplicate task entry for metabolite {metabolite.id!r} "
            f"(would redefine reaction {reaction_id!r}); "
            "check for a repeated IN/OUT row for this metabolite in the task."
        )
    reaction = Reaction(reaction_id)
    reaction.add_metabolites({metabolite: coefficient})
    reaction.bounds = (lb, ub)
    model.add_reactions([reaction])
    return reaction


def _parse_equation_term(term: str) -> tuple[float, str]:
    match = _EQU_TERM_PATTERN.match(term)
    if not match:
        raise TaskModelError(f"Could not parse equation term {term!r}")
    coef = float(match.group("coef")) if match.group("coef") else 1.0
    return coef, match.group("met").strip()


def _add_equation_reaction(model: Model, lookup: dict[str, str], equation: str, lb: float, ub: float) -> Reaction:
    arrow_match = _EQU_ARROW_PATTERN.search(equation)
    if not arrow_match:
        raise TaskModelError(f"Equation {equation!r} does not contain '=>' or '<=>'")
    lhs, rhs = equation[: arrow_match.start()], equation[arrow_match.end():]

    stoichiometry: dict[Metabolite, float] = {}
    for side, sign in ((lhs, -1.0), (rhs, 1.0)):
        # Split additive terms only on a "+" with whitespace on both sides
        # (e.g. "ATP[c] + H2O[c]"), never on one glued directly to adjacent
        # characters -- a naive `side.split("+")` shatters charged species
        # like "H+[c]" or "NAD+[c]" into "H"/"NAD" and "[c]", since their
        # own name contains a "+" with no surrounding space.
        for term in _EQU_TERM_SPLIT_PATTERN.split(side):
            term = term.strip()
            if not term:
                continue
            coef, token = _parse_equation_term(term)
            metabolite = resolve_metabolite(model, lookup, token)
            stoichiometry[metabolite] = stoichiometry.get(metabolite, 0.0) + sign * coef

    slug = re.sub(r"[^A-Za-z0-9]+", "_", equation).strip("_")[:60]
    reaction_id = f"{EQUATION_PREFIX}{slug}"
    if reaction_id in model.reactions:
        raise TaskModelError(f"Duplicate EQU entry {equation!r} (would redefine reaction {reaction_id!r})")
    reaction = Reaction(reaction_id)
    reaction.add_metabolites(stoichiometry)
    reaction.bounds = (lb, ub)
    model.add_reactions([reaction])
    return reaction


def build_task_model(model: Model, task: MetabolicTask, met_lookup: dict[str, str] | None = None) -> Model:
    """Return a copy of `model` constrained to evaluate `task`.

    Parameters
    ----------
    model:
        The base metabolic model.
    task:
        The task to translate into model constraints.
    met_lookup:
        Optional precomputed `build_metabolite_lookup(model)` result, to
        avoid rebuilding it for every task when processing many tasks.
    """
    tmodel = model.copy()
    lookup = met_lookup if met_lookup is not None else build_metabolite_lookup(model)

    for reaction in tmodel.boundary:
        reaction.bounds = (0.0, 0.0)
    for reaction in tmodel.reactions:
        if reaction.objective_coefficient != 0:
            reaction.objective_coefficient = 0

    pseudo_reactions: list[Reaction] = []

    for bm in task.inputs:
        metabolite = resolve_metabolite(tmodel, lookup, bm.metabolite)
        pseudo_reactions.append(
            _add_boundary_reaction(tmodel, metabolite, SOURCE_PREFIX, 1.0, bm.lower_bound, bm.upper_bound)
        )
    for bm in task.outputs:
        metabolite = resolve_metabolite(tmodel, lookup, bm.metabolite)
        pseudo_reactions.append(
            _add_boundary_reaction(tmodel, metabolite, DEMAND_PREFIX, -1.0, bm.lower_bound, bm.upper_bound)
        )
    for eq in task.equations:
        pseudo_reactions.append(_add_equation_reaction(tmodel, lookup, eq.equation, eq.lower_bound, eq.upper_bound))

    for cb in task.changed_bounds:
        try:
            reaction = tmodel.reactions.get_by_id(cb.reaction_id)
        except KeyError:
            raise TaskModelError(f"CHANGED RXN {cb.reaction_id!r} not found in model")
        reaction.bounds = (cb.lower_bound, cb.upper_bound)

    for reaction in pseudo_reactions:
        reaction.objective_coefficient = -1

    return tmodel


def set_min_total_flux_objective(tmodel: Model) -> None:
    """Point `tmodel`'s objective at minimizing total flux over its real reactions.

    Explicitly excludes the task's own pseudo source/demand/equation
    reactions from the sum: those are usually fixed or already tightly
    bounded by the task itself, so minimizing them (as `build_task_model`'s
    own default objective does, for feasibility-checking purposes) adds
    little to no information and does not actually minimize how much of
    the real network the task recruits. This is what defines a meaningful
    *reference* (pFBA) flux distribution for a task, and the objective a
    route enumerator progressively re-solves under cardinality cuts.
    """
    # Two performance pitfalls confirmed by profiling on a 13000-reaction
    # model, each turning this into an ~8-minute call instead of the ~0.1s
    # it should be:
    # 1. Building the objective as `sum(r.forward_variable + r.reverse_variable
    #    for r in real_reactions)` -- Python's sum() over thousands of sympy
    #    terms -- is catastrophically slow: each `+` builds a new symbolic
    #    expression tree instead of appending to a flat sum.
    #    set_linear_coefficients sets them all in one bulk call instead.
    # 2. Looping `reaction.objective_coefficient = 0` over every reaction to
    #    clear the old objective is itself slow per-call at this scale --
    #    and redundant besides: assigning a whole new Objective to
    #    tmodel.objective already discards every previous coefficient.
    tmodel.objective = tmodel.problem.Objective(0, direction="min")
    coefficients = {}
    for reaction in tmodel.reactions:
        if is_pseudo_reaction(reaction.id):
            continue
        coefficients[reaction.forward_variable] = 1.0
        coefficients[reaction.reverse_variable] = 1.0
    tmodel.solver.objective.set_linear_coefficients(coefficients)


def check_task(model: Model, task: MetabolicTask, met_lookup: dict[str, str] | None = None) -> str:
    """Return the optimization status of `task` on `model` ('optimal', 'infeasible', ...).

    Returns ``'inconsistent'`` if the task itself could not be translated
    into a valid model (e.g. an unresolvable metabolite reference, or a
    duplicated IN/OUT/EQU entry). Returns ``'error'`` for any other
    unexpected failure, so that one malformed task never aborts a batch of
    many.
    """
    try:
        tmodel = build_task_model(model, task, met_lookup)
        solution = tmodel.optimize()
        return solution.status
    except TaskModelError:
        return "inconsistent"
    except Exception:
        return "error"


_worker_model = None
_worker_lookup = None


def _init_check_worker(model, lookup):
    global _worker_model, _worker_lookup
    _worker_model = model
    _worker_lookup = lookup


def _check_worker_task(task):
    return task.id, task.description, task.should_fail, check_task(_worker_model, task, _worker_lookup)


def check_task_list(model: Model, tasks: list[MetabolicTask], verbose: bool = False, processes: int | None = None):
    """Check feasibility of every task in `tasks` against `model`.

    Each task is independent (its own `model.copy()` inside `check_task`),
    so this is embarrassingly parallel across tasks -- pass `processes` to
    split the task list across that many worker processes, each getting
    its own (forked) copy of `model` once, not once per task.

    Returns a pandas DataFrame indexed by task id with columns
    ``description``, ``status`` and ``passed`` (whether the observed
    status matches the task's ``should_fail`` expectation).
    """
    import pandas as pd

    lookup = build_metabolite_lookup(model)

    if processes is None or processes <= 1:
        results = [(task.id, task.description, task.should_fail, check_task(model, task, lookup)) for task in tasks]
    else:
        from multiprocessing import Pool

        with Pool(processes=processes, initializer=_init_check_worker, initargs=(model, lookup)) as pool:
            results = pool.map(_check_worker_task, tasks)

    rows = []
    for task_id, description, should_fail, status in results:
        expected_status = "infeasible" if should_fail else "optimal"
        passed = status == expected_status
        if verbose:
            print(f"[{'OK' if passed else 'FAIL'}] task {task_id} ({description}): {status}")
        rows.append((task_id, description, status, passed))

    return pd.DataFrame(rows, columns=["id", "description", "status", "passed"]).set_index("id")
