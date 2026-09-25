"""Enumerate alternate optimal 'routes' for a metabolic task.

Some tasks admit more than one, equally optimal way to satisfy their
IN/OUT/EQU constraints -- alternate routes that can correspond to
genuinely different functional ways to perform the task (e.g. two
parallel pathways that are interchangeable). This uses the recursive MILP
integer-cut method of Lee, Phalakornkule, Domach & Grossmann (2000),
"Recursive MILP model for finding all the alternate optima in LP models
for metabolic networks":

1. Find a *reference* route: minimize the sum of absolute fluxes over the
   model's own (non-pseudo) reactions (see
   `task_model.set_min_total_flux_objective`), subject to the task's
   pseudo source/demand/equation reactions (and any `CHANGED RXN`
   overrides) as hard bounds. This is deliberately **not** the same as
   minimizing the pseudo reactions themselves -- those are usually fixed
   or already tightly bounded by the task, so minimizing them adds little
   to no information and does not actually minimize how much of the real
   network the task recruits. The reference route's support is the set of
   real reactions carrying nonzero flux in this solution. This reference
   route is exactly the single flux distribution a plain pFBA-based
   scoring approach (e.g. CellFie, TIDE) would use -- everything else this
   function finds is invisible to that approach.
2. Add a binary y_r for each reaction r in the current support, linked so
   that y_r = 0 forces v_r = 0 (`lb_r*y_r <= v_r <= ub_r*y_r`).
3. Add one cardinality cut over that support: `sum(y_r) <= |support| - 1`
   -- i.e. at least one currently-active reaction must switch off.
4. Re-solve the same minimize-total-flux problem. If infeasible, no
   further alternative exists and enumeration stops. If feasible, the new
   solution's actual support is the next route; recurse from step 2 using
   it (reusing any y_r already introduced for a reaction that reappears).

Because each cut only ever concerns reactions that have *actually*
appeared in some previously found route (not a pre-computed superset of
"everything that could possibly be active", the way flux variability
analysis would give), the number of binaries stays tied to how many
distinct reactions really get used across the found routes, rather than
to the model's size -- which is what keeps this tractable on a
genome-scale model. Demonstrated on Human-GEM task 5 (ATP regeneration
from glucose): found 10 distinct, biologically meaningful routes (e.g.
swapping ATP- for dATP-dependent hexokinase/pyruvate kinase) in ~9
minutes.

This is the "build" side of context-aware task scoring: slow and
solver-dependent, meant to be run occasionally (once per model/task-list
version) to populate a routes database that `mteapy.routes` then reads
from at scoring time -- not something the typical scoring workflow runs
live.
"""

from __future__ import annotations

from typing import Iterator

import pandas as pd
from cobra.core import Model
from cobra.core.solution import get_solution

from mteapy.task_model import build_metabolite_lookup, build_task_model, is_pseudo_reaction, set_min_total_flux_objective
from mteapy.tasks import MetabolicTask


def _support(solution, flux_threshold: float) -> frozenset[str]:
    return frozenset(
        rid for rid, flux in solution.fluxes.items()
        if not is_pseudo_reaction(rid) and abs(flux) > flux_threshold
    )


def _solve_or_none(tmodel: Model):
    """Solve `tmodel` and return its Solution, or None if infeasible.

    Uses `slim_optimize` to check feasibility before ever asking for
    primal values: `Model.optimize()`/`get_solution` fetch primal values
    unconditionally, and while GLPK tolerates that on an infeasible
    problem, CPLEX's backend raises a hard `CplexSolverError` instead
    (there is genuinely no solution vector to report), which would
    otherwise crash enumeration instead of just ending it.
    """
    tmodel.slim_optimize()
    if tmodel.solver.status != "optimal":
        return None
    return get_solution(tmodel)


def iter_alternate_routes(
    model: Model,
    task: MetabolicTask,
    max_routes: int = 10,
    flux_threshold: float = 1e-7,
    met_lookup: dict[str, str] | None = None,
    solver_threads: int | None = None,
) -> Iterator[frozenset[str]]:
    """Lazily yield up to `max_routes` distinct minimal-total-flux routes for `task`.

    Each yielded `frozenset[str]` is a reaction-id support: the first is
    the reference (pFBA) route (see module docstring), and each subsequent
    one is guaranteed to omit at least one reaction used by the route right
    before it (and, since cuts accumulate, effectively differs from every
    earlier route). Yields nothing if the task itself is infeasible.

    This is a generator specifically so a caller can stop consuming it
    early (e.g. because it is tracking its own wall-clock budget) without
    losing the routes already found -- unlike `enumerate_alternate_routes`,
    which only ever returns once the whole search has ended one way or
    another.

    `solver_threads`, if given, caps how many threads the underlying
    solver may use for *this* task's MILP -- important for CPLEX/Gurobi
    specifically, which otherwise default to auto-detecting and using
    every core on the node. `model.copy()` (used internally to build the
    per-task model) does not preserve a thread-count set on `model` itself,
    so this has to be applied fresh per task; it is a no-op for solvers
    (e.g. GLPK) with no such native parameter.
    """
    lookup = met_lookup if met_lookup is not None else build_metabolite_lookup(model)
    tmodel = build_task_model(model, task, lookup)
    if solver_threads is not None:
        try:
            tmodel.solver.problem.parameters.threads.set(solver_threads)
        except AttributeError:
            pass  # solver has no native thread-count knob (e.g. GLPK)
    set_min_total_flux_objective(tmodel)

    solution = _solve_or_none(tmodel)
    if solution is None:
        return

    support = _support(solution, flux_threshold)
    yield support

    problem = tmodel.problem
    y_vars: dict[str, object] = {}

    for i in range(max_routes - 1):
        new_ids = [rid for rid in support if rid not in y_vars]
        new_vars = [problem.Variable(f"y_{rid}", type="binary") for rid in new_ids]
        for rid, y in zip(new_ids, new_vars):
            y_vars[rid] = y
        tmodel.add_cons_vars(new_vars)

        new_cons = []
        for rid in new_ids:
            reaction = tmodel.reactions.get_by_id(rid)
            y = y_vars[rid]
            new_cons.append(
                problem.Constraint(reaction.flux_expression - reaction.upper_bound * y, ub=0, name=f"mtask_link_ub_{rid}")
            )
            new_cons.append(
                problem.Constraint(reaction.flux_expression - reaction.lower_bound * y, lb=0, name=f"mtask_link_lb_{rid}")
            )
        tmodel.add_cons_vars(new_cons)

        cut_expression = sum(y_vars[rid] for rid in support)
        tmodel.add_cons_vars([problem.Constraint(cut_expression, ub=len(support) - 1, name=f"mtask_cut_{i}")])

        solution = _solve_or_none(tmodel)
        if solution is None:
            return
        support = _support(solution, flux_threshold)
        yield support


def enumerate_alternate_routes(
    model: Model,
    task: MetabolicTask,
    max_routes: int = 10,
    flux_threshold: float = 1e-7,
    met_lookup: dict[str, str] | None = None,
    solver_threads: int | None = None,
) -> list[frozenset[str]]:
    """Enumerate up to `max_routes` distinct minimal-total-flux routes for `task`.

    Returns a list of `frozenset[str]` reaction-id supports (see
    `iter_alternate_routes`). Returns an empty list if the task itself is
    infeasible.
    """
    return list(iter_alternate_routes(model, task, max_routes, flux_threshold, met_lookup, solver_threads))


def compute_task_alternate_routes(
    model: Model,
    tasks: list[MetabolicTask],
    max_routes: int = 10,
    flux_threshold: float = 1e-7,
    verbose: bool = True,
) -> tuple[dict[str, list[frozenset[str]]], pd.DataFrame]:
    """Enumerate alternate routes for every task in `tasks`.

    Returns
    -------
    routes_by_task:
        A dict mapping task id -> list of routes (each a frozenset of
        reaction ids), for every feasible task.
    summary_df:
        A DataFrame indexed by task id with columns ``status``,
        ``n_routes`` and ``included`` (False for infeasible/errored tasks).
    """
    lookup = build_metabolite_lookup(model)
    routes_by_task: dict[str, list[frozenset[str]]] = {}
    summary_rows = []

    for task in tasks:
        if verbose:
            print(f"- Enumerating routes for task {task.id} ({task.description})", end=" ")
        try:
            routes = enumerate_alternate_routes(
                model, task, max_routes=max_routes, flux_threshold=flux_threshold, met_lookup=lookup,
            )
        except Exception as exc:  # noqa: BLE001
            if verbose:
                print(f"ERROR: {exc}")
            summary_rows.append((task.id, "error", 0, False))
            continue

        if not routes:
            if verbose:
                print("skipped (infeasible)")
            summary_rows.append((task.id, "infeasible", 0, False))
            continue

        routes_by_task[task.id] = routes
        summary_rows.append((task.id, "optimal", len(routes), True))
        if verbose:
            print(f"OK ({len(routes)} route(s), sizes {[len(r) for r in routes]})")

    summary_df = pd.DataFrame(summary_rows, columns=["id", "status", "n_routes", "included"]).set_index("id")
    return routes_by_task, summary_df
