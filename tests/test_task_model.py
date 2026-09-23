import pytest
from cobra.core import Metabolite

from mteapy.task_model import TaskModelError, build_task_model, check_task, check_task_list
from mteapy.tasks import BoundedMetabolite, ChangedBound, EquationConstraint, MetabolicTask


def test_equation_handles_charged_species_name_with_plus(toy_model):
    # Real Human-GEM tasks write EQU entries like
    # "ATP[c] + H2O[c] => ADP[c] + Pi[c] + H+[c]" -- the proton's own name
    # contains a "+" with no surrounding space, which must NOT be treated
    # as an additive separator the way "ATP[c] + H2O[c]" is (that one does
    # have spaces on both sides of its "+"). A naive `side.split("+")`
    # shatters "H+[c]" into "H" and "[c]" and fails to resolve either.
    model = toy_model.copy()
    proton = Metabolite("H_c", name="H+", compartment="c")
    model.add_metabolites([proton])
    # No separate sink reaction for the proton: a single-metabolite
    # reaction is itself a boundary reaction, which build_task_model
    # correctly closes as part of setting up the task -- it must instead
    # leave the model via the task's own OUT entry, exactly like a real
    # EQU byproduct would in practice (the real model already has
    # non-boundary reactions that handle H+[c] balance).

    # Block both native A[c]->C[c] routes (R3 directly, R1 as the first
    # step of R1+R2), so the task is only satisfiable through the new EQU
    # pseudo-reaction -- otherwise this test would pass even with EQU
    # parsing completely broken.
    task = MetabolicTask(
        id="1", description="equation with a charged species",
        inputs=[BoundedMetabolite("A[c]", 1, 1)],
        outputs=[BoundedMetabolite("C[c]", 1, 1), BoundedMetabolite("H+[c]", 0, 1000)],
        changed_bounds=[ChangedBound("R1", 0, 0), ChangedBound("R3", 0, 0)],
        equations=[EquationConstraint("A[c] => C[c] + H+[c]", 1, 1)],
    )
    tmodel = build_task_model(model, task)
    solution = tmodel.optimize()
    assert solution.status == "optimal"

    without_equation = MetabolicTask(
        id="2", description="no equation, R1/R3 blocked",
        inputs=[BoundedMetabolite("A[c]", 1, 1)],
        outputs=[BoundedMetabolite("C[c]", 1, 1), BoundedMetabolite("H+[c]", 0, 1000)],
        changed_bounds=[ChangedBound("R1", 0, 0), ChangedBound("R3", 0, 0)],
    )
    tmodel2 = build_task_model(model, without_equation)
    assert tmodel2.optimize().status == "infeasible"


def test_build_task_model_closes_boundary_and_clears_objective(toy_model):
    task = MetabolicTask(
        id="1", description="A to C",
        inputs=[BoundedMetabolite("A[c]", 1, 1)],
        outputs=[BoundedMetabolite("C[c]", 1, 1)],
    )
    tmodel = build_task_model(toy_model, task)

    assert tmodel.reactions.get_by_id("EX_X").bounds == (0, 0)
    assert tmodel.reactions.get_by_id("R_TR").objective_coefficient == 0

    solution = tmodel.optimize()
    assert solution.status == "optimal"

    # original model is untouched
    assert toy_model.reactions.get_by_id("EX_X").bounds == (-1000, 1000)


def test_task_infeasible_without_route(toy_model):
    # E[c] does not exist -> unresolvable reference
    task = MetabolicTask(
        id="1", description="bad",
        inputs=[BoundedMetabolite("A[c]", 1, 1)],
        outputs=[BoundedMetabolite("E[c]", 1, 1)],
    )
    with pytest.raises(TaskModelError):
        build_task_model(toy_model, task)

    assert check_task(toy_model, task) == "inconsistent"


def test_check_task_feasible_and_infeasible(toy_model):
    feasible = MetabolicTask(
        id="1", description="A to C",
        inputs=[BoundedMetabolite("A[c]", 1, 1)],
        outputs=[BoundedMetabolite("C[c]", 1, 1)],
    )
    assert check_task(toy_model, feasible) == "optimal"

    infeasible = MetabolicTask(
        id="2", description="no route",
        inputs=[BoundedMetabolite("D[c]", 1, 1)],
        outputs=[BoundedMetabolite("C[c]", 1, 1)],
    )
    assert check_task(toy_model, infeasible) == "infeasible"


def test_check_task_list_reports_should_fail(toy_model):
    should_fail_task = MetabolicTask(
        id="1", description="expected to fail",
        inputs=[BoundedMetabolite("D[c]", 1, 1)],
        outputs=[BoundedMetabolite("C[c]", 1, 1)],
        should_fail=True,
    )
    df = check_task_list(toy_model, [should_fail_task])
    assert df.loc["1", "status"] == "infeasible"
    assert df.loc["1", "passed"] is True or bool(df.loc["1", "passed"]) is True


def test_equation_constraint_provides_alternate_route(toy_model):
    # Block both native B[c]->C[c] paths (R2) and the direct shortcut (R3);
    # the task should only be feasible because the EQU pseudo-reaction
    # supplies an equivalent B[c] -> C[c] conversion.
    task = MetabolicTask(
        id="1", description="equation provides B->C route",
        inputs=[BoundedMetabolite("A[c]", 1, 1)],
        outputs=[BoundedMetabolite("C[c]", 1, 1)],
        changed_bounds=[ChangedBound("R2", 0, 0), ChangedBound("R3", 0, 0)],
        equations=[EquationConstraint("B[c] => C[c]", 1, 1)],
    )
    assert check_task(toy_model, task) == "optimal"

    without_equation = MetabolicTask(
        id="2", description="no equation, R2/R3 blocked",
        inputs=[BoundedMetabolite("A[c]", 1, 1)],
        outputs=[BoundedMetabolite("C[c]", 1, 1)],
        changed_bounds=[ChangedBound("R2", 0, 0), ChangedBound("R3", 0, 0)],
    )
    assert check_task(toy_model, without_equation) == "infeasible"


def test_changed_bound_restricts_reaction(toy_model):
    # Force R3 shut via CHANGED RXN, leaving only the D[c] output unreachable
    task = MetabolicTask(
        id="1", description="force R4 route closed",
        inputs=[BoundedMetabolite("A[c]", 1, 1)],
        outputs=[BoundedMetabolite("D[c]", 1, 1)],
        changed_bounds=[ChangedBound("R4", 0, 0)],
    )
    assert check_task(toy_model, task) == "infeasible"

    unrestricted = MetabolicTask(
        id="2", description="R4 route open",
        inputs=[BoundedMetabolite("A[c]", 1, 1)],
        outputs=[BoundedMetabolite("D[c]", 1, 1)],
    )
    assert check_task(toy_model, unrestricted) == "optimal"


def test_duplicate_output_metabolite_raises(toy_model):
    # Same metabolite listed twice as an OUT (e.g. with two different bound
    # windows) would otherwise silently produce a detached pseudo-reaction
    # and crash later -- it must instead raise a clear, catchable error.
    task = MetabolicTask(
        id="1", description="duplicate output",
        inputs=[BoundedMetabolite("A[c]", 1, 1)],
        outputs=[
            BoundedMetabolite("C[c]", 0, 1),
            BoundedMetabolite("C[c]", 0, 2),
        ],
    )
    with pytest.raises(TaskModelError):
        build_task_model(toy_model, task)
    assert check_task(toy_model, task) == "inconsistent"


def test_duplicate_input_metabolite_raises(toy_model):
    task = MetabolicTask(
        id="1", description="duplicate input",
        inputs=[
            BoundedMetabolite("A[c]", 1, 1),
            BoundedMetabolite("A[c]", 0, 1),
        ],
        outputs=[BoundedMetabolite("C[c]", 1, 1)],
    )
    with pytest.raises(TaskModelError):
        build_task_model(toy_model, task)


def test_unknown_changed_rxn_raises(toy_model):
    task = MetabolicTask(
        id="1", description="bad changed rxn",
        inputs=[BoundedMetabolite("A[c]", 1, 1)],
        outputs=[BoundedMetabolite("C[c]", 1, 1)],
        changed_bounds=[ChangedBound("NOT_A_REACTION", 0, 0)],
    )
    with pytest.raises(TaskModelError):
        build_task_model(toy_model, task)
