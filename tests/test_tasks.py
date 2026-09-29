import os
import textwrap

import pytest

from mteapy.tasks import (
    BoundedMetabolite,
    ChangedBound,
    EquationConstraint,
    MetabolicTask,
    TaskParseError,
    parse_task_file,
    task_definition_hash,
)


def _task(**overrides):
    defaults = dict(id="T1", description="a task", inputs=[], outputs=[], equations=[], changed_bounds=[])
    defaults.update(overrides)
    return MetabolicTask(**defaults)


def test_task_definition_hash_ignores_cosmetic_fields():
    a = _task(description="original description")
    b = _task(description="a totally different description", comments="also different")
    assert task_definition_hash(a) == task_definition_hash(b)


def test_task_definition_hash_changes_with_inputs():
    a = _task(inputs=[BoundedMetabolite("glucose[e]", 0, 1000)])
    b = _task(inputs=[BoundedMetabolite("glucose[e]", 0, 500)])  # different upper bound
    assert task_definition_hash(a) != task_definition_hash(b)


def test_task_definition_hash_changes_with_outputs_equations_and_changed_bounds():
    base = _task()
    with_output = _task(outputs=[BoundedMetabolite("CO2[e]", 0, 1000)])
    with_equation = _task(equations=[EquationConstraint("A[c] => B[c]", 0, 1000)])
    with_changed_bound = _task(changed_bounds=[ChangedBound("MAR00001", 0, 500)])

    hashes = {task_definition_hash(t) for t in (base, with_output, with_equation, with_changed_bound)}
    assert len(hashes) == 4  # all four are genuinely distinct


def test_task_definition_hash_is_order_independent_within_a_list():
    a = _task(inputs=[BoundedMetabolite("A[c]", 0, 1), BoundedMetabolite("B[c]", 0, 1)])
    b = _task(inputs=[BoundedMetabolite("B[c]", 0, 1), BoundedMetabolite("A[c]", 0, 1)])
    # json.dumps(sort_keys=True) only sorts dict keys, not list order -- this
    # documents that IN/OUT order currently DOES matter to the hash (a
    # reordering-only edit to a task file would count as "changed"), which
    # is a conservative, safe default: it can only cause an unnecessary
    # re-enumeration, never a silently-stale resume.
    assert task_definition_hash(a) != task_definition_hash(b)


def _write(tmp_path, name, content):
    path = tmp_path / name
    path.write_text(textwrap.dedent(content))
    return str(path)


def test_parse_simple_in_out_continuation_rows(tmp_path):
    content = """\
ID\tDESCRIPTION\tSHOULD FAIL\tIN\tIN LB\tIN UB\tOUT\tOUT LB\tOUT UB
1\tTask one\t\tglucose[c]\t1\t1\tATP[c]\t30\t32
\t\t\tO2[c]\t6\t6\tCO2[c]\t6\t6
2\tTask two\tTRUE\tH2O[e]\t\t\tO2[e]\t1\t
"""
    path = _write(tmp_path, "tasks.tsv", content)
    tasks = parse_task_file(path)

    assert [t.id for t in tasks] == ["1", "2"]

    t1 = tasks[0]
    assert t1.description == "Task one"
    assert t1.should_fail is False
    assert [bm.metabolite for bm in t1.inputs] == ["glucose[c]", "O2[c]"]
    assert t1.inputs[0].lower_bound == 1 and t1.inputs[0].upper_bound == 1
    assert t1.inputs[1].lower_bound == 6 and t1.inputs[1].upper_bound == 6
    assert [bm.metabolite for bm in t1.outputs] == ["ATP[c]", "CO2[c]"]
    assert t1.outputs[0].lower_bound == 30 and t1.outputs[0].upper_bound == 32

    t2 = tasks[1]
    assert t2.should_fail is True
    # missing IN LB/UB fall back to defaults
    assert t2.inputs[0].lower_bound == 0
    assert t2.inputs[0].upper_bound == 1000
    # OUT UB missing falls back to default, OUT LB given as 1
    assert t2.outputs[0].lower_bound == 1
    assert t2.outputs[0].upper_bound == 1000


def test_parse_semicolon_separated_single_cell(tmp_path):
    content = """\
ID\tDESCRIPTION\tSHOULD FAIL\tIN\tIN LB\tIN UB\tOUT\tOUT LB\tOUT UB\tEQU\tEQU LB\tEQU UB
1\tAerobic ATP\t\tO2[e];glucose[e]\t\t\tH2O[e];CO2[e]\t\t\tATP[c] + H2O[c] => ADP[c] + Pi[c]\t1\t
"""
    path = _write(tmp_path, "tasks.tsv", content)
    tasks = parse_task_file(path)
    assert len(tasks) == 1
    t = tasks[0]
    assert [bm.metabolite for bm in t.inputs] == ["O2[e]", "glucose[e]"]
    assert [bm.metabolite for bm in t.outputs] == ["H2O[e]", "CO2[e]"]
    assert len(t.equations) == 1
    eq = t.equations[0]
    assert eq.equation == "ATP[c] + H2O[c] => ADP[c] + Pi[c]"
    assert eq.lower_bound == 1
    assert eq.upper_bound == 1000  # default UB when only LB given


def test_reversible_equation_default_bounds(tmp_path):
    content = """\
ID\tDESCRIPTION\tEQU\tEQU LB\tEQU UB
1\tReversible\tA[c] <=> B[c]\t\t
"""
    path = _write(tmp_path, "tasks.tsv", content)
    tasks = parse_task_file(path)
    eq = tasks[0].equations[0]
    assert eq.lower_bound == -1000
    assert eq.upper_bound == 1000


def test_changed_rxn_requires_bounds(tmp_path):
    content = """\
ID\tDESCRIPTION\tCHANGED RXN\tCHANGED LB\tCHANGED UB
1\tChange bounds\tR1\t\t
"""
    path = _write(tmp_path, "tasks.tsv", content)
    with pytest.raises(TaskParseError):
        parse_task_file(path)


def test_changed_rxn_parsed(tmp_path):
    content = """\
ID\tDESCRIPTION\tCHANGED RXN\tCHANGED LB\tCHANGED UB
1\tChange bounds\tR1;R2\t0;0\t5;10
"""
    path = _write(tmp_path, "tasks.tsv", content)
    tasks = parse_task_file(path)
    cb = tasks[0].changed_bounds
    assert [c.reaction_id for c in cb] == ["R1", "R2"]
    assert cb[0].lower_bound == 0 and cb[0].upper_bound == 5
    assert cb[1].lower_bound == 0 and cb[1].upper_bound == 10


def test_extra_columns_become_annotations(tmp_path):
    content = """\
ID\tDESCRIPTION\tSYSTEM\tSUBSYSTEM\tREFERENCE
1\tTask\tENERGY\tGLYCOLYSIS\tSome paper
"""
    path = _write(tmp_path, "tasks.tsv", content)
    tasks = parse_task_file(path)
    t = tasks[0]
    assert t.system == "ENERGY"
    assert t.subsystem == "GLYCOLYSIS"
    assert t.annotations == {"REFERENCE": "Some paper"}


def test_missing_id_column_raises(tmp_path):
    content = "DESCRIPTION\tIN\nfoo\tbar\n"
    path = _write(tmp_path, "tasks.tsv", content)
    with pytest.raises(TaskParseError):
        parse_task_file(path)


def test_hash_prefixed_row_is_skipped_entirely(tmp_path):
    # Human-GEM's own task list convention: a "#" in the leading (normally
    # blank) index column disables a task while documenting why, instead of
    # deleting the row -- e.g. tasks 37/38 are disabled this way with the
    # note "The metabolite NA[c] in this task is conflicting with HumanGEM".
    content = """\
\tID\tDESCRIPTION\tIN\tIN LB\tIN UB\tOUT\tOUT LB\tOUT UB
#\t1\tDisabled task\tNA[c]\t1\t1\tfoo[c]\t1\t1\tconflicts with model
\t2\tNormal task\tglucose[c]\t1\t1\tATP[c]\t1\t1
"""
    path = _write(tmp_path, "tasks.tsv", content)
    tasks = parse_task_file(path)
    assert [t.id for t in tasks] == ["2"]


def test_hash_prefixed_row_continuation_rows_also_skipped(tmp_path):
    content = """\
\tID\tDESCRIPTION\tIN\tIN LB\tIN UB\tOUT\tOUT LB\tOUT UB
#\t1\tDisabled task\tNA[c]\t1\t1\tfoo[c]\t1\t1
\t\t\tO2[c]\t1\t1\tCO2[c]\t1\t1
\t2\tNormal task\tglucose[c]\t1\t1\tATP[c]\t1\t1
"""
    path = _write(tmp_path, "tasks.tsv", content)
    tasks = parse_task_file(path)
    assert [t.id for t in tasks] == ["2"]
    assert len(tasks[0].inputs) == 1  # only its own IN, none leaked from the disabled task above


def test_hash_marker_ignored_without_a_leading_index_column(tmp_path):
    # If the file has no leading blank column before ID (offset stays 0),
    # a literal "#" could only ever land inside a real data cell -- there
    # is no separate "index" slot for it to disable a row from, so nothing
    # should be treated as disabled.
    content = """\
ID\tDESCRIPTION\tIN\tIN LB\tIN UB\tOUT\tOUT LB\tOUT UB
1\tTask\tglucose[c]\t1\t1\tATP[c]\t1\t1
"""
    path = _write(tmp_path, "tasks.tsv", content)
    tasks = parse_task_file(path)
    assert [t.id for t in tasks] == ["1"]
