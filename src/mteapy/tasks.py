"""Parsing and data model for metabolic task lists.

Metabolic task lists follow the tab-separated format used across the
Human-GEM / RAVEN ecosystem (see ``checkTasks.m`` / ``parseTaskList.m`` in
RAVEN): one row per task, with optional continuation rows for tasks that
need more than one input/output/equation/changed-bound entry. A cell can
also hold several entries separated by ``;``.

Recognized columns (all optional except ``ID``):

``ID``
    Task identifier. Kept as a string since task lists in the wild use
    both numeric ids ("1", "2", ...) and short codes ("ER", ...).
``DESCRIPTION``
    Free-text description of the task.
``SHOULD FAIL``
    Whether the task is expected to be infeasible (truthy string -> True).
``SYSTEM`` / ``SUBSYSTEM``
    Optional metabolic system/subsystem annotation.
``IN`` / ``IN LB`` / ``IN UB``
    Metabolites that must be taken up (as ``name[compartment]``), with
    bounds on how much can be consumed.
``OUT`` / ``OUT LB`` / ``OUT UB``
    Metabolites that must be produced, with bounds on how much can be
    produced.
``EQU`` / ``EQU LB`` / ``EQU UB``
    A stoichiometric equation (e.g. ``ATP[c] + H2O[c] => ADP[c] + Pi[c]``)
    that must proceed within the given net-flux bounds.
``CHANGED RXN`` / ``CHANGED LB`` / ``CHANGED UB``
    Existing model reactions whose bounds are temporarily overridden for
    the duration of the task.
``PRINT FLUX``
    Whether the task's flux distribution should be reported (informational
    only, not used to build task models).
``COMMENTS``
    Free-text comments.

Any other column present in the file is kept verbatim in
``MetabolicTask.annotations``.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field

DEFAULT_LOWER_BOUND = 0.0
DEFAULT_UPPER_BOUND = 1000.0
DEFAULT_EQUATION_UPPER_BOUND = 1000.0

# Columns that hold repeatable entries, each paired with its LB/UB columns.
_LIST_COLUMNS = {
    "IN": ("IN LB", "IN UB"),
    "OUT": ("OUT LB", "OUT UB"),
    "EQU": ("EQU LB", "EQU UB"),
    "CHANGED RXN": ("CHANGED LB", "CHANGED UB"),
}

_SCALAR_COLUMNS = {"ID", "DESCRIPTION", "SHOULD FAIL", "SYSTEM", "SUBSYSTEM", "PRINT FLUX", "COMMENTS"}

_TRUTHY = {"1", "true", "yes", "y", "t"}


class TaskParseError(ValueError):
    """Raised when a task list file cannot be parsed."""


@dataclass(frozen=True)
class BoundedMetabolite:
    """A metabolite (as ``name[compartment]``) with lower/upper bounds."""

    metabolite: str
    lower_bound: float
    upper_bound: float


@dataclass(frozen=True)
class EquationConstraint:
    """A stoichiometric equation with bounds on its net flux."""

    equation: str
    lower_bound: float
    upper_bound: float


@dataclass(frozen=True)
class ChangedBound:
    """An override of an existing model reaction's bounds."""

    reaction_id: str
    lower_bound: float
    upper_bound: float


@dataclass
class MetabolicTask:
    """A single metabolic task."""

    id: str
    description: str = ""
    should_fail: bool = False
    system: str = ""
    subsystem: str = ""
    inputs: list[BoundedMetabolite] = field(default_factory=list)
    outputs: list[BoundedMetabolite] = field(default_factory=list)
    equations: list[EquationConstraint] = field(default_factory=list)
    changed_bounds: list[ChangedBound] = field(default_factory=list)
    print_flux: bool = False
    comments: str = ""
    annotations: dict[str, str] = field(default_factory=dict)

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return f"MetabolicTask(id={self.id!r}, description={self.description!r})"


def _to_bool(value: str) -> bool:
    return value.strip().lower() in _TRUTHY


def _expand_bounds(raw: str, n: int, default: float) -> list[float]:
    """Expand a (possibly ``;``-separated) bound cell to ``n`` values."""
    raw = raw.strip()
    if not raw:
        return [default] * n
    parts = [p.strip() for p in raw.split(";") if p.strip()]
    if len(parts) == 1:
        return [float(parts[0])] * n
    if len(parts) == n:
        return [float(p) for p in parts]
    raise TaskParseError(f"Expected 1 or {n} bound value(s), got {len(parts)} in {raw!r}")


def _equation_default_bounds(equation: str) -> tuple[float, float]:
    reversible = "<=>" in equation
    lb = -DEFAULT_EQUATION_UPPER_BOUND if reversible else DEFAULT_LOWER_BOUND
    return lb, DEFAULT_EQUATION_UPPER_BOUND


def _new_accumulator() -> dict[str, list]:
    return {
        "IN": [], "IN LB": [], "IN UB": [],
        "OUT": [], "OUT LB": [], "OUT UB": [],
        "EQU": [], "EQU LB": [], "EQU UB": [],
        "CHANGED RXN": [], "CHANGED LB": [], "CHANGED UB": [],
    }


def _accumulate_row(acc: dict[str, list], row: dict[str, str]) -> None:
    for column, (lb_col, ub_col) in _LIST_COLUMNS.items():
        raw = row.get(column, "")
        if not raw.strip():
            continue
        entries = [e.strip() for e in raw.split(";") if e.strip()]
        if column == "EQU":
            # bounds default per-equation depending on reversibility
            lb_values, ub_values = [], []
            for entry in entries:
                default_lb, default_ub = _equation_default_bounds(entry)
                lb_values.extend(_expand_bounds(row.get(lb_col, ""), 1, default_lb))
                ub_values.extend(_expand_bounds(row.get(ub_col, ""), 1, default_ub))
        elif column == "CHANGED RXN":
            lb_raw, ub_raw = row.get(lb_col, ""), row.get(ub_col, "")
            if not lb_raw.strip() or not ub_raw.strip():
                raise TaskParseError(
                    f"CHANGED RXN entry {entries!r} requires both CHANGED LB and CHANGED UB"
                )
            lb_values = _expand_bounds(lb_raw, len(entries), 0.0)
            ub_values = _expand_bounds(ub_raw, len(entries), 0.0)
        else:
            lb_values = _expand_bounds(row.get(lb_col, ""), len(entries), DEFAULT_LOWER_BOUND)
            ub_values = _expand_bounds(row.get(ub_col, ""), len(entries), DEFAULT_UPPER_BOUND)

        acc[column].extend(entries)
        acc[lb_col].extend(lb_values)
        acc[ub_col].extend(ub_values)


def _finalize_task(header_row: dict[str, str], acc: dict[str, list]) -> MetabolicTask:
    task_id = header_row.get("ID", "").strip()
    if not task_id:
        raise TaskParseError("Task is missing an ID")

    annotations = {
        key: value
        for key, value in header_row.items()
        if key not in _SCALAR_COLUMNS and key not in _LIST_COLUMNS and value.strip()
    }

    return MetabolicTask(
        id=task_id,
        description=header_row.get("DESCRIPTION", "").strip(),
        should_fail=_to_bool(header_row.get("SHOULD FAIL", "")),
        system=header_row.get("SYSTEM", "").strip(),
        subsystem=header_row.get("SUBSYSTEM", "").strip(),
        inputs=[
            BoundedMetabolite(m, lb, ub)
            for m, lb, ub in zip(acc["IN"], acc["IN LB"], acc["IN UB"])
        ],
        outputs=[
            BoundedMetabolite(m, lb, ub)
            for m, lb, ub in zip(acc["OUT"], acc["OUT LB"], acc["OUT UB"])
        ],
        equations=[
            EquationConstraint(e, lb, ub)
            for e, lb, ub in zip(acc["EQU"], acc["EQU LB"], acc["EQU UB"])
        ],
        changed_bounds=[
            ChangedBound(r, lb, ub)
            for r, lb, ub in zip(acc["CHANGED RXN"], acc["CHANGED LB"], acc["CHANGED UB"])
        ],
        print_flux=_to_bool(header_row.get("PRINT FLUX", "")),
        comments=header_row.get("COMMENTS", "").strip(),
        annotations=annotations,
    )


def parse_task_file(path: str) -> list[MetabolicTask]:
    """Parse a tab-separated metabolic task list into `MetabolicTask` objects.

    Handles both styles of multi-valued entries found in the wild: several
    ``;``-separated values in a single cell, and/or continuation rows (rows
    with an empty ``ID``) that add more inputs/outputs/equations/changed
    bounds to the task started by the previous non-empty-ID row.

    A task is skipped entirely (not returned) if its leading index column
    (the blank column before ``ID``, when the file has one) holds ``#``
    instead of being empty -- the convention Human-GEM's own task list uses
    to disable a task while documenting why, e.g. a metabolite reference
    that conflicts with the current model, rather than deleting the row.
    Any continuation rows belonging to a disabled task are skipped too.
    """
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise TaskParseError(f"{path} is empty")

        # Some exports have a leading blank/index column before ID.
        if header and header[0].strip() == "":
            header = header[1:]
            offset = 1
        else:
            offset = 0

        if "ID" not in header:
            raise TaskParseError(f"{path}: could not find an ID column in header {header!r}")

        tasks: list[MetabolicTask] = []
        header_row: dict[str, str] | None = None
        acc: dict[str, list] | None = None

        for raw_row in reader:
            is_disabled_row = offset == 1 and len(raw_row) > 0 and raw_row[0].strip().startswith("#")
            row = raw_row[offset:]
            # Pad/truncate defensively in case of ragged rows.
            row = row + [""] * (len(header) - len(row))
            row_dict = dict(zip(header, row))

            if row_dict.get("ID", "").strip():
                if header_row is not None:
                    tasks.append(_finalize_task(header_row, acc))
                if is_disabled_row:
                    header_row = None  # continuation rows below are skipped too, same as a stray blank row
                    continue
                header_row = row_dict
                acc = _new_accumulator()
                _accumulate_row(acc, row_dict)
            else:
                if header_row is None:
                    continue  # stray blank row before the first task, or inside a disabled task
                _accumulate_row(acc, row_dict)

        if header_row is not None:
            tasks.append(_finalize_task(header_row, acc))

    return tasks
