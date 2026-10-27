"""Sanity checks on the data handed to a scoring method.

TIDE's input is a differential-expression result (signed log fold changes),
the others' is expression (non-negative). Feeding one the other's kind runs
without error and gives meaningless scores, so every front end (the GUI
before a run, the CLI before computing) calls these and shows what they say.

Each check returns ``(errors, warnings)`` as lists of messages: an error means
the data cannot be scored at all, a warning means it can but probably should
not be.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# Fold changes of a real comparison are rarely beyond +-20 (log2): a column
# whose upper tail is far past that is more likely counts or TPM.
LARGE_LFC = 20.0


def check_lfc(values: pd.Series, name: str = "log-FC column") -> tuple[list[str], list[str]]:
    """Checks on a column that is supposed to hold log fold changes."""
    errors: list[str] = []
    warnings: list[str] = []
    finite = pd.to_numeric(values, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    if finite.empty:
        return [f"'{name}' has no numeric values"], warnings
    if finite.nunique() == 1:
        return [f"every value of '{name}' is {finite.iloc[0]:g}: there is nothing to score"], warnings

    low, high = float(finite.min()), float(finite.max())
    if low >= 0:
        warnings.append(
            f"'{name}' has no negative values (range {low:g} to {high:g}). Log fold changes go both ways, so "
            f"this looks like expression (or another non-fold-change column), not differential expression: TIDE "
            f"would read every gene as up-regulated. Pick the log-FC column, or use TAS/CellFie for expression."
        )
    elif float(finite.abs().quantile(0.99)) > LARGE_LFC:
        warnings.append(
            f"'{name}' has very large values (99% of |values| are within {float(finite.abs().quantile(0.99)):.3g}, "
            f"max {max(abs(low), high):.3g}); log2 fold changes are rarely beyond +-{LARGE_LFC:g}. Check that it "
            f"is not counts or expression."
        )
    return errors, warnings


def check_pvalues(values: pd.Series, name: str = "p-value column") -> tuple[list[str], list[str]]:
    """Checks on a column used as significance (p-value or adjusted p-value)."""
    finite = pd.to_numeric(values, errors="coerce").dropna()
    if finite.empty:
        return [f"'{name}' has no numeric values"], []
    if float(finite.min()) < 0 or float(finite.max()) > 1:
        return [f"'{name}' has values outside [0, 1] ({finite.min():g} to {finite.max():g}): it is not a "
                f"p-value column"], []
    return [], []


def check_expression(frame: pd.DataFrame, method: str = "this method") -> tuple[list[str], list[str]]:
    """Checks on a genes x samples matrix that is supposed to hold expression."""
    values = frame.apply(pd.to_numeric, errors="coerce")
    negative = [c for c in values.columns if (values[c] < 0).any()]
    if not negative:
        return [], []
    low = float(np.nanmin(values[negative].to_numpy()))
    shown = ", ".join(f"'{c}'" for c in negative[:4]) + (" ..." if len(negative) > 4 else "")
    return [], [
        f"{len(negative)} column(s) ({shown}) have negative values (min {low:g}). Expression is non-negative, so "
        f"this looks like log fold changes (a TIDE input) and {method} will not give meaningful scores on it."
    ]


# A warning that starts with this marker is one the front ends must make hard to miss
# (large, coloured); the marker is stripped before showing it.
ALERT_PREFIX = "ALERT: "
ALERT_FRACTION = 0.10       # this share of the genes the routes need missing from the data raises an alert


def check_gene_coverage(present, needed, missing_genes: str = "exclude") -> tuple[list[str], list[str]]:
    """How many of the genes the scoring needs (those in the GPRs of the task routes) the data has.

    The data lacking a gene is not the same as measuring it at 0: with ``missing_genes="exclude"``
    the reactions that depend on it are left out of the task scores, so a low coverage makes the
    scores rest on fewer reactions. Up to `ALERT_FRACTION` missing is a plain warning; more than that
    (or none of them present) is an alert, marked with `ALERT_PREFIX`.
    """
    needed = set(needed)
    if not needed:
        return [], []
    missing = len(needed - set(present))
    if not missing:
        return [], []
    fraction = missing / len(needed)
    effect = ("the reactions that depend on them are left out of the task scores"
              if missing_genes == "exclude" else "they are counted as 0")
    message = f"{missing} of the {len(needed)} genes the task routes need ({fraction:.0%}) are missing from the data; {effect}."
    if fraction < ALERT_FRACTION:
        return [], [message]
    if fraction >= 0.5:
        message += " Most of the model's metabolism cannot be scored with this data: check that the gene ids match the model's."
    else:
        message += " Check the gene ids, and read the scores with that in mind."
    return [], [ALERT_PREFIX + message]
