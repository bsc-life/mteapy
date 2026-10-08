"""Declarative descriptions of the scoring methods (name, parameters, defaults).

One definition per method, shared by every front end: the CLI builds its
argparse options from it (`add_arguments`) and the GUI renders its parameter
form from it (`describe_methods`), so defaults, choices and help text can't
drift apart -- `tests/test_methods.py` checks the CLI against it.

Only methods that take an expression matrix and score it against enumerated
routes are described here (TAS for now; CellFie next). TIDE needs a
differential-expression table instead and is not part of this registry yet.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Param:
    name: str                       # also the CLI option (``--name``) and the JSON key
    label: str
    type: str                       # "choice" | "float" | "int" | "bool"
    default: object
    help: str
    choices: tuple = ()
    min: float | None = None
    max: float | None = None
    dest: str | None = None         # argparse dest when it differs from `name` (the CLI's historical attribute)
    when: tuple = ()                # ((param, value), ...): only relevant while all of these hold; the GUI hides it otherwise


@dataclass(frozen=True)
class Method:
    id: str
    label: str
    description: str
    params: tuple[Param, ...]


TAS = Method(
    id="TAS",
    label="TAS - Task Activity Score",
    description=(
        "Projects expression straight through each task's GPR rules (AND = min, OR = the chosen function) and "
        "aggregates a route's reaction scores into one score; a task's score is its best-supported route's. No "
        "percentile threshold and no permutation test -- the plain quantity CellFie's and TIDE's own activity "
        "scores are built from."
    ),
    params=(
        Param("aggregation", "Aggregation", "choice", "min",
              "How a task's (or one route's) reaction scores combine into one task score: min is the strict "
              "weakest link, median is permissive, mean matches the TIDE methodology.",
              choices=("min", "median", "mean")),
        Param("or_func", "OR function", "choice", "max",
              "How OR relationships in gene-protein-reaction rules are resolved: max returns the maximum value "
              "(non-negative expression), absmax the absolute maximum (signed signals such as log-fold-changes).",
              choices=("max", "absmax")),
    ),
)

CELLFIE = Method(
    id="CellFie",
    label="CellFie",
    description=(
        "CellFie (Richelle et al., 2021): converts expression into gene activity levels, 5*ln(1 + expr/threshold), "
        "with a per-gene (local) or one global threshold, projects them through the GPR rules (AND = min, OR = max) "
        "and scores each task as the mean over its best-supported route. The thresholds are computed over the whole "
        "loaded dataset, so the score of a sample depends on the other samples."
    ),
    params=(
        Param("threshold_type", "Threshold", "choice", "local",
              "A local approach uses a different threshold for each gene, a global one the same threshold for all "
              "genes when computing gene activity levels.", choices=("local", "global"), dest="thresh_type"),
        Param("local_threshold_type", "Local threshold", "choice", "minmaxmean",
              "minmaxmean: a gene's threshold is its mean expression across samples, kept between a lower and an "
              "upper bound. mean: the gene's mean expression.", choices=("minmaxmean", "mean"),
              dest="local_thresh_type", when=(("threshold_type", "local"),)),
        Param("minmaxmean_threshold_type", "Bounds are", "choice", "percentile",
              "Whether the lower and upper bounds are percentiles of the distribution of all expression values or "
              "plain values.", choices=("percentile", "value"), dest="minmaxmean_thresh_type",
              when=(("threshold_type", "local"), ("local_threshold_type", "minmaxmean"))),
        Param("upper_bound", "Upper bound", "float", 0.75,
              "Upper bound of the minmaxmean threshold. Percentiles are between 0 and 1.", min=0,
              when=(("threshold_type", "local"), ("local_threshold_type", "minmaxmean"))),
        Param("lower_bound", "Lower bound", "float", 0.25,
              "Lower bound of the minmaxmean threshold. Percentiles are between 0 and 1.", min=0,
              when=(("threshold_type", "local"), ("local_threshold_type", "minmaxmean"))),
        Param("global_threshold_type", "Global threshold is", "choice", "percentile",
              "Whether the global threshold is a percentile of the distribution of all expression values or a "
              "plain value.", choices=("percentile", "value"), dest="global_thresh_type",
              when=(("threshold_type", "global"),)),
        Param("global_value", "Global value", "float", 0.75,
              "The global threshold; percentiles are between 0 and 1.", min=0,
              when=(("threshold_type", "global"),)),
        Param("log_transformed", "Input is log-transformed", "bool", False,
              "Tick when the expression is already log-transformed: percentile thresholds are then taken directly "
              "on the given values instead of in log10 space (the original algorithm assumes raw, linear TPM/FPKM)."),
    ),
)

METHODS: dict[str, Method] = {m.id: m for m in (TAS, CELLFIE)}


def get_method(method_id: str) -> Method:
    try:
        return METHODS[method_id]
    except KeyError:
        raise KeyError(f"unknown method {method_id!r}; available: {list(METHODS)}") from None


def describe_methods() -> list[dict]:
    """JSON-serializable form of every method, for a front end to render."""
    return [
        {"id": m.id, "label": m.label, "description": m.description,
         "params": [{"name": p.name, "label": p.label, "type": p.type, "default": p.default, "help": p.help,
                     "choices": list(p.choices), "min": p.min, "max": p.max,
                     "when": [list(c) for c in p.when]} for p in m.params]}
        for m in METHODS.values()
    ]


def _coerce(param: Param, value):
    if param.type == "choice":
        if value not in param.choices:
            raise ValueError(f"{param.name}: {value!r} is not one of {list(param.choices)}")
        return value
    if param.type == "bool":
        if not isinstance(value, bool):
            raise ValueError(f"{param.name}: expected true or false, got {value!r}")
        return value
    try:
        number = int(value) if param.type == "int" and float(value) == int(float(value)) else float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{param.name}: expected a number, got {value!r}") from None
    if param.type == "int" and not isinstance(number, int):
        raise ValueError(f"{param.name}: expected a whole number, got {value!r}")
    if (param.min is not None and number < param.min) or (param.max is not None and number > param.max):
        raise ValueError(f"{param.name}: {number} is outside [{param.min}, {param.max}]")
    return number


def validate_params(method: Method, params: dict | None) -> dict:
    """`params` with defaults filled in and values checked; unknown names are an error."""
    params = dict(params or {})
    known = {p.name: p for p in method.params}
    unknown = sorted(set(params) - set(known))
    if unknown:
        raise ValueError(f"unknown parameter(s) for {method.id}: {unknown}; known: {list(known)}")
    return {name: _coerce(p, params[name]) if name in params else p.default for name, p in known.items()}


def add_arguments(parser, method: Method, group=None) -> None:
    """Add one ``--<name>`` option per parameter, with the spec's default, choices and help."""
    target = group or parser
    for p in method.params:
        kwargs = {"dest": p.dest or p.name, "default": p.default, "help": p.help}
        if p.type == "choice":
            kwargs.update(type=str, choices=list(p.choices))
        elif p.type == "bool":
            kwargs.update(action="store_true")
        else:
            kwargs.update(type=int if p.type == "int" else float)
        target.add_argument(f"--{p.name}", **kwargs)
