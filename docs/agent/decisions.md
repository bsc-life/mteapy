# Decisions and why

Newest first. Each entry is here because the reasoning isn't obvious
from the diff alone -- if you're about to "simplify" or "fix" one of
these, read the entry first.

## TAS (Task Activity Score) added as a third scoring method

`mteapy.context_scoring.score_tasks_matrix` already did raw
expression-through-GPR projection with configurable aggregation
(min/median/mean) -- it just wasn't exposed as a first-class, CLI-facing
method the way CellFie/TIDE are. It was promoted to one (`compute_TAS`,
`run-mtea analyze TAS`) specifically to answer a methodological question
that came up in the GTEx analysis: when a task scores as inactive under
CellFie, is that because the expression signal is genuinely low, or
because CellFie's own percentile-threshold transform is responsible?
TAS has no threshold step at all, so re-scoring with it directly tests
that. (It genuinely resolved one such case in downstream use -- a set of
amino-acid-degradation tasks that turned out to be always-zero for a
reason unrelated to either method: a task-definition bounds issue in
the task list being scored, not a scoring artifact in either method.)

Named after `TaskActivityReport` (the dataclass `score_task` already
returned) rather than inventing new vocabulary -- considered "RTS" (Raw
Task Score) and "TGES" (Task Gene Expression Score) first; TAS won
because it's literally the CLI-facing name for data that already existed
under that name internally.

## `task_sources.task_list` grouping added; `cellfie_consensus` retired

Found via a double-counting bug: `full`, `cellfie_consensus`, and
`cellfie_consensus_gurobi` all coexisted as sources in the routes
database, and callers' own `load_multiroute_tasks`-style helpers
summed across all of them with no source filter -- 220(full) +
145(cellfie_consensus) + 145(cellfie_consensus_gurobi) = 510+
"multi-route tasks", double-counting the same CellFie list twice under
two different solvers and blending incompatible task_id namespaces from
`full` and `cellfie` together.

Fix was two-part: (1) add `task_list` as an explicit grouping column on
`task_sources` so a caller can select "the cellfie task list" by name
without knowing which specific `source` backs it right now, with
`get_sources_for_task_list` deliberately raising rather than silently
picking one if more than one source is ever registered under the same
`task_list` -- ambiguity there means the database's classification needs
fixing, not a silent pick. (2) The user's explicit call on the duplicate
CellFie variants: "we don't want duplicate cellfie sub families, either
consolidate them or keep gurobi ones" -- kept Gurobi, retired the other
via `reset_source` (which deletes routes/route_reactions/
enumeration_runs for that source+model but leaves `tasks`/`models`
alone, in case anything else still references the task rows).

## `register_task_source`'s optional fields use COALESCE, not overwrite

Self-inflicted bug, caught in the same session it was introduced: while
backfilling `task_list` for existing sources, a call to
`register_task_source` without recomputing `origin_repo`/`origin_ref`
(implicit `None`) wiped those fields, because the original SQL used
unconditional `excluded.origin_repo` on conflict instead of
`COALESCE(excluded.origin_repo, task_sources.origin_repo)`. Fixed for
both those fields and the new `task_list` field from the start. **Any
new optional column added to `task_sources` (or similar upsert-pattern
tables) needs the same COALESCE treatment**, or a caller that
legitimately doesn't have a value for it will erase a previously-good one.

## 18 long-stalled `full` tasks got a pFBA-only route, not a deeper fix

These 18 tasks (amino-acid/nucleotide uptake, membrane lipid de novo
synthesis, a few others) hit a hard enumeration timeout at ~9000s with
zero routes recorded. Confirmed via an independent MN5 run (a dedicated
2-hour-per-task retry pass, run separately from this session) that more
time alone doesn't fix it -- this is a genuine combinatorial bottleneck
in the full alternate-route search for these specific tasks, not a
stale-code artifact (an earlier hypothesis -- a known performance fix in
`set_min_total_flux_objective` -- was checked and ruled out: the single
reference-route LP alone solves in under a second for these tasks, the
blow-up is specifically in the cardinality-cut recursion).

Decision: give each one just its single pFBA reference route (the same
solve TIDE/CellFie-style scoring would use regardless), recorded with
`truncated=True` so `is_exhaustive` stays `False` -- a future `--mode
resume` can still search harder. This is explicitly a partial result,
not a claim that no alternates exist. Raised "full" task coverage from
220/257 to 238/257 tasks with >=1 route.

## CellFie's two bugs: fixed after the user set the certainty bar

Two divergences from the original MATLAB CellFie were found by reading
`CellFie.m`/`selectGeneFromGPR_all.m`/`findUsedGenesLevels_all.m`
directly:

1. **Percentile-threshold space**: original MATLAB computes
   `10^prctile(log10(x), p)` (log-space then exponentiate); this
   codebase computed `np.quantile(x, p)` directly (linear space). Fixed
   with a `log_transformed` flag (default `False`, matching MATLAB's
   assumption that input is raw/linear expression) rather than silently
   reinterpreting existing callers' data -- the user's own instruction:
   "for the fix I think the safe [option] would be to have a flag
   indicating whether that expression data is already log-transformed."
2. **Missing-gene propagation**: original MATLAB uses an explicit `-1`
   "no data" sentinel that propagates through AND(min)/OR(max)
   naturally, and excludes a task's all-missing reactions from its score
   average entirely. This codebase defaulted missing genes to `0.0` and
   zero-filled missing reactions into task averages. Fixed to propagate
   `None`/`NaN` instead. The user asked for 100% certainty before this
   one specifically (it's the riskier fix, changing more call sites) --
   confirmed by re-reading the MATLAB source precisely before touching
   `map_gpr_w_names`/`calculate_RAL`/`calculate_CellFie_scores`.

## Git LFS adopted for `routes_human2.db`; full history rewrite, with consent

The database grew large enough (tens of MB, climbing) to need LFS.
`git lfs migrate import --everything` rewrites **every** commit
reachable from **every** ref it finds `.gitattributes` relevant to from
the root commit forward, not just the commits that touch the tracked
file -- this silently rewrote all of `main`'s history (57 commits), not
just the 9 commits on the feature branch being worked on at the time.
Caught before pushing; the user was given the full consequence (new
hashes for everything, tags becoming historically detached, anyone with
an existing clone needing to re-clone) and explicitly chose to proceed
anyway ("I manage this repo and I'm leader and currently only active
developer"). Local tags got silently remapped by the same migration and
had to be force-reset to match the untouched remote tags afterward.

## CellFie model/task import from the original repo: deferred, not decided

The original goal (let a user pick which model "flavor" -- Human-GEM vs.
CellFie's own bundled model -- to run analysis against, importing
CellFie's model/tasks the same way Human-GEM's are handled) stalled on a
real licensing question: CellFie's original repo is GPLv3; mteapy is
MIT. Converting `.mat` -> SBML/`.txt` does not escape GPL copyleft --
copyleft attaches to the derivative content, not the file syntax. The
user asked to pause and think about this rather than proceed ("let me
think about this and let's move back to the scientific part"). **This is
still open**, not abandoned -- don't restart the import work without
first resolving the licensing question with the user.
