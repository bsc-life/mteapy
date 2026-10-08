# Decisions and why

Newest first. Each entry is here because the reasoning isn't obvious
from the diff alone -- if you're about to "simplify" or "fix" one of
these, read the entry first.

## Model folders + manifest registry; method specs; GUI runs

Nothing model-specific is hard-coded any more. A model is a folder with a
`manifest.json` (name, version, model/database files with sha256, task lists,
annotation tables, gene-id type, license) -- `mteapy.registry` finds,
validates and integrity-checks them (modified model file, database built
against a different model, a task list whose hash differs from the one the
database was imported from). A model that fails is *listed* with its
`problems`, not hidden, so the CLI (`models list`) and GUI can say why. This
is for Recon2.2/Recon3D later (the user keeps their own copies precisely to
avoid the CellFie GPL question) -- adding a model should be a folder, not a
code change. Manifest paths must be relative and stay inside the folder.
The remote `index.json` / `models fetch` and forked-repo model sources were
discussed and deliberately deferred; the bundled model stays in the wheel.
The `routes.db` is now at `data/models/HumanGEM/routes.db` (LFS pattern in
`.gitattributes` updated; the old `routes_human2.db` is removed). The
classic/TIDE path's files (`HumanGEM.xml.gz`, `task_structure_matrix.tsv`,
essential genes, `task_metadata.tsv`) stay in `data/` for now. The user
suggested deriving "classic" from each task's first (pFBA) route instead of
the matrix; measured on the CellFie list, route 0 equals the matrix row for
only 102/185 tasks (median similarity 1.0, mean 0.78), so it is a good
definition but would not reproduce legacy classic results -- to be done with
an explicit `route_rank` column (0 = reference), since "first route" is only
implicit (lowest route_id) today.

Scoring methods are described once, declaratively (`mteapy.methods`:
`Param`/`Method`), and both front ends render from it: the CLI builds its
argparse options from the spec (a test pins the CLI to it) and the GUI builds
its parameter form from `GET /api/methods`. The GUI "Run" scores every task
against *every* sample in a worker thread (a matrix, not one sample at a
time), with progress/cancel; results are in memory only so far. TIDE is not in
the registry (it needs a differential-expression table and permutations).
Decided but not built yet: saved results embed the model-relevant expression
subset (~2.8k genes) so a result is viewable and reproducible on its own.

## Schema v2: self-contained task + route database, integer keys

The legacy database stored routes keyed by `(source, task_id, model_id)`
strings and kept task definitions only in external task files, so a
database could not be read (or trusted) without the exact files it was
built from, and "source" vs "task_list" was two overlapping notions. v2
(`mteapy.taskdb`): one DB per model; model reactions/metabolites stored
once and referenced by integer key (route rows are two integers + a flux,
and the importer verifies every task against the model); a named
`task_lists` row per list; every task's full definition in the DB;
`route_reactions.flux NOT NULL`. 102 MB -> ~30 MB with identical content
(all 17,616 migrated routes' reaction sets and fluxes verified equal to the
legacy ones; TAS scores through the CLI match the old code exactly).
The user's reasoning: users should get one file that carries tasks, routes
and fluxes, with the model verified by hash rather than embedded (46 MB of
XML buys nothing). Task lists are **immutable** once imported -- changing a
task means a new list name -- which makes the old resume-time
definition-hash guard unnecessary (routes can't outlive their definition).
`taskdb.connect` refuses legacy databases; `cmds/migrate_db.py` converts
them. Not yet done: the `HumanGEM/` model-folder layout + manifest/registry
(config instead of hardcoded names), and the remote model/DB fetch.

The legacy `tasks.definition_hash` was NULL for most tasks (registered
before the column existed), so migration derives the expected hash from the
task file the legacy DB recorded (sha256 verified) rather than treating
NULL as "changed" -- the first attempt did the latter and dropped ~16k
routes by mistake.

**Future idea (user):** the enumerator should always save the flux, not
just the reaction support. `mteapy.enumeration` already returns
`RouteResult.fluxes` and the v2 recorder requires them; what still emits
supports only is the older MN5 greasy pipeline (`csmemo`'s
`run_one_task.py`), whose output `cmds/import_routes.py` has to re-solve.
Porting that job to mteapy's enumerator removes the re-solve.

## 24 "Full" degradation tasks had all-zero lower bounds: corrected

Task 31/32 (fructose/galactose) and 77-97, 99 (amino-acid degradation)
in `metabolicTasks_Full.txt` had every IN/OUT lower bound 0, so "do
nothing" was a feasible zero-reaction solution and no expression data
could ever make them score above zero (the "always-zero tasks" case in
`debugging.md`). The CellFie versions of the same tasks force the
substrate (`1 1`). Testing on the full model: forcing the substrate alone
is infeasible (O2 UB was 1, one alanine needs several), substrate `1 1` +
O2 UB 1000 works for 23/25, and also swapping the odd `urate[e]` output
for `urea[e]` makes 24/25 feasible (arginine needs it; I judged the urate
a template typo from the CellFie versions -- an inference, not
documented). Applied in the Human-GEM checkout's task file (uncommitted
edit at the time of writing; tasks' COMMENTS record the change). **Task
246 (PAP degradation)** stays infeasible when forced -- its outputs don't
balance for PAPS -- so it is flagged `INVALID:` in its COMMENTS (becomes
`valid = 0` on import) pending curation. Their routes were re-enumerated
on MN5 (greasy, `gp_debug`, CPLEX) and 23 imported (task 31 timed out with no
routes; 12 of the 23 have duplicate supports in the raw output, so fewer
distinct routes than the cap -- see `current-state.md`); the legacy empty
routes for those 24 tasks were not migrated.

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

*(Superseded by schema v2 above -- `task_lists` replaces `task_sources`; kept for the history.)*

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

## Repeated supports mark a run truncated; saved results are bound to the model and routes

Two decisions from the review of the 24-task import. (1) The 12 tasks whose MN5
output repeated supports keep only their distinct routes and are **not**
recomputed for now: a repeating enumerator most likely means solver numerical
noise left a cut unapplied, which is a property of the solve, not of the task.
Instead the *store* path is guarded: `taskdb.record_enumeration_result` warns
and records `truncated=True, hit_cap=False` for any batch with repeated
supports, so the enumerator's own `degenerate_duplicate` stop and imported
output are treated alike and a resume revisits those tasks. (2) Saved GUI
results (`mteapy-results` v1 JSON) embed the scoring signal restricted to the
genes the task list's routes reach, and are accepted on load only against the
same model sha256 and the same `task_list_routes_fingerprint`. CellFie results
embed the gene activity levels, not the expression, because CellFie thresholds
depend on the whole dataset; they load as scores (and networks) only.
