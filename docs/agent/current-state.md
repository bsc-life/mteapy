# Current state

Last updated: 2026-10-06, end of the session that added TAS and moved
this repo to its own independent dev venv.

## Branch status -- read this first

Active development is on **`feature/context-aware-scoring`**, pushed to
`origin/feature/context-aware-scoring`, clean (nothing uncommitted except
a local-only DB backup file, see below). It is **33 commits ahead of
`main`** and **not yet merged**. `main` is at commit `5e3c003` ("Add
__init__.py to fix namespace-package leakage") -- it does *not* have TAS,
the `task_list` grouping, the webapp source-mixing fix, the CellFie bug
fixes, or the route-database schema additions. If you're asked to "check
out mteapy" or "use the latest", that almost certainly means the feature
branch, not `main` -- confirm which one the task actually needs.

No decision has been made yet about merging. Don't merge unilaterally;
ask first.

## Test suite

152 tests passing (`pytest` from repo root with the dev venv active).

## What's new this session (on the feature branch)

- **TAS** (`compute_TAS`, `run-mtea analyze TAS`) -- see `decisions.md`.
- **`task_sources.task_list`** grouping + retirement of the duplicate
  `cellfie_consensus` source -- see `decisions.md`.
- **Webapp**: fixed cross-source task mixing (every lookup now keyed by
  `(source, task_id)`), added the no-data "show topology only" mode,
  fixed a stale hardcoded GTEx example-dataset path.
- **18 `full`-list tasks** that were stuck at zero routes now have a
  pFBA-only reference route (`truncated=True`) -- see `decisions.md`.
- **Two CellFie algorithm fixes** (percentile-space, missing-gene
  propagation) -- see `decisions.md`. Extensive new test coverage for both.
- Git LFS set up for `routes_human2.db`, including a full consented
  history rewrite (both `main` and the feature branch affected).
- Namespace-package bug fixed (`__init__.py` added) on both branches.

## Known gaps / open questions

- **Model-version mismatch between classic and context-aware paths** is
  known and accepted, not a bug -- see `architecture.md`. Don't "fix" it
  without a deliberate decision to do so.
- **CellFie/task import from the original GPL repo is unresolved** -- see
  `decisions.md`'s last entry. Don't restart that work without
  revisiting the licensing question with the user first.
- **No webapp feature yet for comparing two datasets/samples side by
  side** (e.g. two different conditions' expression for the same
  tissue). Came up while using the webapp to explore an external
  analysis's results; not built, no decision made on whether to build it.
- A DB backup file, `src/mteapy/data/routes_human2.db.bak-pre-timeout-rerun`,
  sits locally, untracked, from before the 18-task pFBA fill. Safe to
  delete once that fill is trusted; not currently gitignored or removed.
- This repo now has its own independent dev venv at `venv/` (Python
  3.10, built from `requirements.txt` + `pip install -e . --no-deps` +
  `gurobipy`/`python-libsbml`/`jupyterlab`/`ipykernel`; kernel name
  `mteapy`). It replaces an earlier arrangement where this repo shared
  one venv with another project -- that older, shared venv silently
  broke once (its `python` symlink started resolving to a system Python
  that had been upgraded out from under it after the venv was created).
  See `debugging.md` for the exact symptom (worth knowing even with an
  independent venv now, since the same failure mode can recur) and
  `workflows.md` for the from-scratch setup steps.
