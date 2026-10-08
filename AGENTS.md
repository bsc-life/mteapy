# mteapy -- agent instructions

Python library + CLI for Metabolic Task Enrichment Analysis: CellFie, TIDE,
TIDE-essential, and TAS (Task Activity Score), plus a context-aware
alternate-route/complex-level scoring layer and a local FastAPI visualizer
(`webapp/`). This repo is the **tool** -- a general-purpose scoring
engine, independent of any specific analysis project that uses it. If
you're here to run or interpret a scientific analysis built on top of
mteapy rather than to change how scoring itself works, that work belongs
in whichever project depends on this package (installed via pip, e.g.
`mteapy @ git+ssh://git@github.com/bsc-life/mteapy.git`) -- look there
for that project's own agent docs.

Full context lives in `docs/agent/`, read on demand rather than inline here:

- `docs/agent/architecture.md` -- package layout, the routes database
  schema, the classic/context-aware scoring axis, the webapp.
- `docs/agent/decisions.md` -- why things are built the way they are,
  with the reasoning, not just the what.
- `docs/agent/current-state.md` -- branch status, what's merged vs. not,
  known gaps. **Read this first in any new session** -- the single
  biggest thing to get right here is which branch you're on.
- `docs/agent/debugging.md` -- specific failure modes already hit once
  (namespace-package shadowing, a venv silently running the wrong
  Python, a stale global Jupyter config) and how to recognize them fast.
- `docs/agent/workflows.md` -- concrete steps for the tasks that come up
  repeatedly (dev environment setup, running tests, adding a new
  `analyze` method, route enumeration, the webapp).

## Standing rules (always apply, not just "ask once")

- **Never publish to PyPI without the user's explicit authorization for
  that specific publish.** A prior authorization does not carry forward.
- **Solver choice is environment-dependent, not a free choice:** Gurobi
  for anything running locally; CPLEX on BSC's MN5 cluster. Ask before
  defaulting to GLPK anywhere real work depends on the result.
- **Never credit Claude/an AI as a commit co-author or in commit
  trailers**, on this repo or any repo the user owns.
- **Only commit when explicitly asked.** Staging/diffing freely is fine;
  creating the commit is not, even when the change is obviously correct.
- `src/mteapy/data/models/*/routes.db` (the per-model task/route databases) are **Git LFS**-tracked. A fresh
  clone needs `git lfs install` before the real file resolves (otherwise
  you'll silently get a pointer stub, not the database).
- Don't re-import the original CellFie repo's `.mat` model/task files
  without re-raising the GPL question explicitly -- format conversion
  (`.mat` -> SBML/`.txt`) does not escape GPL copyleft on the content.
  This was deliberately deferred, not resolved; treat it as still open.
