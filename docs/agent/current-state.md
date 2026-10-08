# Current state

Last updated: 2026-10-07, end of the session that moved everything to the
schema-v2 task/route database, added the model-folder/manifest registry, and
started the GUI's analysis panel.

## Branch status -- read this first

Active development is on **`feature/context-aware-scoring`**, pushed to
`origin/feature/context-aware-scoring` as of the previous session's last commit
(`745b0d6`) -- **everything described under "This session" below is
UNCOMMITTED in the working tree** (the user commits; don't commit unasked).
`main` is still at `5e3c003` and has none of the feature branch's work; no
decision has been made about merging (the user floated merging to a
`development` branch first, then branching the refactor off it). Don't merge
unilaterally; ask first.

## Test suite

182 tests passing (`pytest` from repo root with the dev venv active).

## This session (uncommitted)

- **Schema v2** (`mteapy.taskdb`): self-contained per-model DB -- model
  entities by integer key, named task lists with full task definitions,
  routes + fluxes (`flux NOT NULL`). Legacy `routes.py` removed; legacy DBs
  convert with `cmds/migrate_db.py` (verified: all 17,616 routes/fluxes
  identical, TAS scores identical through the CLI). See `architecture.md`.
- **CLI**: `tasks import`, `tasks enumerate-routes --task-list`, `models list`,
  `--model` / `--task-list` on `analyze`; `cmds/import_routes.py` loads MN5
  greasy JSON (supports only) and solves fluxes on route-restricted
  submodels (`network.compute_route_fluxes_submodel`, ~1000x faster).
- **Model folder + manifest registry** (`mteapy.registry`): the bundled model
  is `src/mteapy/data/models/HumanGEM/` (manifest, model, `routes.db`, task
  files, annotations); packaged via `pyproject.toml` package-data (verified by
  building a wheel). `.gitattributes` LFS pattern changed to
  `src/mteapy/data/models/*/*.db`.
- **24 "Full" task bounds corrected** (31, 32, 77-97, 99) in the Human-GEM
  checkout's `metabolicTasks_Full.txt` (uncommitted there too, on branch
  `fix/task-90-148-curation`); task 246 flagged `INVALID:` -- see `decisions.md`.
- **Webapp**: model + task-list selectors (per-model contexts), hideable
  panels, an Analysis bar (method + generated parameter form + Run with
  progress/cancel, scoring all samples in a worker). TAS only; see
  `architecture.md`/`webapp/README.md`.
- `mteapy.methods`: one declarative spec per method shared by CLI and GUI.

## Pending -- do these next

- **Task 31 (fructose degradation)**: MN5 job 47042151 hit the 6000 s cap with
  0 routes (`timeout`, truncated; locally feasible, 46 reactions). The user
  deferred it for later revision -- no routes stored, not to be papered over
  with a pFBA-only route without asking.
- **12 imported tasks have fewer than 100 distinct routes** (77, 78, 81, 83-86,
  88, 90, 91, 93, 95): the MN5 output repeated supports. The DB holds only the
  distinct ones (UNIQUE hash), and their `enumeration_runs` rows are now
  `truncated=1, hit_cap=0`. Decision: do not recompute for now; a repeating
  enumerator is probably solver numerical noise leaving a cut unapplied.
  `record_enumeration_result` now warns and records `truncated` whenever a
  batch contains repeated supports.
- **GUI**: TAS and CellFie in the methods registry; save/load results done
  (see `webapp/README.md`); TIDE later (needs DE input and permutations).
  Classic mapping via route 0 needs an explicit `route_rank` column first (see
  `decisions.md`).
- **Servers**: an older copy of the webapp may still be running on port 8765
  from before this refactor -- restart it to pick up the new code.
- The MN5 pipeline (`csmemo`-based `run_one_task.py`) saves route supports
  only; porting it to `mteapy.enumeration` (which returns fluxes) would
  remove the re-solve in `import_routes`.

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
