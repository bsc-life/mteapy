# Current state

Last updated: 2026-11-11, end of the session that added `mteapy.results`, the CLI `--bundle` option and the
Compute/Results tabs in the webapp.

## Branch status -- read this first

Active development is on **`feature/context-aware-scoring`**; everything is committed and pushed
(`origin/feature/context-aware-scoring` at `8153a3c`, tree clean). `main` is still at `5e3c003` and has none of this
work; no decision has been made about merging (the user floated merging to a `development` branch first). Don't merge
unilaterally; ask first. Don't commit unasked.

## Recently done

- Schema-v2 task/route DB, model-folder registry, method specs shared by CLI and GUI (see `architecture.md`).
- Scoring semantics: no-GPR reactions are left out of aggregates, `missing_genes` (`exclude` default | `zero`) for every
  method (see `decisions.md`). Outputs made before this (e.g. in `metabolic-variability`) must be regenerated.
- `mteapy.results` + `run-mtea analyze TAS|CellFie|TIDE ... --mapping-strategy context-aware --bundle FILE[.gz]`; the
  webapp has Compute and Results tabs and loads those files (`.gz` decompressed in the browser). TIDE-essential is
  legacy (CLI only). 4 CLI-bundle tests added in `tests/test_webapp.py`.
- Size measurement of bundles (TAS, Human-GEM, `HumanGEM-Full`, GTEx v10 single samples): the routes reach only 1,246
  genes; gzip JSON is ~4.7 KB/sample (100 samples = 0.47 MB, extrapolated ~5 MB for 1,000 and ~80 MB for 17,000);
  scoring takes ~0.36 s/sample, so compute time, not file size, is the limit. A float32 binary saves only ~10%.
  Parquet/HDF5 not worth it for now; if tens of thousands of samples ever matter, serve the signal per sample from
  the server instead of changing format.

## Open ideas (nothing decided)

- CLI export of task graphs like the web's SVG (options: Graphviz, Node+dagre, headless browser).
- Vendor d3/dagre into `webapp/static` so the page works offline (today they load from cdnjs/jsdelivr; a blocked CDN
  makes tasks/networks silently not show).
- A bundle option that omits the signal (scores only, ~0.5 KB/sample; the web then could not draw networks).
- Side-by-side comparison of two datasets/samples in the webapp.

## Test suite

233 passing (`venv/bin/python -m pytest`).

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
- **GUI**: TAS, CellFie and TIDE (context-aware; TIDE-essential stays CLI-only)
  in the methods registry; save/load results and the pre-run input checks
  (`mteapy.inputcheck`, shared with the CLI) done (see `webapp/README.md`).
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
