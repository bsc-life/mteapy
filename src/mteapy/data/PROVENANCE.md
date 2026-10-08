# Provenance of bundled data

mteapy ships two independent generations of Human-GEM data, for its two
independent scoring methods. They are deliberately different model
versions -- see `run-mtea tasks enumerate-routes`'s docstring for why
the newer generation doesn't try to reproduce the older one's numbers.

## TIDE / CellFie generation (older)

`HumanGEM.xml.gz`, `HumanGEM_genes.tsv`, `HumanGEM_essential_genes_matrix.tsv`,
`task_structure_matrix.tsv`, `task_metadata.tsv`.

- sha256(`HumanGEM.xml.gz`) = `9120ded91ff2830593d5a2f160aa313a321fa8378252a6ee6a99ab6fd86ee08d`
- Exact originating Human-GEM commit/tag is not recorded anywhere in this
  repo's history; these files predate the provenance-tracking convention
  below. Known only from mteapy's own docstrings: ~13085 reactions, the
  model version the original TIDE/CellFie reproduction and the AGS-TIDE
  paper (github.com/bsc-life/ags-paper) were validated against.

## Context-aware scoring generation (newer)

`models/HumanGEM/HumanGEM_v201.xml`, `models/HumanGEM/routes.db` (the whole folder is described by its `manifest.json`).

- sha256(`models/HumanGEM/HumanGEM_v201.xml`) = `6ce49b620391f0ad76be24fdbdfc884fa376ff534dd3813c95abbe9d8c66fa5e`
  (12877 reactions, 2848 genes -- this exact hash is also recorded in
  `models/HumanGEM/routes.db`'s own `models` table, so the database's routes are
  independently verifiable as having been computed against this file and
  no other).
- Filename intent: corresponds to Human-GEM's `v2.0.1` release. This
  **cannot be verified against a clean tag**, though: Human-GEM ships no
  SBML export in its own repo (this file was generated via a COBRApy
  export at some earlier point, not committed here), and the sibling
  `Human-GEM/` checkout this project otherwise uses sits 72 commits past
  the `v2.0.0` tag -- past `v2.0.1` too, but including this project's own
  task-list curation fixes on top (e.g. "Fix task 187 duplicate-OUT
  structural parse bug", "Fix task 90 and 148 structural parse bugs in
  CellFie task list"). The exact commit this XML was originally exported
  from is not retrievable after the fact.
- `models/HumanGEM/routes.db` is a schema-v2 task + route database
  (`mteapy.taskdb`): its `task_lists` table records, per list (e.g.
  `HumanGEM-Full`, `CellFie`), the exact task-list file path, its sha256,
  and (when available) the origin repo URL and git commit it was read
  from -- populated automatically by `run-mtea tasks import`. The task
  definitions themselves are stored in the database.

## Going forward

Any future regeneration of the bundled model/routes DB should record, at
minimum, the exporting script's own `git rev-parse HEAD` (of whichever
repo the source model came from) at export time -- exactly what
`run-mtea tasks import`'s `_git_provenance()` helper does
automatically for the task-list file already. Retrofit the model export
step the same way if it's ever automated.
