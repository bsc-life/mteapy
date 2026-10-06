# Debugging notes

Specific failure modes already hit once, with the actual symptom so
they're recognizable fast next time, not just the fix.

## `import mteapy` silently resolves to the wrong thing (namespace package)

**Symptom**: `import mteapy` succeeds, but `mteapy.__file__` is `None`
and/or `mteapy.routes`/other submodules fail with `ModuleNotFoundError`
even though the package is pip-installed correctly. `mteapy.__path__`
shows a `_NamespacePath` pointing somewhere unexpected, e.g. the
workspace root rather than the installed package's `src/mteapy`.

**Cause**: Python's `''` (cwd) `sys.path` entry, combined with running
from a directory that has a *sibling* directory literally named
`mteapy/` with no `__init__.py` of its own (the outer git-repo-root
`mteapy/` directory, one level up from `src/mteapy/`, counts). Python
treats that as a namespace package and it can shadow or merge with the
real installed package depending on exact `sys.path` ordering.

**Fix already applied**: `__init__.py` added to `src/mteapy/` and
`src/mteapy/cmds/` on both `main` and the feature branch. If this
resurfaces, check whether a *new* directory literally named `mteapy`
without an `__init__.py` has been created somewhere on the Python path
(e.g. a careless `mkdir mteapy` while scaffolding something), not just
whether the installed package itself is intact.

## The dev venv silently runs the wrong Python version

**Symptom**: `cobra`/`gurobipy`/`mteapy` imports fail with
`ModuleNotFoundError`, but `pip list` in the "same" venv shows almost
nothing installed, or `pip` itself is missing
(`ModuleNotFoundError: No module named 'pip'` when running
`<venv>/bin/pip`). Everything *looks* activated correctly.

**Cause**: the venv's `bin/python` is a symlink to a *system* Python
(e.g. `/usr/bin/python3`), not a private copy. If the OS's system Python
gets upgraded (e.g. 3.10 -> 3.12) after the venv was created, the
symlink silently starts resolving to the new version, while all the
venv's actual installed packages still live under the *old* version's
`lib/pythonX.Y/site-packages/` -- which the new interpreter never looks
in. No error at venv-creation or activation time; it only surfaces when
something tries to import a package.

**Diagnosis**: `cat <venv>/pyvenv.cfg` and compare its recorded
`version_info` against `<venv>/bin/python --version`. If they disagree,
or if `ls <venv>/lib/` shows *two* different `python3.X` directories
(one with real packages, one nearly empty), this is it.

**Fix**: rebuild the venv from scratch against a Python binary matching
whatever the *existing* installed packages (especially anything
requiring a prebuilt binary egg, like a CPLEX install) were actually
built for -- see `workflows.md`'s "Dev environment setup" for the exact
steps used last time (Python 3.10, matching the CPLEX egg found in the
old venv's site-packages).

## `jupyter nbconvert` crashes with an unrelated `ModuleNotFoundError`

**Symptom**: `jupyter nbconvert --execute ...` fails immediately with
something like `ModuleNotFoundError: No module named
'jupyter_contrib_nbextensions'`, even for a notebook that has nothing to
do with that package.

**Cause**: a **global**, user-level `~/.jupyter/jupyter_nbconvert_config.json`
(left over from a *different*, unrelated project on the same machine)
registers `jupyter_contrib_nbextensions.nbconvert_support.*`
preprocessors unconditionally, and that package isn't installed in this
project's venv.

**Fix**: do **not** edit or delete that global config -- it belongs to
another project. Override just for this invocation:
`JUPYTER_CONFIG_DIR=<empty/throwaway dir> jupyter nbconvert ...`.

## A task scores exactly zero in every tissue/sample -- don't assume it's a scoring bug

Seen twice, two different real causes, neither a scoring bug:

1. **A gene-id mismatch**: the GPR's gene id is absent from the
   expression reference entirely (not just low-expressed) -- check
   whether the HGNC symbol exists in the expression data under a
   *different* id before suspecting the model or the scoring code.
   (Real example: Human-GEM's GPR for `SLC37A4` used a stale Ensembl id;
   the same gene has substantial real expression under a different id.)
2. **A task-definition bounds issue**: if the task's single enumerated
   route has **zero real reactions** in it at all (check
   `load_task_routes`'s result size directly, not just the score), the
   task's own `lower_bound`s on its IN/OUT metabolites may all be 0.0 --
   meaning "do nothing" is a trivially feasible, zero-cost solution, and
   no expression data under any method could ever make it score
   non-zero. This is a property of the task list, not of the tissue/sex
   being scored.

Always check which of these it is (gene-id lookup vs. route reaction
count) before concluding the scoring pipeline itself is broken.

## CPLEX install needs `setup.py`, not pip

The CPLEX Python bindings matching a local IBM CPLEX Studio install
(`/opt/ibm/ILOG/CPLEX_Studio221/cplex/python/<pyver>/x86-64_linux/`) are
installed by running `python setup.py install` from inside that
directory -- there's no pip package for this. Deprecation warnings from
`setup.py`/`easy_install` during this are expected and harmless.
`gurobipy`, by contrast, installs cleanly via plain `pip install
gurobipy` and auto-detects a license at `~/gurobi.lic` with no extra steps.
