"""Background analysis runs for the webapp.

A run scores every task of a task list against every sample of a dataset
with one method and one parameter set. It executes in a worker thread so
the (single-process) server keeps answering requests while it computes, and
reports progress that the page polls. Runs are serialized (one worker): the
work is CPU-bound Python, so overlapping runs would only slow each other
down. Finished runs are kept in memory (the most recent `max_runs`) for the
lifetime of the server.
"""

from __future__ import annotations

import threading
import time
import traceback
import uuid
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field

from mteapy.context_scoring import RunCancelled


@dataclass
class Run:
    run_id: str
    config: dict                     # model, task_list, dataset_id, method, params
    status: str = "queued"           # queued | running | done | error | cancelled
    done: int = 0
    total: int = 0
    error: str | None = None
    genes_matched: int | None = None
    created: float = field(default_factory=time.time)
    started: float | None = None
    finished: float | None = None
    result: dict | None = None       # {"samples", "tasks", "signal", "aggregation", "or_func"} once done
    cancel: threading.Event = field(default_factory=threading.Event)

    def info(self) -> dict:
        return {
            "run_id": self.run_id, "status": self.status, "done": self.done, "total": self.total,
            "error": self.error, "genes_matched": self.genes_matched, "config": self.config,
            "seconds": (None if self.started is None else round((self.finished or time.time()) - self.started, 2)),
        }


class RunManager:
    def __init__(self, max_runs: int = 20):
        self._runs: OrderedDict[str, Run] = OrderedDict()
        self._lock = threading.Lock()
        self._executor = ThreadPoolExecutor(max_workers=1, thread_name_prefix="mteapy-run")
        self._max_runs = max_runs

    def submit(self, config: dict, job, genes_matched: int | None = None) -> Run:
        """Queue `job(run)` -- it must return the run's result dict, call
        `check_cancelled(run)` / update `run.done`/`run.total` as it goes."""
        run = Run(uuid.uuid4().hex[:12], config, genes_matched=genes_matched)
        with self._lock:
            self._runs[run.run_id] = run
            while len(self._runs) > self._max_runs:
                oldest = next(iter(self._runs))
                if self._runs[oldest].status in ("queued", "running"):
                    break
                self._runs.popitem(last=False)
        self._executor.submit(self._execute, run, job)
        return run

    def add_finished(self, config: dict, result: dict, genes_matched: int | None = None) -> Run:
        """Register an already-computed run (results loaded from a file)."""
        run = Run(uuid.uuid4().hex[:12], config, genes_matched=genes_matched)
        run.status, run.result = "done", result
        run.started = run.finished = time.time()
        with self._lock:
            self._runs[run.run_id] = run
            while len(self._runs) > self._max_runs:
                oldest = next(iter(self._runs))
                if self._runs[oldest].status in ("queued", "running"):
                    break
                self._runs.popitem(last=False)
        return run

    def _execute(self, run: Run, job) -> None:
        if run.cancel.is_set():
            run.status, run.finished = "cancelled", time.time()
            return
        run.status, run.started = "running", time.time()
        try:
            run.result = job(run)
            run.status = "done"
        except RunCancelled:
            run.status = "cancelled"
        except Exception as exc:  # noqa: BLE001 -- surfaced to the page, traceback to the server log
            traceback.print_exc()
            run.status, run.error = "error", f"{type(exc).__name__}: {exc}"
        finally:
            run.finished = time.time()

    def get(self, run_id: str) -> Run | None:
        with self._lock:
            return self._runs.get(run_id)

    def cancel(self, run_id: str) -> Run | None:
        run = self.get(run_id)
        if run is not None and run.status in ("queued", "running"):
            run.cancel.set()
        return run


def check_cancelled(run: Run) -> None:
    if run.cancel.is_set():
        raise RunCancelled()
