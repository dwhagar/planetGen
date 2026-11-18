"""OPS.244, slice B: what the tests of the Generate page's jobs share. A job is a run of the job ledger (`planetgen/queue/pagerun.py`), so a test that needs a job in a state writes
the rows the runner and the launcher would, and a test that needs a job to run runs its steps in a thread of its own (`runs_in_threads`), as the supervisor's forked job process would."""

import contextlib
import json
import os
import threading
import time

from planetgen.queue import controls, ledger, pagerun, supervisor  # noqa: F401 -- `pagerun` and `supervisor` are used through this module by the tests
from planetgen.web import job_runner, jobs

DEAD_PID = 2_000_000_000
"""int: A process id no process has."""


@contextlib.contextmanager
def ledger_rows():
    """A connection to the control database of the test session."""
    control = pagerun.connect()
    try:
        yield control
    finally:
        control.close()


def forget_all_jobs():
    """Deletes every Generate page job of the ledger (the session's control database is shared by the tests of a process)."""
    with ledger_rows() as control:
        control.execute("DELETE FROM engine_runs WHERE kind = ?", (pagerun.KIND,))
        control.commit()


def backdate(job_id, seconds=3600):
    """Makes the job `seconds` older, as if it had been started long ago."""
    with ledger_rows() as control:
        found = control.execute("SELECT j.id, j.arguments FROM engine_jobs j WHERE j.stage = ? AND j.command = ? AND j.item_key = ?", (pagerun.STAGE, pagerun.JOB, job_id)).fetchone()
        named = json.loads(found["arguments"])
        named["created_at"] -= seconds
        control.execute("UPDATE engine_jobs SET arguments = ? WHERE id = ?", (json.dumps(named), found["id"]))
        control.commit()


def _steps(control, run):
    return [row["id"] for row in control.execute("SELECT id FROM engine_jobs WHERE run_id = ? AND command = ? ORDER BY item_key", (run, pagerun.STEP)).fetchall()]


def end(job_id, status="succeeded", step=1, error="the step failed"):
    """Writes the end of the job as its runner and the launcher would (`succeeded`, `failed` at step number `step`, `canceled`, or `interrupted`: step `step` running in a process that is gone), and frees the job's lock."""
    with ledger_rows() as control:
        run = pagerun.run_of(control, job_id)
        rows = _steps(control, run)
        if status == "succeeded":
            for row in rows:
                ledger.begin_attempt(control, row, pid=1)
                pagerun.step_done(control, row, run)
        elif status == "canceled":
            controls.end_rows(control, run)
            control.commit()
        else:
            for row in rows[:step - 1]:
                ledger.begin_attempt(control, row, pid=1)
                pagerun.step_done(control, row, run)
            ledger.begin_attempt(control, rows[step - 1], pid=DEAD_PID if status == "interrupted" else 1)
            if status == "failed":
                pagerun.step_failed(control, rows[step - 1], run, error)
            elif status == "interrupted":
                pass
            else:
                raise ValueError(status)
    jobs.release_lock(job_id)


def running(job_id, step=1, pid=None):
    """Makes step number `step` of the job the one that runs, in the process `pid` (this test's own by default)."""
    import os

    with ledger_rows() as control:
        run = pagerun.run_of(control, job_id)
        rows = _steps(control, run)
        for row in rows[:step - 1]:
            ledger.begin_attempt(control, row, pid=1)
            pagerun.step_done(control, row, run)
        pagerun.begin(control, run)
        ledger.begin_attempt(control, rows[step - 1], pid=pid or os.getpid())


def run_to_end(job_id, **kwargs):
    """Runs the job's steps in this thread, as the launcher would hand them out; returns 0 when every step succeeded."""
    return job_runner.run_to_end(job_id, **kwargs)


def runs_in_threads(monkeypatch):
    """Makes `jobs.start_job` hand a job to a thread of the test that runs its steps (the supervisor's forked process, in a test that has none); returns the list of threads started."""
    threads = []

    def hand_over(_connection, job_id):
        thread = threading.Thread(target=run_to_end, args=(job_id,), daemon=True, name=f"page-job-{job_id}")
        thread.start()
        threads.append(thread)

    monkeypatch.setattr(jobs, "hand_over", hand_over)
    return threads


def run_steps(root, steps, kind="galaxy", title="A job", database=None):
    """Starts a job of `steps` with its folder under `root`, runs it to its end in this thread and returns `(job id, its folder, 0 when every step succeeded)`."""
    forget_all_jobs()
    job_id = jobs.start_job(kind, title, steps, database=database, root=str(root), spawn=False)
    code = run_to_end(job_id)
    return job_id, os.path.join(str(root), job_id), code


def wait_finished(job_id, root, timeout=60):
    """Waits until the job is over and its lock is free; returns `jobs.get_job`'s dict."""
    deadline = time.time() + timeout
    while time.time() < deadline:
        job = jobs.get_job(job_id, root)
        if job and job["finished"] and jobs._read_lock(jobs._connection()) != job_id:
            return job
        time.sleep(0.05)
    raise AssertionError(f"job {job_id} did not finish: {jobs.get_job(job_id, root)}")
