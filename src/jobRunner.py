#!/usr/bin/env python3
# src/jobRunner.py

"""
Runs one background job the web interface's admin Generate page started
(`html/web/jobs.py` spawns `python3 src/jobRunner.py <job dir>` detached
from the web server, then returns right away).

A job is a list of steps, each one command line (`resetDb.py`, then
`generate.py plan`, then `generate.py galaxy ...`). This runs them in
order, appends every step's output to `<job dir>/output.log`, and keeps
`<job dir>/state.json` current, so the page can show which step is
running and how the job ended:

    {"status": "running" | "succeeded" | "failed" | "cancelled",
     "pid": 1234, "step": 2, "started_at": ..., "finished_at": ...,
     "exit_code": 0, "error": "..."}

It stops at the first step that fails. SIGTERM (the page's Cancel
button) stops the running step and marks the job cancelled. When it is
done it removes the jobs directory's `active` lock, if the lock is still
this job's, so the next job can start.

Standard library only: this starts before anything else is imported, so a
broken install still leaves a readable `state.json` behind.
"""

import json
import os
import signal
import subprocess
import sys
import tempfile
import time

LOCK_NAME = "active"
"""str: The jobs directory's lock file, holding the running job's id."""


def _write_json(path, body):
    directory = os.path.dirname(path)
    fd, tmp = tempfile.mkstemp(prefix=".state-", dir=directory)
    with os.fdopen(fd, "w", encoding="utf-8") as f:
        json.dump(body, f)
    os.replace(tmp, path)


def _release_lock(jobs_dir, job_id):
    lock = os.path.join(jobs_dir, LOCK_NAME)
    try:
        with open(lock, "r", encoding="utf-8") as f:
            holder = f.read().strip()
        if holder == job_id:
            os.remove(lock)
    except OSError:
        pass


def run(job_dir):
    """
    Runs the job in `job_dir` to completion.

    Returns:
        int: 0 when every step succeeded, else 1.
    """
    job_dir = os.path.abspath(job_dir)
    jobs_dir = os.path.dirname(job_dir)
    with open(os.path.join(job_dir, "job.json"), "r", encoding="utf-8") as f:
        job = json.load(f)

    state_path = os.path.join(job_dir, "state.json")
    progress_path = os.path.join(job_dir, "progress.json")
    state = {"status": "running", "pid": os.getpid(), "step": 0, "started_at": time.time(),
             "finished_at": None, "exit_code": None, "error": None}
    _write_json(state_path, state)

    current = {"proc": None, "cancelled": False}

    def _on_term(signum, frame):
        current["cancelled"] = True
        proc = current["proc"]
        if proc is not None and proc.poll() is None:
            try:
                os.killpg(proc.pid, signal.SIGTERM)
            except OSError:
                pass

    signal.signal(signal.SIGTERM, _on_term)
    signal.signal(signal.SIGINT, _on_term)

    env = dict(os.environ)
    env.update(job.get("env") or {})
    env["PYTHONUNBUFFERED"] = "1"
    env["PLANETGEN_PROGRESS_FILE"] = progress_path

    steps = job["steps"]
    try:
        with open(os.path.join(job_dir, "output.log"), "ab", buffering=0) as log:
            for index, step in enumerate(steps, start=1):
                if current["cancelled"]:
                    break
                state["step"] = index
                _write_json(state_path, state)
                try:
                    os.remove(progress_path)
                except OSError:
                    pass
                header = f"\n=== Step {index} of {len(steps)}: {step['label']} ===\n"
                log.write(header.encode("utf-8"))
                started = time.time()
                # Its own process group, so Cancel stops the step's worker
                # processes (`plan`'s multiprocessing pool) too.
                current["proc"] = subprocess.Popen(
                    step["argv"], stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                    cwd=job.get("cwd") or None, env=env, start_new_session=True,
                )
                if current["cancelled"]:
                    _on_term(signal.SIGTERM, None)
                code = current["proc"].wait()
                current["proc"] = None
                log.write(f"=== Step {index} exited with status {code} after "
                          f"{time.time() - started:.0f} s ===\n".encode("utf-8"))
                state["exit_code"] = code
                if current["cancelled"]:
                    break
                if code != 0:
                    state["status"] = "failed"
                    state["error"] = f"{step['label']} exited with status {code}."
                    break
        if current["cancelled"]:
            state["status"] = "cancelled"
            state["error"] = "Cancelled by an admin."
        elif state["status"] == "running":
            state["status"] = "succeeded"
    except Exception as exc:  # noqa: BLE001 -- recorded for the page, never lost
        state["status"] = "failed"
        state["error"] = f"The job runner failed: {exc}"
    finally:
        state["finished_at"] = time.time()
        _write_json(state_path, state)
        _release_lock(jobs_dir, job["id"])
    return 0 if state["status"] == "succeeded" else 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: jobRunner.py <job dir>")
    sys.exit(run(sys.argv[1]))
