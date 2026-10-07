#!/usr/bin/env python3
# src/jobRunner.py

"""
Runs one background job the web interface's admin Generate page started
(`html/web/jobs.py` spawns `python3 src/jobRunner.py <job dir>` detached
from the web server, then returns right away).

A job is a list of steps, each one command line (`planetgen.cli.reset`, then
`generate.py plan`, then `generate.py galaxy ...`). This runs them in
order, appends every step's output to `<job dir>/output.log`, and keeps
`<job dir>/state.json` current, so the page can show which step is
running and how the job ended:

    {"status": "running" | "succeeded" | "failed" | "cancelled",
     "pid": 1234, "step": 2, "started_at": ..., "finished_at": ...,
     "exit_code": 0, "error": "..."}

It stops at the first step that fails. The page's Cancel button writes a
`cancel` file into the job directory; this checks for it while a step
runs, stops the step and every process it started, and marks the job
cancelled. SIGTERM does the same (a server shutting down, on POSIX).
When it is done it removes the jobs directory's `active` lock, if the lock is still
this job's, so the next job can start.

The job is also the root of a job tree (ADM.12, `workQueue.open_node`):
a "web-job" node with one "step" node per step, each step's
`generate.py` run hanging its own nodes under its step
(`workQueue.PARENT_ENV_VAR`), so the admin queue page shows the whole
job with timings. That part is best effort.

Standard library only at the top: this starts before anything else is
imported, so a broken install still leaves a readable `state.json`
behind (`_JobTree` imports the rest only once that is written).
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

CANCEL_NAME = "cancel"
"""str: The file in a job's directory that asks it to stop (must match
`jobs.CANCEL_NAME`)."""

CANCELLED_EXIT_CODE = 130
"""int: A step's exit status when an admin cancelled it from the queue
page (must match `workQueue.CANCELLED_EXIT_CODE`)."""

POLL_SECONDS = 0.25
"""float: How often a running step is checked for a cancel request."""

WINDOWS = os.name == "nt"


def _write_json(path, body):
    directory = os.path.dirname(path)
    fd, tmp = tempfile.mkstemp(prefix=".state-", dir=directory)
    with os.fdopen(fd, "w", encoding="utf-8") as f:
        json.dump(body, f)
    # On Windows, replacing a file another process has open (the page
    # reading state.json) fails with PermissionError until it's closed,
    # which is a moment later.
    for attempt in range(40):
        try:
            os.replace(tmp, path)
            return
        except PermissionError:
            if attempt == 39:
                os.remove(tmp)
                raise
            time.sleep(0.05)


def _step_process_options():
    """Popen options that give a step its own process group, so stopping
    it can stop the worker processes it starts (`plan`'s pool) too."""
    if WINDOWS:
        return {"creationflags": subprocess.CREATE_NEW_PROCESS_GROUP | subprocess.CREATE_NO_WINDOW}
    return {"start_new_session": True}


def _stop_tree(proc):
    """Stops a step and every process it started."""
    if proc is None or proc.poll() is not None:
        return
    if WINDOWS:
        subprocess.run(["taskkill", "/T", "/F", "/PID", str(proc.pid)],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       creationflags=subprocess.CREATE_NO_WINDOW, check=False)
        return
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except OSError:
        pass


def _release_lock(jobs_dir, job_id):
    lock = os.path.join(jobs_dir, LOCK_NAME)
    try:
        with open(lock, "r", encoding="utf-8") as f:
            holder = f.read().strip()
        if holder == job_id:
            os.remove(lock)
    except OSError:
        pass


def _mysql_config(job):
    """The job's database (`jobs.mysql_env`'s variables), as planetGen's
    `MySQLConfig`; imports planetGen, so only called once `state.json` is
    written."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from planetgen.db import store

    env = job.get("env") or {}
    return store.MySQLConfig(
        host=env.get("PLANETGEN_MYSQL_HOST"), port=env.get("PLANETGEN_MYSQL_PORT"),
        user=env.get("PLANETGEN_MYSQL_USER"), password=env.get("PLANETGEN_MYSQL_PASSWORD"),
        database=env.get("PLANETGEN_MYSQL_DATABASE"),
    )


def _run_line(job):
    """
    The job log's first line (OPS.10): the galaxy seed when the job
    starts, the PlanetGen version with its key, and the job's title
    (`version_key.run_line`). Each `generate.py` step then writes its own
    as its first line. Never raises: without planetGen's modules or the
    database, the line says what it can.
    """
    run = job.get("title") or job["id"]
    try:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from planetgen.db import store
        from planetgen.galaxy import version_key
    except Exception:  # noqa: BLE001 -- the job runs without the line's details
        return f"Galaxy seed unknown, PlanetGen unknown, run: {run}"
    seed = None
    try:
        conn = store.get_connection(_mysql_config(job), ensure_schema=False)
        try:
            seed = store.get_galaxy_seed(conn)
        finally:
            conn.close()
    except Exception:  # noqa: BLE001 -- no database or no galaxy yet
        pass
    return version_key.run_line(seed, run)


class _JobTree:
    """The job's nodes in the control database's job tree, or nothing at
    all when planetGen's modules or the control database aren't there."""

    def __init__(self, job):
        self.queue = None
        self.root = None
        try:
            from stellarObjects import workQueue
            from planetgen.db import store

            base = _mysql_config(job)
            self.queue = workQueue
            self.root = workQueue.open_node(
                "web-job", job.get("title") or job["id"], store.control_mysql_config(base),
                web_job_id=job["id"], database=job.get("database"),
            )
        except Exception:  # noqa: BLE001 -- the job runs without its tree
            self.queue = self.root = None

    def open_step(self, label):
        if self.queue is None:
            return None
        try:
            return self.queue.open_node("step", label)
        except Exception:  # noqa: BLE001
            return None

    def close(self, node, state):
        if self.queue is None or node is None:
            return
        try:
            self.queue.close_node(node, state)
        except Exception:  # noqa: BLE001
            pass


_TREE_STATES = {"succeeded": "done", "failed": "failed", "cancelled": "cancelled"}
"""dict: `state.json` status -> job tree node state."""


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
    cancel_path = os.path.join(job_dir, CANCEL_NAME)

    def _on_term(signum, frame):
        current["cancelled"] = True
        _stop_tree(current["proc"])

    def _cancel_requested():
        if not current["cancelled"] and os.path.exists(cancel_path):
            current["cancelled"] = True
        return current["cancelled"]

    signal.signal(signal.SIGTERM, _on_term)
    signal.signal(signal.SIGINT, _on_term)

    env = dict(os.environ)
    env.update(job.get("env") or {})
    env["PYTHONUNBUFFERED"] = "1"
    env["PLANETGEN_PROGRESS_FILE"] = progress_path

    steps = job["steps"]
    tree = _JobTree(job)
    step_node = None
    try:
        with open(os.path.join(job_dir, "output.log"), "ab", buffering=0) as log:
            log.write((_run_line(job) + "\n").encode("utf-8"))
            for index, step in enumerate(steps, start=1):
                if _cancel_requested():
                    break
                state["step"] = index
                _write_json(state_path, state)
                try:
                    os.remove(progress_path)
                except OSError:
                    pass
                header = f"\n=== Step {index} of {len(steps)}: {step['label']} ===\n"
                log.write(header.encode("utf-8"))
                step_node = tree.open_step(step["label"])
                if step_node is not None:
                    env[tree.queue.PARENT_ENV_VAR] = step_node.id
                started = time.time()
                current["proc"] = subprocess.Popen(
                    step["argv"], stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                    cwd=job.get("cwd") or None, env=env, **_step_process_options(),
                )
                while True:
                    if _cancel_requested():
                        _stop_tree(current["proc"])
                    try:
                        code = current["proc"].wait(timeout=POLL_SECONDS)
                        break
                    except subprocess.TimeoutExpired:
                        pass
                current["proc"] = None
                tree.close(step_node, "cancelled" if current["cancelled"] or code == CANCELLED_EXIT_CODE
                           else ("done" if code == 0 else "failed"))
                step_node = None
                log.write(f"=== Step {index} exited with status {code} after "
                          f"{time.time() - started:.0f} s ===\n".encode("utf-8"))
                state["exit_code"] = code
                if code == CANCELLED_EXIT_CODE and not current["cancelled"]:
                    # Cancelled from the admin queue page (ADM.10).
                    current["cancelled"] = True
                    current["by_queue"] = True
                if current["cancelled"]:
                    break
                if code != 0:
                    state["status"] = "failed"
                    state["error"] = f"{step['label']} exited with status {code}."
                    break
        if current["cancelled"]:
            state["status"] = "cancelled"
            state["error"] = ("Cancelled from the admin queue page." if current.get("by_queue")
                              else "Cancelled by an admin.")
        elif state["status"] == "running":
            state["status"] = "succeeded"
    except Exception as exc:  # noqa: BLE001 -- recorded for the page, never lost
        state["status"] = "failed"
        state["error"] = f"The job runner failed: {exc}"
    finally:
        state["finished_at"] = time.time()
        _write_json(state_path, state)
        _release_lock(jobs_dir, job["id"])
        tree.close(step_node, "failed")
        tree.close(tree.root, _TREE_STATES.get(state["status"], "failed"))
    return 0 if state["status"] == "succeeded" else 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: jobRunner.py <job dir>")
    sys.exit(run(sys.argv[1]))
