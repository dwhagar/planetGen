# planetgen/web/job_runner.py

"""
Runs one background job the web interface's admin Generate page started,
as an RQ job (PERF.24 step 3): `planetgen.web.jobs.start_job` queues
`run(<job dir>)` on a queue of the job's own and starts a burst worker
(`planetgen.cli.worker`) for it, detached from the web server, then
returns right away.

A job is a list of steps, each one command line (`planetgen.cli.reset`, then
`planetgen plan`, then `planetgen galaxy ...`). This runs them in
order, appends every step's output to `<job dir>/output.log`, and keeps
`<job dir>/state.json` current, so the page can show which step is
running and how the job ended:

    {"status": "running" | "succeeded" | "failed" | "cancelled",
     "pid": 1234, "step": 2, "started_at": ..., "finished_at": ...,
     "exit_code": 0, "error": "..."}

`pid` is the job's worker (`runner.pid`), whose command line carries the
job's id, so the page can tell a live job from one whose worker died.

It stops at the first step that fails. The page's Cancel button writes a
`cancel` file into the job directory; this checks for it while a step
runs, stops the step and every process it started, and marks the job
cancelled. (A file rather than RQ's stop command: that kills the worker's
work horse and would leave the step and its processes running.) SIGTERM
does the same (a server shutting down, on POSIX).
When it is done it removes the jobs directory's `active` lock, if the lock is still
this job's, so the next job can start.

The job is also the root of a job tree (ADM.12, `work.open_node`):
a "web-job" node with one "step" node per step, each step's
`planetgen` run hanging its own nodes under its step
(`work.PARENT_ENV_VAR`), so the admin queue page shows the whole
job with timings. That part is best effort.

Standard library only at the top, so a broken install still leaves a
readable `state.json` behind (`_JobTree` imports the rest only once that
is written).
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

PID_NAME = "runner.pid"
"""str: The job's worker's pid, written by `jobs.start_job`."""

CANCELLED_EXIT_CODE = 130
"""int: A step's exit status when an admin cancelled it from the queue
page (must match `work.CANCELLED_EXIT_CODE`)."""

POLL_SECONDS = 0.25
"""float: How often a running step is checked for a cancel request."""


def _write_json(path, body):
    directory = os.path.dirname(path)
    fd, tmp = tempfile.mkstemp(prefix=".state-", dir=directory)
    with os.fdopen(fd, "w", encoding="utf-8") as f:
        json.dump(body, f)
    os.replace(tmp, path)


def _step_process_options():
    """Popen options that give a step its own process group, so stopping
    it can stop the worker processes it starts (`plan`'s pool) too."""
    return {"start_new_session": True}


def _stop_tree(proc):
    """Stops a step and every process it started."""
    if proc is None or proc.poll() is not None:
        return
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except OSError:
        pass


def _worker_pid(job_dir):
    """The job's worker's pid (`PID_NAME`), else this process's."""
    try:
        with open(os.path.join(job_dir, PID_NAME), "r", encoding="utf-8") as f:
            return int(f.read().strip())
    except (OSError, ValueError):
        return os.getpid()


def _drop_queue(job_id):
    """Forgets the job's own RQ queue once it's done (best effort); its
    burst worker then finds it empty and exits."""
    try:
        from planetgen.queue import redisqueue
        from planetgen.web import jobs

        connection = redisqueue.connect()
        redisqueue.queue(jobs.QUEUE_PREFIX + job_id, connection).delete(delete_jobs=False)
    except Exception:  # noqa: BLE001 -- a leftover empty queue is harmless
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
    (`version_key.run_line`). Each `planetgen` step then writes its own
    as its first line. Never raises: without planetGen's modules or the
    database, the line says what it can.
    """
    run = job.get("title") or job["id"]
    try:
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


STAGE_OFFSET_ENV = "PLANETGEN_STAGE_OFFSET"
STAGE_TOTAL_ENV = "PLANETGEN_STAGE_TOTAL"
"""str: The variables a step's process reads to number its stages across the whole job
(`generation.stages.STAGE_OFFSET_ENV`; the runner does not import planetGen's modules at load)."""


def _staged_command(argv):
    """`(command, args)` of a `galaxy` or `plan` command line, else `None`."""
    from planetgen.cli import generate as generate_cli
    from planetgen.web import jobs

    start = len(jobs.GENERATE_COMMAND) + 1
    if len(argv) <= start or argv[1:start] != list(jobs.GENERATE_COMMAND) or argv[start] not in ("galaxy", "plan"):
        return None
    try:
        _parser, parsers = generate_cli.build_parser()
        return argv[start], parsers[argv[start]].parse_args(argv[start + 1:])
    except (SystemExit, Exception):  # noqa: BLE001 -- the step then uses its recorded time
        return None


def _stage_estimates(conn, argv):
    """The expected seconds of each stage of a `galaxy` or `plan` command line (`generation.stages.stage_estimates`),
    `[]` for any other."""
    from planetgen.generation import stages

    parsed = _staged_command(argv)
    return stages.stage_estimates(conn, *parsed) if parsed else []


def _stage_estimate(conn, argv):
    """The summed stage times of a `galaxy` or `plan` command line (`generation.stages.estimate_seconds`), or `None`."""
    from planetgen.generation import stages

    parsed = _staged_command(argv)
    return stages.estimate_seconds(conn, *parsed) if parsed else None


def _step_stage_counts(steps):
    """How many stages each step holds (a step the page gave no stage list is one stage), for numbering a job's steps
    and stages across the whole job."""
    return [max(len(step.get("stages") or []), 1) for step in steps]


def _step_heading(first, count, total):
    """`"Step 4 of 12"`, or `"Steps 4 to 12 of 12"` for a step with several stages."""
    return f"Step {first} of {total}" if count == 1 else f"Steps {first} to {first + count - 1} of {total}"


class _JobTree:
    """The job's nodes in the control database's job tree, or nothing at
    all when planetGen's modules or the control database aren't there."""

    def __init__(self, job):
        self.queue = None
        self.root = None
        try:
            from planetgen.queue import work
            from planetgen.db import store

            base = _mysql_config(job)
            self.queue = work
            self.root = work.open_node(
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

    def step_estimates(self, steps):
        """
        `(per step, per stage)`: the seconds each step and each stage of each step is expected to take (`None` where
        there is no record), for the overall bar (PERF.55): a `galaxy` or `plan` step from the stored time of each
        of its stages with the settings it runs with (PERF.56), any other step, or one with a stage not yet
        recorded, from what the step took in earlier runs. The per-stage list holds one entry for each stage of the
        job, in order (a step with no stage list counts as one).
        """
        labels = [step["label"] for step in steps]
        counts = _step_stage_counts(steps)
        nothing = ([None] * len(labels), [None] * sum(counts))
        if self.queue is None or not self.root.store.available:
            return nothing
        try:
            conn = self.root.store._connect()
            try:
                recorded = self.queue.recorded_step_seconds(conn)
                per_step, per_stage = [], []
                for step, count in zip(steps, counts):
                    argv = step.get("argv") or []
                    per_step.append(_stage_estimate(conn, argv) or recorded.get(step["label"]))
                    found = _stage_estimates(conn, argv)
                    per_stage += found if len(found) == count else [per_step[-1]] + [None] * (count - 1)
            finally:
                conn.close()
        except Exception:  # noqa: BLE001 -- the bar then estimates from the running step alone
            return nothing
        return per_step, per_stage

    def close(self, node, state):
        if self.queue is None or node is None:
            return
        try:
            self.queue.close_node(node, state)
        except Exception:  # noqa: BLE001
            pass


_TREE_STATES = {"succeeded": "done", "failed": "failed", "cancelled": "cancelled"}
"""dict: `state.json` status -> job tree node state."""


def _log_traceback(job_dir):
    """Appends the traceback of the exception being handled to the job's
    output log, so the page's log shows the runner's own failure in full
    (ADM.25). Never raises."""
    try:
        import traceback

        with open(os.path.join(job_dir, "output.log"), "ab", buffering=0) as log:
            log.write(("\n" + traceback.format_exc()).encode("utf-8"))
    except Exception:  # noqa: BLE001 -- the error is in state.json regardless
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
    state = {"status": "running", "pid": _worker_pid(job_dir), "step": 0, "started_at": time.time(),
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
    state["step_estimates"], state["stage_estimates"] = tree.step_estimates(steps)
    _write_json(state_path, state)
    try:
        with open(os.path.join(job_dir, "output.log"), "ab", buffering=0) as log:
            log.write((_run_line(job) + "\n").encode("utf-8"))
            counts = _step_stage_counts(steps)
            total = sum(counts)
            for index, step in enumerate(steps, start=1):
                if _cancel_requested():
                    break
                state["step"] = index
                state["step_started_at"] = time.time()
                _write_json(state_path, state)
                try:
                    os.remove(progress_path)
                except OSError:
                    pass
                first = sum(counts[:index - 1]) + 1
                heading = _step_heading(first, counts[index - 1], total)
                log.write(f"\n=== {heading}: {step['label']} ===\n".encode("utf-8"))
                env[STAGE_OFFSET_ENV] = str(first - 1)
                env[STAGE_TOTAL_ENV] = str(total)
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
                log.write(f"=== {heading.split(' of ')[0]} exited with status {code} after "
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
        _log_traceback(job_dir)
    finally:
        state["finished_at"] = time.time()
        _write_json(state_path, state)
        # The queue goes before the lock: a job that has released its lock
        # is over, queue included (a test or the next job saw the queue
        # a moment after).
        _drop_queue(job["id"])
        _release_lock(jobs_dir, job["id"])
        tree.close(step_node, "failed")
        tree.close(tree.root, _TREE_STATES.get(state["status"], "failed"))
    return 0 if state["status"] == "succeeded" else 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: python -m planetgen.web.job_runner <job dir>")
    sys.exit(run(sys.argv[1]))
