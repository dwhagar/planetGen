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
When it is done it frees the job lock in Redis (`jobs.release_lock`), if the lock is still
this job's, so the next job can start.

How long each step took is kept by its label in the jobs directory (`STEP_TIMES_NAME`), for the next job's estimates.

Standard library only at the top, so a broken install still leaves a
readable `state.json` behind (the estimates import the rest only once that
is written).
"""

import json
import os
import signal
import subprocess
import sys
import tempfile
import time

from planetgen.util import logrecord, memwatch

CANCEL_NAME = "cancel"
"""str: The file in a job's directory that asks it to stop (must match
`jobs.CANCEL_NAME`)."""

PID_NAME = "runner.pid"
"""str: The job's worker's pid, written by `jobs.start_job`."""

STEP_TIMES_NAME = "step-times.json"
"""str: The file in the jobs directory with how long each kind of step took in the jobs that succeeded (PERF.33): `{label:
[how many times, total seconds]}`, read for the next job's estimates (`recorded_step_seconds`)."""

MEMORY_NAME = "memory.json"
"""str: The file in a job's directory with the memory peak of each step that ran (PERF.183): a list of
`{"step", "label", "first_stage", "stages", "seconds", "peaks"}` (`memwatch.PeakTracker.peaks`)."""

CANCEL_KILL_SECONDS = 20
"""float: How long a step that was asked to stop (SIGTERM to its process group) may keep running before the whole group is
killed (SIGKILL, PERF.186). A step that catches SIGTERM and waits for its units, or hangs in a database call, used to
run on for ever after Cancel."""

POLL_SECONDS = 0.25
"""float: How often a running step is checked for a cancel request."""

STOP_WAIT_SECONDS = 10
"""float: How long a step that RQ's time limit ended gets to exit after SIGTERM before it is killed."""


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


def _kill_tree(proc):
    """Kills a step and every process in its group, for a step that ignored the request to stop."""
    if proc is None or proc.poll() is not None:
        return
    try:
        os.killpg(proc.pid, signal.SIGKILL)
    except OSError:
        pass


def _reap(proc):
    """Waits a few seconds for a step `_stop_tree` stopped to be gone, so it is not left a zombie of this process."""
    if proc is None:
        return
    try:
        proc.wait(timeout=STOP_WAIT_SECONDS)
    except subprocess.TimeoutExpired:
        proc.kill()


def _record_memory(job_dir, log, watch, index, label, first, count):
    """Writes the memory line after a step into the log and the step's entry into `MEMORY_NAME` (PERF.183). Never fails the job."""
    try:
        log.write((watch.line() + "\n").encode("utf-8"))
        path = os.path.join(job_dir, MEMORY_NAME)
        try:
            with open(path, "r", encoding="utf-8") as f:
                entries = json.load(f)
        except (OSError, ValueError):
            entries = []
        entries.append({"step": index, "label": label, "first_stage": first, "stages": count,
                        "seconds": round(watch.peaks["seconds"], 1), "peaks": watch.peaks})
        _write_json(path, entries)
    except Exception:  # noqa: BLE001 -- a job is never lost to its memory note
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


def _release_lock(job_id):
    """Frees the job lock if it names this job (best effort; the next job clears a stale one)."""
    try:
        from planetgen.web import jobs

        jobs.release_lock(job_id)
    except Exception:  # noqa: BLE001 -- a lock left behind is cleared by the next job
        pass


def _start_waiting_job(jobs_dir, job_dir):
    """Starts the job waiting behind this one (ADM.123, best effort: an admin page asking for the running job starts
    it too)."""
    try:
        from planetgen.web import jobs

        jobs.promote_queued(jobs_dir)
    except Exception:  # noqa: BLE001 -- the job stays queued
        _log_traceback(job_dir)


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


JOB_ENV = "PLANETGEN_JOB_ID"
"""str: The variable naming the job a step belongs to, so a run it starts records the job (PERF.99;
`planetgen.queue.wrapper.JOB_ENV`)."""

RECORDS_FILE = "records.log"
"""str: The file in a job's folder that holds the job's log as records, one per line (UX.122): the runner's stage marks
and every line its steps log (`log.set_record_file`), the same lines `output.log` holds as the page shows them."""

RECORD_FILE_ENV = "PLANETGEN_RECORD_FILE"
"""str: Set to the job's records file for every step, so the step's command writes its log lines there too."""

STAGE_OFFSET_ENV = "PLANETGEN_STAGE_OFFSET"
STAGE_TOTAL_ENV = "PLANETGEN_STAGE_TOTAL"
"""str: The variables a step's process reads to number its stages across the whole job
(`generation.stages.STAGE_OFFSET_ENV`; the runner does not import planetGen's modules at load)."""


def _staged_command(argv):
    """`(command, args)` of a command line of a command that has stages (`commands.STAGED_COMMANDS`), else `None`."""
    from planetgen.cli import generate as generate_cli
    from planetgen.queue import commands
    from planetgen.web import jobs

    start = len(jobs.GENERATE_COMMAND) + 1
    if len(argv) <= start or argv[1:start] != list(jobs.GENERATE_COMMAND) \
            or argv[start] not in commands.STAGED_COMMANDS:
        return None
    try:
        _parser, parsers = generate_cli.build_parser()
        return argv[start], parsers[argv[start]].parse_args(argv[start + 1:])
    except (SystemExit, Exception):  # noqa: BLE001 -- the step then uses its recorded time
        return None


def _job_entries(argv):
    """`[(command, args)]` of a `planetgen run-job` command line (`generation.job.parse_entries`), else `None`."""
    from planetgen.cli import generate as generate_cli
    from planetgen.generation import job
    from planetgen.web import jobs

    start = len(jobs.GENERATE_COMMAND) + 1
    if len(argv) <= start or argv[1:start] != list(jobs.GENERATE_COMMAND) or argv[start] != "run-job":
        return None
    try:
        _parser, parsers = generate_cli.build_parser()
        shared = parsers["run-job"].parse_args(argv[start + 1:])
        return job.parse_entries(shared.spec, shared, parsers, generate_cli.validate_command)
    except (SystemExit, Exception):  # noqa: BLE001 -- the step then uses its recorded time
        return None


def _stage_estimates(conn, argv):
    """The expected seconds of each stage of a `galaxy`, `plan` or `run-job` command line (`generation.stages.stage_estimates`,
    `generation.job.stage_estimates`), `[]` for any other."""
    from planetgen.generation import job, stages

    entries = _job_entries(argv)
    if entries is not None:
        return job.stage_estimates(conn, entries)
    parsed = _staged_command(argv)
    return stages.stage_estimates(conn, *parsed) if parsed else []


def _stage_estimate(conn, argv):
    """The summed stage times of a `galaxy`, `plan` or `run-job` command line (`generation.stages.estimate_seconds`), or `None`
    unless every stage has one."""
    from planetgen.generation import stages

    if _job_entries(argv) is not None:
        found = _stage_estimates(conn, argv)
        return None if not found or any(seconds is None for seconds in found) else sum(found)
    parsed = _staged_command(argv)
    return stages.estimate_seconds(conn, *parsed) if parsed else None


def _step_stage_counts(steps):
    """How many stages each step holds (a step the page gave no stage list is one stage), for numbering a job's steps
    and stages across the whole job."""
    return [max(len(step.get("stages") or []), 1) for step in steps]


def recorded_step_seconds(jobs_dir):
    """
    How long each kind of step took in the jobs that ran before (PERF.33): the mean seconds of every step that succeeded, by its
    label, from `STEP_TIMES_NAME` in `jobs_dir`; empty where there is no record.

    Returns:
        dict: `{label: mean seconds}`.
    """
    try:
        with open(os.path.join(jobs_dir, STEP_TIMES_NAME), "r", encoding="utf-8") as handle:
            kept = json.load(handle)
        return {label: total / count for label, (count, total) in kept.items() if count}
    except (OSError, ValueError, TypeError):
        return {}


def _record_step_seconds(jobs_dir, label, seconds):
    """Adds one succeeded step's `seconds` to the mean the jobs directory keeps for `label`. Best effort: a lost sample
    only makes the next estimate a little less sure."""
    path = os.path.join(jobs_dir, STEP_TIMES_NAME)
    try:
        try:
            with open(path, "r", encoding="utf-8") as handle:
                kept = json.load(handle)
        except (OSError, ValueError):
            kept = {}
        count, total = kept.get(label, [0, 0.0])
        kept[label] = [count + 1, total + seconds]
        _write_json(path, kept)
    except (OSError, ValueError, TypeError):
        pass


def step_estimates(conn, steps, recorded):
    """
    `(per step, per stage)`: the seconds each step and each stage of each step is expected to take (`None` where
    there is no record), for the overall bar (PERF.55): a `galaxy` or `plan` step from the stored time of each
    of its stages with the settings it runs with (PERF.56), any other step, or one with a stage not yet
    recorded, from what the step took in earlier runs. The per-stage list holds one entry for each stage of the
    job, in order (a step with no stage list counts as one).

    Args:
        conn: A control database connection (the stored stage times).
        steps (list[dict]): The job's steps.
        recorded (dict): What each step took in earlier jobs (`recorded_step_seconds`).
    """
    counts = _step_stage_counts(steps)
    per_step, per_stage = [], []
    for step, count in zip(steps, counts):
        argv = step.get("argv") or []
        per_step.append(_stage_estimate(conn, argv) or recorded.get(step["label"]))
        found = _stage_estimates(conn, argv)
        per_stage += found if len(found) == count else [per_step[-1]] + [None] * (count - 1)
    return per_step, per_stage


def expected_seconds(job):
    """
    How long `job` (as `job.json` holds it) is expected to take (PERF.144): the sum of its steps' expected seconds
    (`step_estimates`: the recorded stage times, else what the step took in earlier runs); `None` unless every step
    has one, or without a control database. The web server calls this to give the job's RQ entry a time limit.
    """
    try:
        from planetgen.db import store

        from planetgen.web import jobs

        conn = store.get_control_connection(store.control_mysql_config(_mysql_config(job)))
        try:
            per_step, _per_stage = step_estimates(conn, job["steps"], recorded_step_seconds(jobs.configured_jobs_dir()))
        finally:
            conn.close()
    except Exception:  # noqa: BLE001 -- no estimate: the job gets the fixed limit
        return None
    return None if not per_step or any(seconds is None for seconds in per_step) else sum(per_step)


def _run_estimates(job, jobs_dir):
    """
    `(per step, per stage)` for `step_estimates`, or `None`s where there is no control database or the planetGen modules are
    not there: the bar then estimates from the running step alone.
    """
    steps = job["steps"]
    nothing = ([None] * len(steps), [None] * sum(_step_stage_counts(steps)))
    try:
        from planetgen.db import store

        conn = store.get_control_connection(store.control_mysql_config(_mysql_config(job)))
        try:
            return step_estimates(conn, steps, recorded_step_seconds(jobs_dir))
        finally:
            conn.close()
    except Exception:  # noqa: BLE001
        return nothing


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
    env[JOB_ENV] = job["id"]

    steps = job["steps"]
    state["step_estimates"], state["stage_estimates"] = _run_estimates(job, jobs_dir)
    _write_json(state_path, state)
    try:
        records_path = os.path.join(job_dir, RECORDS_FILE)
        env[RECORD_FILE_ENV] = records_path          # every step's command also writes its lines here, as records (UX.122)
        with open(os.path.join(job_dir, "output.log"), "ab", buffering=0) as log, \
                open(records_path, "ab", buffering=0) as records:
            def mark(record, blank_before=False):
                """A stage mark: the display line in the output the page shows and the download holds (the logging standard's
                web form, with its time), the record in the records file."""
                log.write((("\n" if blank_before else "") + logrecord.display_line(record) + "\n").encode("utf-8"))
                records.write((logrecord.file_line(record) + "\n").encode("utf-8"))

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
                count = counts[index - 1]
                mark(logrecord.stage_record(logrecord.STAGE_START, "job", first, total, step["label"], count=count),
                     blank_before=True)
                env[STAGE_OFFSET_ENV] = str(first - 1)
                env[STAGE_TOTAL_ENV] = str(total)
                started = time.time()
                with memwatch.PeakTracker() as watch:
                    current["proc"] = subprocess.Popen(
                        step["argv"], stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                        cwd=job.get("cwd") or None, env=env, **_step_process_options(),
                    )
                    watch.set_group(current["proc"].pid)
                    state["step_pid"] = current["proc"].pid     # so Cancel can end the step if this runner dies (PERF.186)
                    _write_json(state_path, state)
                    asked = None
                    while True:
                        if _cancel_requested():
                            _stop_tree(current["proc"])
                            asked = asked or time.monotonic()
                            if time.monotonic() - asked > CANCEL_KILL_SECONDS:
                                _kill_tree(current["proc"])
                        try:
                            code = current["proc"].wait(timeout=POLL_SECONDS)
                            break
                        except subprocess.TimeoutExpired:
                            pass
                current["proc"] = None
                state.pop("step_pid", None)
                _record_memory(job_dir, log, watch, index, step["label"], first, count)
                took = time.time() - started
                if code == 0:
                    _record_step_seconds(jobs_dir, step["label"], took)
                    mark(logrecord.stage_record(logrecord.STAGE_FINISH, "job", first, total, seconds=took, count=count))
                else:
                    mark(logrecord.stage_record(
                        logrecord.STAGE_FAILED, "job", first, total, seconds=took, count=count,
                        reason="cancelled" if current["cancelled"] else f"exited with status {code}"))
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
        # RQ's time limit (`jobs.job_time_limit`) raises in this thread: the step it was waiting on goes with it.
        _stop_tree(current["proc"])
        _reap(current["proc"])
        state["status"] = "failed"
        state["error"] = (f"The job ran past its time limit and RQ ended it: {exc}"
                          if type(exc).__name__ == "JobTimeoutException" else f"The job runner failed: {exc}")
        _log_traceback(job_dir)
    finally:
        state["finished_at"] = time.time()
        _write_json(state_path, state)
        # The queue goes before the lock: a job that has released its lock
        # is over, queue included (a test or the next job saw the queue
        # a moment after).
        _drop_queue(job["id"])
        _release_lock(job["id"])
        _start_waiting_job(jobs_dir, job_dir)
    return 0 if state["status"] == "succeeded" else 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: python -m planetgen.web.job_runner <job dir>")
    sys.exit(run(sys.argv[1]))
