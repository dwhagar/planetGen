# planetgen/web/job_runner.py

"""
Runs one step of a background job the web interface's admin Generate page started (OPS.244, slice B).

A job is a list of steps, each one command line (`planetgen.cli.reset`, then `planetgen plan`, then `planetgen galaxy ...`), and a run of the job ledger (`queue/pagerun.py`): `planetgen.web.jobs.start_job`
writes a row for each step and the supervisor's launcher forks one job process for the step that is `ready` (`queue/launcher.py`, `queue/forkedjob.py`), which calls `run_step`. The step's command line is run as
a process of its own, in its own session, its output appended to `<job dir>/output.log`; the job's log as records goes to `<job dir>/records.log` and the step's progress to `<job dir>/progress.json`
(`PLANETGEN_PROGRESS_FILE`), where the page reads them. What the page shows of the job (which step, how it ended) is read from the ledger's rows.

A step that ends well makes the next one `ready` (`pagerun.step_done`); the first that fails stops the job (`pagerun.step_failed`). The page's Cancel is a request to the supervisor (`controls.stop_engine_run`) that
ends the rows and sends this process the `stop` request on its pipe (with none alive it ends the rows, which this process reads every `ROW_CHECK_SECONDS`); `run_step` then stops the step and every process it started (SIGTERM to the step's process group, SIGKILL after `CANCEL_KILL_SECONDS`).
The process group of a step is kept in `step.pid` while it runs, so Cancel can end a step whose runner is gone (`jobs.cancel_job`). When the job is over its lock in Redis is freed
(`jobs.release_lock`) and the job waiting behind it is started (`jobs.promote_queued`).

How long each step took is kept by its label in the jobs directory (`STEP_TIMES_NAME`), for the next job's estimates.
"""

import json
import os
import signal
import subprocess
import sys
import tempfile
import time

from planetgen.util import log as planetgen_log
from planetgen.util import logrecord, memwatch

STEP_PID_NAME = "step.pid"
"""str: The file in a job's directory that holds the process group of the step that runs, for as long as it runs."""

STEP_TIMES_NAME = "step-times.json"
"""str: The file in the jobs directory with how long each kind of step took in the jobs that succeeded (PERF.33): `{label:
[how many times, total seconds]}`, read for the next job's estimates (`recorded_step_seconds`)."""

MEMORY_NAME = "memory.json"
"""str: The file in a job's directory with the memory peak of each step that ran (PERF.183): a list of
`{"step", "label", "first_stage", "stages", "seconds", "peaks"}` (`memwatch.PeakTracker.peaks`)."""

CANCEL_KILL_SECONDS = 5
"""float: How long a step that was asked to stop (SIGTERM to its process group) may keep running before the whole group is killed (SIGKILL, PERF.186). Less than the supervisor's grace
(`controls.GRACE_SECONDS`), after which it kills this process, and the step's group with it is then no one's."""

POLL_SECONDS = 0.25
"""float: How often a running step is checked for a request to stop."""

ROW_CHECK_SECONDS = 1.0
"""float: How often the step's row is read to see whether it was canceled: a Cancel with no supervisor alive ends the rows and has no pipe to send this process the request."""

STOP_WAIT_SECONDS = 10
"""float: How long a step that was stopped gets to exit before it is killed (`_reap`)."""

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


class NoLink:
    """The pipe of a step that runs in the caller's own process (`run_to_end`, a test): nobody asks it to stop."""

    gone = False
    stop_requested = False

    def requests(self):
        return []

    def send(self, message):
        return None


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


def _signal_tree(proc, sig):
    """Sends `sig` to a running step's process group (the step and every process it started); nothing for a step that is gone."""
    if proc is None or proc.poll() is not None:
        return
    try:
        os.killpg(proc.pid, sig)
    except OSError:
        pass


def _row_is_canceled(control, row_id):
    """
    Whether the job's `engine_jobs` row says Cancel. A row that is not there is no Cancel: nobody asked to stop the job, so it runs on and the runner does not stop a step that may still be
    working. (OPS.296: the runner read a field of no row, failed with a `TypeError` and stopped the step it was following; why a row went missing there is not known: the control schema is
    a database of its own that a reset does not touch.)
    """
    row = control.execute("SELECT state FROM engine_jobs WHERE id = ?", (row_id,)).fetchone()
    return row is not None and row["state"] == "canceled"


def _stop_tree(proc):
    """Stops a step and every process it started."""
    _signal_tree(proc, signal.SIGTERM)


def _kill_tree(proc):
    """Kills a step and every process in its group, for a step that ignored the request to stop."""
    _signal_tree(proc, signal.SIGKILL)


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
    except Exception as exc:  # noqa: BLE001 -- a job is never lost to its memory note
        planetgen_log.swallowed("web/job_runner", "_record_memory", exc)


def _release_lock(job_id):
    """Frees the job lock if it names this job (best effort; the next job clears a stale one)."""
    try:
        from planetgen.web import jobs

        jobs.release_lock(job_id)
    except Exception as exc:  # noqa: BLE001 -- a lock left behind is cleared by the next job
        planetgen_log.swallowed("web/job_runner", "_release_lock", exc)


def _start_waiting_job(jobs_dir, job_dir):
    """Starts the job waiting behind this one (ADM.123, best effort: an admin page asking for the running job starts
    it too)."""
    try:
        from planetgen.web import jobs

        jobs.promote_queued(jobs_dir)
    except Exception as exc:  # noqa: BLE001 -- the job stays queued
        planetgen_log.swallowed("web/job_runner", "_start_waiting_job", exc)
        _log_traceback(job_dir)


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
    except (SystemExit, Exception) as exc:  # noqa: BLE001 -- the step then uses its recorded time
        planetgen_log.swallowed("web/job_runner", "_staged_command", exc)
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
    except (SystemExit, Exception) as exc:  # noqa: BLE001 -- the step then uses its recorded time
        planetgen_log.swallowed("web/job_runner", "_job_entries", exc)
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


def _step_output(path, start, limit=8192):
    """The last `limit` bytes a step wrote to the job's output file after offset `start`, as text (empty when it cannot be read)."""
    try:
        with open(path, "rb") as handle:
            end = handle.seek(0, os.SEEK_END)
            handle.seek(max(start, end - limit))
            return handle.read().decode("utf-8", errors="replace")
    except OSError:
        return ""


def _log_traceback(job_dir):
    """Appends the traceback of the exception being handled to the job's
    output log, so the page's log shows the runner's own failure in full
    (ADM.25). Never raises."""
    try:
        import traceback

        with open(os.path.join(job_dir, "output.log"), "ab", buffering=0) as log:
            log.write(("\n" + traceback.format_exc()).encode("utf-8"))
    except Exception as exc:  # noqa: BLE001 -- the error is in the ledger regardless
        planetgen_log.swallowed("web/job_runner", "_log_traceback", exc)


def _mysql_config(job):
    """The database the job writes (`job["database"]`) on the configured server, as planetGen's `MySQLConfig`."""
    from planetgen.db import store

    return store.MySQLConfig(database=job.get("database"))


def _step_env(job, job_dir, progress_path):
    """The environment of a step: this process's, the database the job writes (`jobs.mysql_env`, built from the site's settings when the step starts, never stored), and the files the step writes to."""
    from planetgen.web import jobs

    env = dict(os.environ)
    if job.get("database"):
        env.update(jobs.mysql_env(_mysql_config(job), job["database"]))
    env["PYTHONUNBUFFERED"] = "1"
    env["PLANETGEN_PROGRESS_FILE"] = progress_path
    env[JOB_ENV] = job["id"]
    env[RECORD_FILE_ENV] = os.path.join(job_dir, RECORDS_FILE)      # every step's command also writes its lines here, as records (UX.122)
    return env


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
    except Exception as exc:  # noqa: BLE001 -- the job runs without the line's details
        planetgen_log.swallowed("web/job_runner", "_run_line", exc)
        return f"Galaxy seed unknown, PlanetGen unknown, run: {run}"
    seed = None
    try:
        conn = store.get_connection(_mysql_config(job), ensure_schema=False)
        try:
            seed = store.get_galaxy_seed(conn)
        finally:
            conn.close()
    except Exception as exc:  # noqa: BLE001 -- no database or no galaxy yet
        planetgen_log.swallowed("web/job_runner", "_run_line", exc)
    return version_key.run_line(seed, run)


def _run_estimates(control, job, jobs_dir):
    """
    `(per step, per stage)` for `step_estimates`, or `None`s where the stage times cannot be read: the bar then estimates from the running step alone.
    """
    steps = job["steps"]
    nothing = ([None] * len(steps), [None] * sum(_step_stage_counts(steps)))
    try:
        return step_estimates(control, steps, recorded_step_seconds(jobs_dir))
    except Exception as exc:  # noqa: BLE001
        planetgen_log.swallowed("web/job_runner", "_run_estimates", exc)
        return nothing


def _say_failed(job_dir, job, number, reason, started):
    """
    Writes into the job's log, and its records, that step `number` failed and why (OPS.255): a runner that failed with nothing in the log left a job that went silent. The stage the job was at is the one that failed.
    """
    try:
        counts = _step_stage_counts(job["steps"])
        first = sum(counts[:number - 1]) + 1
        record = logrecord.stage_record(logrecord.STAGE_FAILED, "job", first, sum(counts), seconds=time.time() - started, reason=reason, count=counts[number - 1])
        with open(os.path.join(job_dir, "output.log"), "ab", buffering=0) as log:
            log.write(("\n" + logrecord.display_line(record) + "\n").encode("utf-8"))
        with open(os.path.join(job_dir, RECORDS_FILE), "ab", buffering=0) as records:
            records.write((logrecord.file_line(record) + "\n").encode("utf-8"))
    except Exception as exc:  # noqa: BLE001 -- the ledger already says it; the line is a courtesy to the log
        planetgen_log.swallowed("web/job_runner", "_say_failed", exc)


def _forget(job_id):
    """Takes the job off the supervisor's keep-alive list (`supervisor.pages_key`). Best effort."""
    try:
        from planetgen.queue import redisqueue, supervisor

        redisqueue.connect().hdel(supervisor.pages_key(), job_id)
    except Exception as exc:  # noqa: BLE001 -- `pagerun.tidy` forgets it at the supervisor's next look
        planetgen_log.swallowed("web/job_runner", "_forget", exc)


def _write_step_pid(job_dir, pid):
    """Keeps the process group of the running step (`pid`, or `None` once it is over) where Cancel finds it if this process is gone."""
    path = os.path.join(job_dir, STEP_PID_NAME)
    try:
        if pid is None:
            os.remove(path)
        else:
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(str(pid))
    except OSError:
        pass


def run_step(spec, link):
    """
    Runs the step of a Generate page job that the launcher's `spec` (`{"row", "run"}`: the step's row of the ledger and its run) names, as the job process the launcher forked, and returns the exit code of
    that process: 0 when the step ended well or was stopped, `template.EXIT_FAILED` when it failed. `link` is the process's pipe to the supervisor (`NoLink` where there is none).
    """
    from planetgen.queue import pagerun

    control = pagerun.connect()
    try:
        return _run_step(control, spec["run"], spec["row"], link)
    finally:
        control.close()


def _run_step(control, run, row_id, link):
    """The body of `run_step` once the ledger is open."""
    from planetgen.queue import controls, pagerun, template

    row = control.execute("SELECT item_key FROM engine_jobs WHERE id = ?", (row_id,)).fetchone()
    job = pagerun.job_of(control, run)
    if row is None or job is None:
        raise LookupError(f"step row {row_id} of ledger run {run} is not a step of a Generate page job")
    number = int(row["item_key"])
    steps = job["steps"]
    step = steps[number - 1]
    jobs_dir = job["root"]
    job_dir = os.path.join(jobs_dir, job["id"])
    os.makedirs(job_dir, mode=0o750, exist_ok=True)
    progress_path = os.path.join(job_dir, "progress.json")
    outcome = {"canceled": False, "failed": None, "proc": None, "peak": None, "finished": False}
    started = time.time()
    try:
        if number == 1:
            pagerun.begin(control, run)
            per_step, per_stage = _run_estimates(control, job, jobs_dir)
            pagerun.set_estimates(control, run, per_step, per_stage)
        env = _step_env(job, job_dir, progress_path)
        counts = _step_stage_counts(steps)
        total = sum(counts)
        first = sum(counts[:number - 1]) + 1
        count = counts[number - 1]
        with open(os.path.join(job_dir, "output.log"), "ab", buffering=0) as log, \
                open(os.path.join(job_dir, RECORDS_FILE), "ab", buffering=0) as records:
            def mark(record, blank_before=False):
                """A stage mark: the display line in the output the page shows and the download holds (the logging standard's
                web form, with its time), the record in the records file."""
                log.write((("\n" if blank_before else "") + logrecord.display_line(record) + "\n").encode("utf-8"))
                records.write((logrecord.file_line(record) + "\n").encode("utf-8"))

            if number == 1:
                log.write((_run_line(job) + "\n").encode("utf-8"))
            try:
                os.remove(progress_path)
            except OSError:
                pass
            mark(logrecord.stage_record(logrecord.STAGE_START, "job", first, total, step["label"], count=count), blank_before=True)
            env[STAGE_OFFSET_ENV] = str(first - 1)
            env[STAGE_TOTAL_ENV] = str(total)
            began = time.time()
            output_start = os.fstat(log.fileno()).st_size
            with memwatch.PeakTracker() as watch:
                proc = subprocess.Popen(step["argv"], stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, cwd=job.get("cwd") or None, env=env, **_step_process_options())
                outcome["proc"] = proc
                watch.set_group(proc.pid)
                _write_step_pid(job_dir, proc.pid)          # so Cancel can end the step if this process dies (PERF.186)
                asked = None
                looked = 0.0
                while True:
                    link.requests()
                    if link.stop_requested:
                        outcome["canceled"] = True
                    elif time.monotonic() - looked >= ROW_CHECK_SECONDS:         # a Cancel ends the rows, and reaches this process by the supervisor's request when one is alive
                        looked = time.monotonic()
                        control.commit()
                        if _row_is_canceled(control, row_id):
                            outcome["canceled"] = True
                    if outcome["canceled"]:
                        _stop_tree(proc)
                        asked = asked or time.monotonic()
                        if time.monotonic() - asked > CANCEL_KILL_SECONDS:
                            _kill_tree(proc)
                    try:
                        code = proc.wait(timeout=POLL_SECONDS)
                        break
                    except subprocess.TimeoutExpired:
                        pass
            outcome["proc"] = None
            _write_step_pid(job_dir, None)
            outcome["peak"] = watch.peaks.get("step")
            _record_memory(job_dir, log, watch, number, step["label"], first, count)
            took = time.time() - began
            if code == 0:
                _record_step_seconds(jobs_dir, step["label"], took)
                mark(logrecord.stage_record(logrecord.STAGE_FINISH, "job", first, total, seconds=took, count=count))
            else:
                mark(logrecord.stage_record(logrecord.STAGE_FAILED, "job", first, total, seconds=took, count=count,
                                            reason="canceled" if outcome["canceled"] else f"exited with status {code}"))
                if not outcome["canceled"]:
                    planetgen_log.process_ended("job", f"Step {step['label']!r} of web job {job['id']}", code,
                                                output=_step_output(os.path.join(job_dir, "output.log"), output_start), job=job["id"])
        if outcome["canceled"]:
            controls.end_rows(control, run)             # a Cancel ends the rows before the step is stopped; a stop with nobody to end them ends them here
            control.commit()
            outcome["finished"] = True
            return 0
        if code == 0:
            outcome["finished"] = pagerun.step_done(control, row_id, run, peak_mb=outcome["peak"])
            return 0
        outcome["failed"] = f"{step['label']} exited with status {code}."
        pagerun.step_failed(control, row_id, run, outcome["failed"], peak_mb=outcome["peak"])
        outcome["finished"] = True
        return template.EXIT_FAILED
    except Exception as exc:  # noqa: BLE001 -- recorded for the page, never lost
        planetgen_log.swallowed("web/job_runner", "_run_step", exc)
        _stop_tree(outcome["proc"])
        _reap(outcome["proc"])
        _write_step_pid(job_dir, None)
        text = f"The job runner failed: {type(exc).__name__}: {exc}"
        _say_failed(job_dir, job, number, text, started)
        _log_traceback(job_dir)
        planetgen_log.error("The job runner failed (job %s): %s", job["id"], exc, exc_info=True)
        try:
            control.rollback()
            pagerun.step_failed(control, row_id, run, text)
        except Exception as again:  # noqa: BLE001 -- the launcher settles a row this process left `running`
            planetgen_log.swallowed("web/job_runner", "_run_step", again)
        outcome["finished"] = True
        return template.EXIT_FAILED
    finally:
        if outcome["finished"]:
            _forget(job["id"])
            _release_lock(job["id"])
            _start_waiting_job(jobs_dir, job_dir)


def run_to_end(job_id, jobs_dir=None, link=None):
    """
    Runs the steps of the job `job_id` one after the other in this process, as the launcher would hand them out (each one's attempt begun here), until the job is over, and returns 0 when every step succeeded,
    else 1. For the one place that has no supervisor (`planetgen smoke`) and the tests.

    Raises:
        LookupError: The ledger has no such job.
    """
    from planetgen.queue import ledger, pagerun

    link = link or NoLink()
    control = pagerun.connect()
    try:
        run = pagerun.run_of(control, job_id)
        if run is None:
            raise LookupError(f"no Generate page job {job_id!r} in the ledger")
        while True:
            row = control.execute("SELECT id FROM engine_jobs WHERE run_id = ? AND command = ? AND state = 'ready' ORDER BY item_key LIMIT 1", (run, pagerun.STEP)).fetchone()
            if row is None:
                break
            ledger.begin_attempt(control, row["id"], pid=os.getpid())
            if _run_step(control, run, row["id"], link):
                break
        states = [found["state"] for found in control.execute("SELECT state FROM engine_jobs WHERE run_id = ? AND command = ?", (run, pagerun.STEP)).fetchall()]
        control.commit()
        return 0 if all(state == "done" for state in states) else 1
    finally:
        control.close()
