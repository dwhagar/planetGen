# planetgen/queue/runcontrol.py

"""
What an admin can do to a run in Redis (PERF.99, docs/design/operations/redis-and-the-work-queue.md 5.5 and 5.6): cancel
it, clear one that has stopped answering, delete a finished one, and keep its log. The Queue page and the Generate page
call these; nothing here knows a page.

* `cancel_run` -- sets the run to `cancelling` (the unit that completes a stage then queues no next one, and a unit that
  has not started finds the run not running and does nothing), cancels the stage's waiting units, sends RQ's stop command
  to the started ones, and ends the run `cancelled`, which frees the machine. The job ids follow from the run id, the
  stage and the unit's index (`units.job_id`), so there is no list to keep; the lookups are one pipeline, so a run of
  several thousand units is stopped in about a second (measured on 75: `tests/test_run_control.py`).
* `clear_run` -- the Queue page's "Clear" for a run that shows no sign of life (`runrecord.is_stale`), or whose process
  is gone (`runrecord.holder_is_dead`): the same stop, ending it `failed` with who cleared it.
* `write_log` -- the run's event stream as `<jobs dir>/runs/<run id>.log`, the last `KEEP_LOGS` kept, for the download
  link. Best effort: a log is a nicety.
"""

import os

import rq
from rq import command

from planetgen.queue import redisqueue, runrecord, units, wrapper
from planetgen.util import log

KEEP_LOGS = 20
"""int: How many run logs the jobs directory keeps."""

LOG_DIRECTORY = "runs"
"""str: The folder of run logs inside the jobs directory."""


def _unfinished_units(conn, run):
    """`[(stage, unit)]` of the units of the run that may still be queued or running: the units of the stages that
    have begun and not ended, less those in their done sets."""
    found = []
    for index, stage in enumerate(run["stages"]):
        if stage["started"] is None or stage["ended"] is not None or not stage["units"]:
            continue
        finished = {int(unit) for unit in conn.smembers(runrecord.done_key(run["id"], index))}
        found += [(index, unit) for unit in range(stage["units"]) if unit not in finished]
    return found


def stop_units(conn, run):
    """
    Takes the units of `run` (`runrecord.read_run`'s dict) that are still waiting off the queue and sends the stop
    command to the started ones. Safe whatever a unit is doing: a unit's save is one transaction, and a stopped unit is
    not retried. The work is a few pipelines whatever the number of units: the ids are worked out from the run, the
    started ones are the run's ids in the queue's started registry, and each other id is removed from the queue's list
    (the ones that were there lose their job record).

    Returns:
        dict: `{"cancelled": waiting units removed, "stopped": started units sent the stop command}`.
    """
    queue = redisqueue.queue(wrapper.queue_of(run["kind"]), conn)
    serializer = redisqueue.serializer_for(queue.name)
    ids = [units.job_id(run["id"], stage, unit) for stage, unit in _unfinished_units(conn, run)]
    started = set(queue.started_job_registry.get_job_ids(cleanup=False)) & set(ids)
    pipe = conn.pipeline()
    waiting = [job_id for job_id in ids if job_id not in started]
    for job_id in waiting:
        pipe.lrem(queue.key, 1, job_id)
    removed = [job_id for job_id, count in zip(waiting, pipe.execute()) if count]
    if removed:
        conn.delete(*(rq.job.Job.key_for(job_id) for job_id in removed))
    stopped = 0
    for job_id in sorted(started):
        try:
            command.send_stop_job_command(conn, job_id, serializer)
            stopped += 1
        except (rq.exceptions.InvalidJobOperation, rq.exceptions.NoSuchJobError):
            pass        # it finished between the two reads
    return {"cancelled": len(removed), "stopped": stopped}


def cancel_run(conn, run_id, who=None, now=None):
    """
    Cancels the live run `run_id` now (see the module docstring).

    Returns:
        dict | None: What was stopped (`stop_units`), or `None` when the run is unknown or already over.
    """
    run = runrecord.read_run(conn, run_id)
    if run is None or run["state"] not in runrecord.LIVE_STATES:
        return None
    runrecord.request_cancel(conn, run_id, now)
    stopped = stop_units(conn, run)
    runrecord.finish_run(conn, run_id, "cancelled", error=f"cancelled by {who}" if who else "cancelled", now=now)
    write_log(conn, run_id)
    return stopped


def is_clearable(conn, run, now=None):
    """Whether `run` is live but shows no sign of life: its record is stale, or the process that follows it is gone."""
    return run is not None and run["state"] in runrecord.LIVE_STATES and (
        runrecord.is_stale(run, now) or runrecord.holder_is_dead(conn, run["id"]))


def clear_run(conn, run_id, who=None, now=None):
    """
    Ends a run that has stopped answering as `failed`, with who cleared it, stops what its workers may still be doing
    for it and frees the machine.

    Returns:
        dict | None: What was stopped, or `None` when the run is unknown, over, or still answering (cancel it instead).
    """
    run = runrecord.read_run(conn, run_id)
    if not is_clearable(conn, run, now):
        return None
    stopped = stop_units(conn, run)
    runrecord.finish_run(conn, run_id, "failed", error=f"cleared by {who}: it showed no sign of life" if who
                         else "cleared: it showed no sign of life", now=now)
    write_log(conn, run_id)
    return stopped


def clear_finished(conn):
    """Deletes every finished run's record. Returns how many."""
    cleared = 0
    for run_id in runrecord.recent_runs(conn, limit=10_000):
        cleared += runrecord.delete_run(conn, run_id)
    return cleared


def log_directory():
    """
    The folder of run logs, created if needed (`<jobs dir>/runs`), or `None` when no jobs directory can be written.
    The jobs directory is the web layer's (`planetgen.web.jobs.jobs_dir`), so it is imported only here.
    """
    try:
        from planetgen.web import jobs
        directory = os.path.join(jobs.jobs_dir(), LOG_DIRECTORY)
        os.makedirs(directory, mode=0o750, exist_ok=True)
        return directory
    except (OSError, ImportError):
        return None


def log_path(run_id, directory=None):
    """Where the run's log is kept, or `None` without a jobs directory."""
    directory = directory or log_directory()
    return None if directory is None else os.path.join(directory, f"{run_id}.log")


def write_log(conn, run_id, directory=None):
    """
    Writes the run's event stream to its log file (`runrecord.write_log_file`), keeping the newest `KEEP_LOGS`.
    Never raises: the log is for reading afterwards, never a reason for a run to fail.

    Returns:
        str | None: The path, or `None` when it could not be written.
    """
    directory = directory or log_directory()
    if directory is None:
        return None
    try:
        return runrecord.write_log_file(conn, run_id, directory, keep=KEEP_LOGS)
    except Exception as exc:  # noqa: BLE001 -- a nicety
        log.debug(f"Run log: not written ({exc}).")
        return None


__all__ = ["cancel_run", "clear_finished", "clear_run", "is_clearable", "log_directory", "log_path", "stop_units",
           "write_log"]
