# planetgen/queue/runcontrol.py

"""
What an admin can do to a run in Redis (PERF.99, docs/design/operations/redis-and-the-work-queue.md 5.5 and 5.6): cancel
it, clear one that has stopped answering, delete a finished one, and keep its log. The Queue page and the Generate page
call these; nothing here knows a page.

* `cancel_run` -- (also for a run held for a person, OPS.221) sets the run to `canceling` (a stage that completes then opens no next one), asks the supervisor to stop the run's jobs in the ledger (`controls.stop_run`: each job
  process it forked stops between two pieces and the run's rows end) and ends the run `canceled`, which frees the machine. OPS.235: the RQ queue's units and its stop command are gone.
* `clear_run` -- the Queue page's "Clear" for a run that shows no sign of life (`runrecord.is_stale`), or whose process
  is gone (`runrecord.holder_is_dead`): the same stop, ending it `failed` with who cleared it.
* `write_log` -- the run's event stream as `<jobs dir>/runs/<run id>.log`, the last `KEEP_LOGS` kept, for the download
  link. Best effort: a log is a nicety.
"""

import os

from planetgen.queue import held, runrecord
from planetgen.util import log

KEEP_LOGS = 20
"""int: How many run logs the jobs directory keeps."""

LOG_DIRECTORY = "runs"
"""str: The folder of run logs inside the jobs directory."""


def cancel_run(conn, run_id, who=None, now=None):
    """
    Cancels the live run `run_id` now (see the module docstring).

    Returns:
        dict | None: What was asked of the ledger (`controls.stop_run`; `{}` for a run with no rows there), or `None` when the run is unknown or already over.
    """
    run = runrecord.read_run(conn, run_id)
    if run is None or run["state"] not in runrecord.LIVE_STATES + (runrecord.HELD,):
        return None
    runrecord.request_cancel(conn, run_id, now)         # a held run is not live: it is ended below with nothing to stop
    ended = _end(conn, run, "canceled", f"canceled by {who}" if who else "canceled", now, who=who)
    if run["state"] == runrecord.HELD:
        held.mark_run(conn, run_id, "canceled", f"canceled by {who}" if who else "canceled")      # OPS.221: the ledger's run row follows
    return ended


def is_clearable(conn, run, now=None):
    """Whether `run` is live but shows no sign of life: its record is stale, or the process that follows it is gone."""
    return run is not None and run["state"] in runrecord.LIVE_STATES and (
        runrecord.is_stale(run, now) or runrecord.holder_is_dead(conn, run["id"]))


def clear_run(conn, run_id, who=None, now=None, reason=None):
    """
    Ends a run that has stopped answering as `failed`, with who cleared it (or `reason`, when the clear has a cause of its own), stops what its workers may still be doing
    for it and frees the machine.

    Returns:
        dict | None: What was asked of the ledger (`{}` for a run with no rows there), or `None` when the run is unknown, over, or still answering (cancel it instead).
    """
    run = runrecord.read_run(conn, run_id)
    if not is_clearable(conn, run, now):
        return None
    return _end(conn, run, "failed", reason or (f"cleared by {who}: it showed no sign of life" if who else "cleared: it showed no sign of life"), now, killed=True, who=who)


def _end(conn, run, state, error, now, killed=False, who=None):
    """Stops the run's jobs, ends the run in `state` with `error` and writes its log; returns what `controls.stop_run` did (`{}` for a run with no rows in the ledger). The jobs of a run are the supervisor's to stop
    (OPS.234): it asks each job process it forked to stop between two pieces and ends the run's rows."""
    from planetgen.queue import controls
    stopped = controls.stop_run(conn, run["id"], who) or {}
    runrecord.finish_run(conn, run["id"], state, error=error, now=now, killed=killed, from_held=run["state"] == runrecord.HELD)
    write_log(conn, run["id"])
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


__all__ = ["cancel_run", "clear_finished", "clear_run", "is_clearable", "log_directory", "log_path", "write_log"]
