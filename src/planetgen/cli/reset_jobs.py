#!/usr/bin/env python3
# planetgen.cli.reset_jobs

"""
Resets every job and the queue: for the day a job is stuck somewhere and nothing else will go on.

    python3 -m planetgen.cli.reset_jobs [--yes] [--dry-run] [--mysql-* ...]

It stops what the update guard lists (Generate jobs, work queue runs, queued jobs, registered command-line runs),
kills the queue workers on this host, then removes what is left:

* **Redis**: every RQ key (`rq:*`: the queues, their registries, the jobs, the workers), every run record
  (`planetgen:run:*`, with the active-run keys and event streams), the counters of workers being started and the
  pause flag of the Queue page (the queues are un-paused). The admin's worker count (`planetgen:workers:wanted`)
  and every cache are kept.
* **The control database**: every row of `work_jobs` (their tasks and arguments go with them) and the work lease,
  which is freed.

The galaxy databases and the Generate page's job directories are not touched. `update.sh --reset-jobs` runs this
first. It asks before it changes anything (`--yes` is the answer for a script); `--dry-run` only prints what
there is. It exits 0 when it finished (or has nothing to do) and 1 when it was cancelled or could not reach a
store it was told to reset.
"""

import argparse
import os
import signal
import socket
import sys
import time

import rq

from planetgen.admin import updateguard
from planetgen.db import store
from planetgen.queue import pause, redisqueue, runrecord, workers
from planetgen.util import log

STOP_SECONDS = 20
"""int: How long the polite cancel of the running work is given before the rest is removed by force."""

KEY_PATTERNS = ("rq:*", f"{runrecord.PREFIX}:run:*", workers.STARTING_PREFIX + "*")
"""tuple: The Redis keys removed: RQ's own, the runs' (records, units' results, event streams, the active-run keys),
and the counters of workers being started; the pause flag (PERF.152) is removed with them (`pause.key`)."""


def _redis(url):
    """The Redis client, or `None` when no server answers."""
    try:
        connection = redisqueue.connect(url)
        connection.ping()
        return connection
    except Exception:  # noqa: BLE001 -- no server: nothing is queued
        return None


def _count_queues(connection):
    """`{queue name: (waiting, running)}` for every queue Redis knows."""
    found = {}
    for queue in rq.Queue.all(connection=connection):
        made = redisqueue.queue(queue.name, connection)
        found[queue.name] = (made.count, len(made.started_job_registry.get_job_ids(cleanup=False)))
    return found


def _redis_keys(connection):
    keys = []
    for pattern in (*KEY_PATTERNS, pause.key()):
        keys.extend(connection.scan_iter(match=pattern, count=500))
    return keys


def _local_workers(connection):
    """The registered workers whose process runs on this host: `(worker, pid)`."""
    here = socket.gethostname()
    found = []
    for worker in rq.Worker.all(connection=connection):
        pid = getattr(worker, "pid", None)
        if worker.hostname == here and isinstance(pid, int) and workers.process_exists(pid):
            found.append((worker, pid))
    return found


def _kill_worker(pid, write):
    """Kills the worker and its work horse (a worker started detached leads its own process group)."""
    try:
        if os.getpgid(pid) == pid:
            os.killpg(pid, signal.SIGKILL)
        else:
            os.kill(pid, signal.SIGKILL)
    except ProcessLookupError:
        return
    except OSError as exc:
        write(f"  could not kill worker process {pid}: {exc}")


def _survey_work_rows(control):
    """`(jobs by state, task count, lease holder)` of the control database, or `None` when it has no
    work tables to reset."""
    conn = store.get_control_connection(control, ensure_schema=False)
    try:
        try:
            jobs = {row["state"]: int(row["n"]) for row in conn.execute(
                "SELECT state, COUNT(*) AS n FROM work_jobs GROUP BY state").fetchall()}
            tasks = int(conn.execute("SELECT COUNT(*) AS n FROM work_tasks").fetchone()["n"])
            lease = conn.execute("SELECT holder FROM work_lease WHERE id = 1").fetchone() or {}
        except Exception:  # noqa: BLE001 -- a control schema without the work queue
            return None
        conn.rollback()
        return jobs, tasks, lease.get("holder")
    finally:
        conn.close()


def _clear_work_rows(control):
    conn = store.get_control_connection(control, ensure_schema=False)
    try:
        conn.execute("DELETE FROM work_jobs")
        conn.execute("UPDATE work_lease SET holder = NULL, job_id = NULL, heartbeat_at = NULL "
                     "WHERE id = 1")
        conn.commit()
    finally:
        conn.close()


def reset_jobs(config, redis_url=None, assume_yes=False, dry_run=False, sources=None, write=print, answer=input,
               stop_seconds=STOP_SECONDS, control=None):
    """
    Resets the jobs and the queue (the module's docstring says what that is).

    Args:
        config (MySQLConfig | None): The galaxy database's connection parameters (the control database is next to
            it); `None` leaves the control database alone.
        redis_url (str): The Redis server (default the configured one).
        assume_yes (bool): Do not ask.
        dry_run (bool): Print what there is, change nothing.
        sources: What the update guard looks at (default `updateguard.default_sources(config)`).
        write, answer: Where the lines go and how a question is asked.
        control (MySQLConfig): The control database (default the one next to `config`).

    Returns:
        bool: Whether it finished (a dry run counts); `False` when the operator said no.
    """
    sources = sources if sources is not None else updateguard.default_sources(config)
    control = control or (store.control_mysql_config(config) if config is not None else None)
    connection = _redis(redis_url)
    running = updateguard.list_running(sources)
    queues = _count_queues(connection) if connection is not None else {}
    local = _local_workers(connection) if connection is not None else []
    keys = _redis_keys(connection) if connection is not None else []
    rows = _survey_work_rows(control) if control is not None else None

    write("Running or queued work the update guard sees:" if running else "The update guard sees no running or queued work.")
    for line in updateguard.describe(running):
        write(line)
    if connection is None:
        write("Redis: no server answers, so there are no queues to clear.")
    else:
        for name, (waiting, started) in sorted(queues.items()):
            write(f"Redis queue {name}: {waiting} waiting, {started} running")
        write(f"Redis: {len(keys)} key(s) of the queues, {len(local)} worker process(es) on this host")
    if rows is None:
        write("Control database: no work queue tables to reset." if control is not None else "Control database: left alone.")
    else:
        jobs, tasks, holder = rows
        states = ", ".join(f"{count} {state}" for state, count in sorted(jobs.items())) or "none"
        write(f"Control database: {sum(jobs.values())} job(s) ({states}), {tasks} task(s); the work lease is "
              f"{'held by ' + holder if holder else 'free'}")
    if dry_run:
        write("Dry run: nothing was changed.")
        return True
    if not (running or keys or local or rows and (rows[0] or rows[2])):
        write("Nothing to reset.")
        return True
    if not assume_yes:
        try:
            reply = answer("Stop all of this and remove it? [y/N] ")
        except EOFError:
            reply = ""
        if str(reply).strip().lower() not in ("y", "yes"):
            write("Cancelled: nothing was changed.")
            return False

    if running:
        write("Asking the running work to stop...")
        updateguard.cancel_all(sources, running, write)
        left = updateguard.wait_until_idle(sources, timeout=stop_seconds, poll=0.5, write=write)
        if left:
            write(f"Still there after {int(stop_seconds)} seconds, removed by force below:")
            for line in updateguard.describe(left):
                write(line)
    if connection is not None:
        for worker, pid in _local_workers(connection):
            write(f"Killing worker {worker.name} (process {pid})")
            _kill_worker(pid, write)
        deadline = time.monotonic() + 5
        while _local_workers(connection) and time.monotonic() < deadline:
            time.sleep(0.2)
        keys = _redis_keys(connection)
        for start in range(0, len(keys), 500):
            connection.delete(*keys[start:start + 500])
        write(f"Redis: removed {len(keys)} key(s)")
    if rows is not None:
        _clear_work_rows(control)
        write(f"Control database: removed {sum(rows[0].values())} job(s) and {rows[1]} task(s), freed the work lease")
    write("Jobs and queue reset. Anything that was cancelled has to be started again.")
    return True


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--yes", "-y", action="store_true", help="Do not ask.")
    parser.add_argument("--dry-run", action="store_true", help="Print what there is, change nothing.")
    store.add_mysql_connection_args(parser)
    log.add_debug_arguments(parser)
    args = parser.parse_args(argv)
    log.apply_debug_arguments(args, parser)
    config = store.mysql_config_from_args(args)
    if not args.yes and not args.dry_run and not sys.stdin.isatty():
        print("error: there is no terminal to ask on; add --yes to reset without asking.", file=sys.stderr)
        return 1
    try:
        done = reset_jobs(config, assume_yes=args.yes, dry_run=args.dry_run)
    except Exception as exc:  # noqa: BLE001 -- said plainly; nothing half-reported
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0 if done else 1


if __name__ == "__main__":
    sys.exit(main())
