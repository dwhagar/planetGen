# planetgen/queue/redisqueue.py

"""
The work queue's Redis (PERF.24, OPS.235): the connection, the JSON serializer the records use, the time limits of a job, and the detached process the supervisor starts.

RQ is gone (OPS.235): nothing queues a job on Redis for a worker to take. A run's units and a page's steps are rows of the job ledger that the supervisor's launcher forks; the queue names
below only name the lanes of those rows.
"""

import datetime
import decimal
import json
import math
import os
import subprocess
import threading

from planetgen.queue import redisclient
from planetgen.util import log

GENERATION_QUEUE = "planetgen-generation"
"""str: The queue of the first design (PERF.24), where the supervising process of a run put its tasks. Nothing is queued on
it since the rebuild: a run's units go to `BULK_QUEUE`. The name stays in `QUEUES`, which keeps the leftovers
of older releases visible to the reset and stop-all commands."""

BULK_QUEUE = "planetgen-bulk"
"""str: A run's units of work (PERF.95): the scatters, sector fills and backfill blocks of the work queue
rebuild. Its workers serve this queue and no other."""

MAINTENANCE_QUEUE = "planetgen-maintenance"
"""str: The units of a `MAINTENANCE` run (the cache warm-up) that has no ledger (PERF.150; OPS.229: a run with a ledger has its units handed out by the launcher). Its workers have a cap of their own
(`workers.maintenance_cap`)."""

CHECK_QUEUE = "planetgen-maintenance-check"
"""str: The units of a `CHECK` run (the database check, PERF.245) that has no ledger (OPS.229). One of the maintenance lane's task types: it has a queue of its own, so the thousands of units a cache
warm-up queued ahead of it never make it wait."""

def maintenance_lane():
    """
    The queues of the maintenance lane, one per task type (PERF.245): a `MAINTENANCE` run's units (the cache warm-up) and a `CHECK` run's units. A run with a ledger does not use them (OPS.229: its units are
    forked jobs of the supervisor's launcher, class `maintenance`); only a run with no ledger does (OPS.271), each queue with workers of its own (`workers.lane_of`). A function, not a constant, so a test that gives
    `MAINTENANCE_QUEUE` a name of its own gets a lane of its own.
    """
    return (MAINTENANCE_QUEUE, CHECK_QUEUE)


QUEUES = (GENERATION_QUEUE, BULK_QUEUE, MAINTENANCE_QUEUE, CHECK_QUEUE)
"""tuple: Every shared queue."""

def redis_url():
    """The Redis server: `PLANETGEN_REDIS_URL`, else `config.json`'s
    `redis.url` (the settings model applies both)."""
    return redisclient.default_url()


def connect(url=None, **options):
    """The shared client for `url` (default `redis_url()`); `options` are `redisclient.client`'s (PERF.94)."""
    return redisclient.client(url or redis_url(), **options)


def _json_default(value):
    """What `JsonSerializer` writes for a value `json` has no form for: a `Decimal` as the number it holds, a
    time as ISO text."""
    if isinstance(value, decimal.Decimal):
        return int(value) if value == value.to_integral_value() else float(value)
    if isinstance(value, (datetime.datetime, datetime.date)):
        return value.isoformat()
    raise TypeError(f"{type(value).__name__} can't be sent to the queue as JSON")


class JsonSerializer:
    """The serializer of every record the queue keeps in Redis (PERF.94, SEC.39, finding F8): arguments and results travel as JSON, not
    pickle, so a process that can reach Redis can write data but not an object that runs code when it is read."""

    @staticmethod
    def dumps(value, *args, **kwargs):
        return json.dumps(value, default=_json_default, *args, **kwargs).encode("utf-8")

    @staticmethod
    def loads(raw, *args, **kwargs):
        return json.loads(raw.decode("utf-8"), *args, **kwargs)


TIMEOUT_MULTIPLE = 10
"""int: How many times its expected seconds a job may run before RQ ends it (PERF.106, PERF.144): work is expected
to take that long on all its workers, so a job that takes ten times as long is hung or starved, not slow."""

MIN_TIMEOUT_SECONDS = 15 * 60
"""int: The shortest limit a job gets, so work that is usually quick is not ended by one slow start-up."""

MAX_TIMEOUT_SECONDS = 6 * 3600
"""int: The longest limit a job gets from `job_timeout`'s defaults, and the one of a job whose work has no
estimate at all (a fixed ceiling, not an estimate); a site names a different floor, ceiling or `unestimated`
limit where that fits its work."""


def job_timeout(expected_seconds=None, ceiling=MAX_TIMEOUT_SECONDS, floor=MIN_TIMEOUT_SECONDS, unestimated=None):
    """
    How long a queued job may run before RQ ends it (PERF.144): `TIMEOUT_MULTIPLE` times the seconds the system
    expects it to take (recorded timings, then the stage ETAs, then the size of the work: each site's own best
    estimate), but at least `floor` and at most `ceiling`. Never -1 (no limit).

    Args:
        expected_seconds (float, optional): The best estimate of the job's run time; `None`, zero, negative or not
            finite counts as no estimate.
        ceiling (int): The most seconds.
        floor (int): The fewest seconds.
        unestimated (int, optional): The limit of a job with no estimate, a fixed ceiling and not an estimate
            (default `ceiling`).

    Returns:
        int: Whole seconds, at least 1.
    """
    if expected_seconds is None or not math.isfinite(expected_seconds) or expected_seconds <= 0:
        return max(int(ceiling if unestimated is None else unestimated), 1)
    return max(int(min(max(TIMEOUT_MULTIPLE * expected_seconds, floor), ceiling)), 1)


def start_detached(argv, cwd, on_exit=None):
    """
    Starts `argv` apart from this process, in its own session, and reaps
    it when it exits, so it doesn't linger as a zombie of the web server
    process; if this process goes first, init adopts it instead.

    Args:
        on_exit (callable, optional): Called with the exit status when the process ends, from the thread that reaps it: the notice of its exit (OPS.251), so a claim
            made for a process that is starting ends with it and not at a time limit.

    Returns:
        subprocess.Popen: The started process.
    """
    proc = subprocess.Popen(
        argv, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        cwd=cwd, close_fds=True, start_new_session=True,
    )
    threading.Thread(target=_reap, args=(proc, argv, on_exit), name=f"planetgen-worker-{proc.pid}", daemon=True).start()
    return proc


def _reap(proc, argv, on_exit=None):
    """Waits for a detached worker and, as its parent, writes down how it ended when that was badly (`log.process_ended`: a crash, a kill by
    the out-of-memory killer or `kill -9`): its own stderr goes nowhere, so nobody else would. A refusal during an update (`EX_TEMPFAIL`) is not one."""
    status = proc.wait()
    try:
        log.process_ended("queue", f"Worker process {proc.pid} ({' '.join(argv[3:]) if len(argv) > 3 else 'worker'})", status,
                          ignore=(0, os.EX_TEMPFAIL), pid=proc.pid)
    except Exception as exc:  # noqa: BLE001 -- a log line is never a reason to lose the reaper thread
        log.swallowed("queue/redisqueue", "_reap", exc)
    if on_exit is not None:
        try:
            on_exit(status)
        except Exception as exc:  # noqa: BLE001 -- a claim that cannot be released here is found stale by the next caller (its holder's process is gone)
            log.swallowed("queue/redisqueue", "_reap on_exit", exc)


class Unavailable(Exception):
    """No Redis server answers at the configured URL."""

