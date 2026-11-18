# planetgen/queue/failures.py

"""
Failure by kind (OPS.220, L6 of `docs/design/operations/job-engine.md`, design section 3.7).

Every way a unit can fail is one of four kinds, and each kind has one policy:

| kind | what it is | what happens |
|---|---|---|
| `passing` | a MySQL deadlock (1213), lock wait timeout (1205) or lost connection (2006, 2013, a command on a closed connection), a Redis connection that failed | the engine runs the unit again after a wait that doubles from `FIRST_WAIT_SECONDS` to `MAX_WAIT_SECONDS`, spread at random, inside the run's budget (`charge`) |
| `killed` | the process was killed (signal 9, out of memory) or timed out | the launcher hands it out again alone (`tokens.ALONE`) with a claim `memorygate.GROW_FACTOR` times as large (`Launcher._grow`, `Launcher.claim`; OPS.236) |
| `idle` | the unit showed no sign of life for longer than its command ever has | the existing path: `unitstuck` dumps its stack, stops it and runs it once more alone |
| `lasting` | any other exception | not retried; the run says so (stage 3's held items, OPS.221, hold the unit and carry on) |

`classify` is the one place that says which kind an exception is, so the three paths that name it (the unit's own `except`, RQ's failure callback and the run's error text) cannot disagree. A passing
failure the budget no longer allows is reported as it happened, with the budget said in the text, and is then as final as a lasting one.

The budget is per run: `RETRIES_PER_RUN` runs again and `WAIT_PER_RUN_SECONDS` of waiting, kept in one Redis hash (`runrecord.failures_key`), so a database that is sick for an hour is not asked the same
question by 500 units with a new wait each. The waits are drawn from this module's own random stream, never from a generation's, so waiting moves no star (GEN.39, GEN.56).
"""

import random

import pymysql
import redis

from planetgen.queue import runrecord
from planetgen.util.format import format_duration_seconds

PASSING, KILLED, IDLE, LASTING = "passing", "killed", "idle", "lasting"
KINDS = (PASSING, KILLED, IDLE, LASTING)

PASSING_MYSQL_ERRORS = (1213, 1205, 2006, 2013)
"""tuple: The MySQL errors that can pass by themselves: a deadlock, a lock wait timeout, the server gone away and the connection lost during a query (the first two are what `db.contention` retries inside a transaction)."""

PASSING_REDIS_ERRORS = (redis.exceptions.ConnectionError, redis.exceptions.TimeoutError, redis.exceptions.BusyLoadingError)
"""tuple: The Redis errors that can pass by themselves: no connection, no answer in time, the server still loading its data."""

KILLED_EXCEPTIONS = ("JobTimeoutException", "WorkerLostError", "AbandonedJobError")
"""tuple: The names of the exceptions RQ gives a job whose process was killed, timed out or lost."""

IDLE_EXCEPTIONS = ("UnitStuck",)
"""tuple: The names of the exceptions that say the unit was stopped for showing no sign of life (`unitstuck`)."""

FIRST_WAIT_SECONDS = 2.0
MAX_WAIT_SECONDS = 300.0
"""float: The wait before the first retry of a unit, and the longest wait there is: the waits double from the first to the cap."""

SPREAD = (0.5, 1.0)
"""tuple: A wait is a random share of its ceiling between these two, so units that failed together do not come back together."""

RETRIES_PER_RUN = 40
"""int: How many times in all the units of one run may be run again for a passing failure."""

WAIT_PER_RUN_SECONDS = 1800.0
"""float: How long in all the units of one run may wait to be run again."""

_rng = random.Random()
"""The waits' own stream, so waiting never moves a draw of the generator."""


def classify(exc):
    """The kind of failure `exc` (an exception, or the name of one as text) is: `PASSING`, `KILLED`, `IDLE` or `LASTING`."""
    name = exc if isinstance(exc, str) else type(exc).__name__
    if name in KILLED_EXCEPTIONS:
        return KILLED
    if name in IDLE_EXCEPTIONS:
        return IDLE
    if isinstance(exc, PASSING_REDIS_ERRORS):
        return PASSING
    if isinstance(exc, pymysql.err.InterfaceError):
        return PASSING
    if isinstance(exc, pymysql.err.OperationalError) and exc.args and exc.args[0] in PASSING_MYSQL_ERRORS:
        return PASSING
    return LASTING


def wait_seconds(retry, rng=None):
    """The wait before retry number `retry` (1 for the first) of one unit: its ceiling doubles from `FIRST_WAIT_SECONDS` to `MAX_WAIT_SECONDS`, and the wait is a random share of it (`SPREAD`)."""
    ceiling = min(FIRST_WAIT_SECONDS * 2.0 ** (max(retry, 1) - 1), MAX_WAIT_SECONDS)
    low, high = SPREAD
    return ceiling * (rng or _rng).uniform(low, high)


def charge(conn, run_id, wait):
    """
    Takes one retry and `wait` seconds from the run's budget.

    Returns:
        bool: `False` when the run has no retry or no waiting left (nothing is taken then).
    """
    key = runrecord.failures_key(run_id)
    used = conn.hmget(key, "retries", "waited")
    retries, waited = int(runrecord._s(used[0]) or 0), float(runrecord._s(used[1]) or 0.0)
    if retries >= RETRIES_PER_RUN or waited + wait > WAIT_PER_RUN_SECONDS:
        return False
    conn.hincrby(key, "retries", 1)
    conn.hincrbyfloat(key, "waited", wait)
    conn.expire(key, runrecord.KEEP_SECONDS)
    return True


def retry_wait(conn, run_id, exc, retries, rng=None):
    """
    What the engine does about the failure `exc` of a unit that has been run again `retries` times: the seconds to wait before it is run again, or `None` for a failure that is not to be retried (not a
    passing one, or the run's budget is used up).
    """
    if not isinstance(exc, Exception) or classify(exc) != PASSING:
        return None
    wait = wait_seconds(retries + 1, rng)
    return wait if charge(conn, run_id, wait) else None


def spent(conn, run_id):
    """`(retries, seconds waited)` the run has taken from its budget so far."""
    used = conn.hmget(runrecord.failures_key(run_id), "retries", "waited")
    return int(runrecord._s(used[0]) or 0), float(runrecord._s(used[1]) or 0.0)


def describe(conn, run_id, exc):
    """What a passing failure that ended the run adds to the run's error text (what the engine had already tried), as ` (...)`; empty for any other kind, whose text is the exception's own."""
    if classify(exc) != PASSING:
        return ""
    retries, waited = spent(conn, run_id)
    return f" (a passing failure the engine could not clear: units of the run were run again {retries} time(s) after {format_duration_seconds(waited)} of waiting)"
