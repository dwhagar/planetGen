"""
planetgen/db/contention.py

Retrying a database transaction that lost a deadlock (MySQL 1213) or timed out on a row lock (1205), without hammering the database, and
asking for the database to itself when plain retrying does not get through (DB.113).

* **Dynamic waits.** A retry waits a random time (decorrelated jitter: between the base and three times the last wait, capped) that is scaled
  by how much trouble this process has had lately (`pressure`: a moving average of its recent attempts that failed), so a deadlock storm
  spreads its writers out instead of having them all come back together.
* **The need.** After `attempts` failed tries the unit raises a `Need` ("exclusive": the contended step to myself) to the gate. The gate grants it
  by taking a MySQL named lock that every other writer of the same gate checks before each of its attempts, so they hold off (and say so) until the
  unit is done. The unit then makes one more try alone. The grant is released when that try ends, either way.
* **The log.** Each retry (error, attempt, wait), each need raised, each grant, wait for a grant and release is logged.

A stage uses it by naming its `gate` (the contended step) and giving `attempt(conn)`, one whole transaction that commits at its end:
`run_with_retries(conn, attempt, gate="scatter-ids", label="Massive stars, shells 5-6")`.
"""
import os
import random
import time

import pymysql

from planetgen.db import connection_loss, store
from planetgen.util import log

RETRYABLE_ERRORS = (1213, 1205)
"""tuple: The MySQL errors that are retried: a deadlock, and a lock wait timeout."""

ATTEMPTS = 8
"""int: Plain tries before the unit raises its `Need`; the try that follows is the granted one."""

BASE_WAIT_S = 0.05
MAX_WAIT_S = 5.0
"""float: The shortest random wait and the longest one a retry can have."""

GRANT_TIMEOUT_S = 120
"""int: How long a need waits for the gate (another unit may hold it) before the unit gives up with the lock error."""

_rng = random.Random()
"""The waits' own stream (not a generation draw), so waiting never moves a unit's draws (GEN.39, GEN.56)."""

_pressure = 0.0
"""float: This process's moving average of failed attempts (0 to 1): see `pressure`."""


def pressure():
    """How much trouble this process has had lately, 0 (none) to 1 (every attempt fails)."""
    return _pressure


def _note(failed):
    global _pressure
    _pressure = 0.8 * _pressure + 0.2 * (1.0 if failed else 0.0)


def next_wait(attempt, last_wait):
    """
    The wait before the retry that follows failed try number `attempt`: random between the base and three times the last wait (the first
    wait has the base to work from), made longer by `pressure()` (up to 3 times), and never above `MAX_WAIT_S`.

    Args:
        attempt (int): The try that just failed, from 1.
        last_wait (float): The wait before it, seconds (0 for the first retry).
    """
    top = max(BASE_WAIT_S * 2, last_wait * 3)
    wait = _rng.uniform(BASE_WAIT_S, top) * (1.0 + 2.0 * _pressure) * (1 + os.getpid() % 4 * 0.05)
    return min(wait, MAX_WAIT_S)


class Need:
    """What a unit asks the gate for: `kind` ("exclusive" now) of the step `gate`, for `label`'s sake, after `attempt` failed tries."""

    def __init__(self, kind, gate, label, attempt, error):
        self.kind, self.gate, self.label, self.attempt, self.error = kind, gate, label, attempt, error

    def __str__(self):
        return f"{self.label} needs {self.kind} use of {self.gate!r} after {self.attempt} tries (MySQL error {self.error})"


class MySQLGate:
    """
    The gate of one contended step, kept in the database so every worker process sees it. The exclusive grant is the named lock
    `planetgen.exclusive.<database>.<gate>`, held by a connection of its own for the grant's length.

    Args:
        conn (Connection): Any connection of the database (its config opens the gate's own).
        gate (str): Names the step.
    """

    def __init__(self, conn, gate):
        self._config = conn._config
        database = self._config.database if self._config is not None else ""
        self.name = f"planetgen.exclusive.{database}.{gate}"[:64]
        self._holder = None

    def is_closed(self):
        """Whether somebody holds the exclusive grant."""
        probe = store.get_connection(self._config)
        try:
            return probe.execute("SELECT IS_FREE_LOCK(?) AS free", (self.name,)).fetchone()["free"] != 1
        finally:
            probe.close()

    def grant(self, need):
        """Takes the exclusive grant, waiting up to `GRANT_TIMEOUT_S` for another unit's. Returns whether it was granted."""
        from planetgen.queue import unitbeat
        holder = store.get_connection(self._config)
        try:
            with unitbeat.waiting("mysql"):         # OPS.136: a unit waiting for another unit's grant is waiting, not stuck
                ok = holder.execute("SELECT GET_LOCK(?, ?) AS ok", (self.name, GRANT_TIMEOUT_S)).fetchone()["ok"] == 1
        except Exception:
            holder.close()
            raise
        if not ok:
            holder.close()
            return False
        self._holder = holder
        return True

    def release(self):
        holder, self._holder = self._holder, None
        if holder is None:
            return
        try:
            holder.execute("SELECT RELEASE_LOCK(?)", (self.name,))
        except pymysql.err.MySQLError:      # a lost session has already dropped it
            pass
        finally:
            holder.close()


def _wait_while_closed(gate, label):
    """Holds off while another unit has the exclusive grant, with growing waits; says so once and again when it is open."""
    waited, wait, said = 0.0, 0.0, False
    while gate.is_closed():
        if not said:
            log.normal(f"{label} holds off: another unit has been granted exclusive use of this step.")
            said = True
        wait = min(1.0, max(0.05, wait * 1.5))
        time.sleep(wait)
        waited += wait
    if said:
        log.normal(f"{label} goes on after {waited:.1f} s.")


def begin_read_committed(conn):
    """
    Starts a transaction at READ COMMITTED (DB.113, DB.114): the Object ID step of the scatter tables (`scatter_uids.assign`) reads every row that
    has no ID yet, and at the default level that read takes locks that meet the other writers'. Call it first in an `attempt`.
    """
    conn.commit()           # a connection that just read has a transaction open, and MySQL refuses the next line inside one
    conn.execute("SET TRANSACTION ISOLATION LEVEL READ COMMITTED")


def run_with_retries(conn, attempt, gate, label, attempts=ATTEMPTS, gate_object=None):
    """
    Runs `attempt(conn)` (one whole transaction, committed at its end) and returns what it returns, trying again after a deadlock or lock wait
    timeout (module note): `attempts` plain tries with dynamic waits, then one more with the gate's exclusive grant.

    Args:
        conn (Connection): The unit's connection.
        attempt (callable): `attempt(conn)`; it starts from a clean transaction (this function rolls back after a failure).
        gate (str): Names the contended step; units of one gate hold off for each other's grants.
        label (str): What the log calls the unit.
        attempts (int): Plain tries before the need.
        gate_object (optional): A gate to use instead of the database's (`MySQLGate`): `is_closed()`, `grant(need)` and `release()`.

    Raises:
        pymysql.err.OperationalError: Any other error at once; a retryable one that persists through the granted try, or a grant that does
            not come in `GRANT_TIMEOUT_S`.
    """
    gate_object = gate_object if gate_object is not None else MySQLGate(conn, gate)
    wait = 0.0
    for number in range(1, attempts + 1):
        _wait_while_closed(gate_object, label)
        try:
            result = attempt(conn)
            _note(False)
            return result
        except pymysql.err.OperationalError as exc:
            conn.rollback()
            if not exc.args or exc.args[0] not in RETRYABLE_ERRORS:
                if exc.args and exc.args[0] in store.LOST_CONNECTION_ERRORS and getattr(conn, "_config", None) is not None:
                    connection_loss.explain_lost(conn._config, exc)       # OPS.171: why the server dropped it, on the way out
                raise
            _note(True)
            error = exc.args[0]
            if number == attempts:
                break
            wait = next_wait(number, wait)
            log.get("generate").notice(f"{label}: MySQL error {error} on try {number} of {attempts}; trying again in {wait:.2f} s "
                                       f"(this process's trouble level {_pressure:.2f}).")      # OPS.167: an expected, notable event
            time.sleep(wait)
    need = Need("exclusive", gate, label, attempts, error)
    log.normal(f"{need}: asking for it.")
    if not gate_object.grant(need):
        raise pymysql.err.OperationalError(1205, f"{label}: the exclusive grant for {gate!r} did not come in {GRANT_TIMEOUT_S} s")
    log.normal(f"{label}: granted exclusive use of {gate!r}; trying once more alone.")
    try:
        result = attempt(conn)
        _note(False)
        return result
    except pymysql.err.OperationalError:
        conn.rollback()
        _note(True)
        log.warning(f"{label}: the granted try failed too.")
        raise
    finally:
        gate_object.release()
        log.normal(f"{label}: released exclusive use of {gate!r}.")
