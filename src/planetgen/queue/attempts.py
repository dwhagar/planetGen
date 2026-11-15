# planetgen/queue/attempts.py

"""
One attempt of a unit, as the ledger sees it (OPS.217, L3 of `docs/design/operations/job-engine.md`): the fenced done mark.

A unit has a row in the job ledger (`engine_jobs`, made by the cut of its stage, `planrows.record_cut`). When a worker takes the unit it raises that row's attempt number and writes the number
to the unit's row in the galaxy's own database (`engine_items`, `ledger.begin_attempt`). Every commit the unit then makes to the galaxy goes through `Attempt.commit`, which updates the row
`WHERE attempt = mine` in the same transaction as the data: no row touched means another attempt took the unit over (a retry started while this one was still alive, a worker thought dead), so this
one rolls back and stops (`ledger.Replaced`), and its data never lands beside the other attempt's. The last commit sets the row's done time, in the same transaction as the last data, and only then
does the unit write its result and join the done set in Redis.

* **A unit that does not go through the engine's `ctx`** (every command today: none gives `pieces` yet, so each runs as its own `unit` function) is fenced once, at the end: its own commits are not
  fenced, so a replaced attempt of one of them can still write its earlier commits twice, as before; its done mark is fenced, so it cannot count twice. A family that converts to `pieces` and `do`
  fences every commit it makes through `ctx.commit()`.
* **The done set stays the barrier** of a stage until OPS.218 (open, conclude and finish become jobs and the supervisor reads the row): this item makes the unit's done mark true first, not Redis's.
* **A unit with no row** (a run cut before the plan was rows, a test's own stage) or no database in its arguments has nothing to fence and runs as before.
"""

import time
import traceback

from planetgen.db import store
from planetgen.queue import attemptlog, ledger, planrows, usage
from planetgen.util import log

Replaced = ledger.Replaced
"""exception: Another attempt took the unit over: the caller rolls back and stops."""


class Attempt:
    """What a worker holds for one attempt: the unit's row, the attempt's number, the control connection that writes its history and where the galaxy is."""

    def __init__(self, control, row, number, config):
        self.control, self.row, self.number, self.config = control, row, number, config
        self._cpu_at_start = usage.cpu_seconds()

    def commit(self, galaxy, cursor_key=None, done=False, results=()):
        """
        Fences the commit `galaxy` (a connection to the unit's database) is about to make, then makes it: the unit's row is updated `WHERE attempt = mine` (with the place reached and, for the last
        commit, the done time) in the transaction that holds the data. `results` are the `(piece_index, result)` pairs of the pieces finished since the last checkpoint, stored in the same
        transaction (OPS.219).

        Raises:
            Replaced: The row is another attempt's; the transaction is rolled back.
        """
        try:
            ledger.commit_item(galaxy, self.row, self.number, cursor_key=cursor_key, done=done, results=results)
        except Replaced:
            galaxy.rollback()
            raise
        galaxy.commit()

    def resume_point(self):
        """Where this attempt starts (OPS.219): `(cursor, results)`, the pieces the unit's earlier attempts finished and what they returned (`ledger.resume_point`); `(0, [])` for a unit with no database."""
        if self.config is None:
            return 0, []
        galaxy = store.get_connection(self.config)
        try:
            return ledger.resume_point(galaxy, self.row, self.number)
        finally:
            galaxy.close()

    def finish(self):
        """The done mark of a unit whose commits were not fenced one by one: opens the unit's database, sets the row done where the attempt is still mine, commits."""
        if self.config is None:
            return
        galaxy = store.get_connection(self.config)
        try:
            self.commit(galaxy, done=True)
        finally:
            galaxy.close()

    def previous(self):
        """The row of the attempt before this one (`ledger.attempt_row`), or `None` for a first attempt or a row that cannot be read."""
        if self.number < 2:
            return None
        try:
            return ledger.attempt_row(self.control, self.row, self.number - 1)
        except Exception as exc:  # noqa: BLE001 -- the header of a debug log is not worth the attempt
            log.swallowed("queue/attempts", "previous", exc)
            return None

    def end(self, outcome, error_kind=None, error_text=None, peak_mb=None, cpu_seconds=None, stack_text=None):
        """
        Writes how the attempt ended and the job's state to the ledger (a failure to write is logged: the unit's own outcome stands), and closes the control connection. The row keeps what
        the attempt cost: the processor seconds this process and its children used since it began (`cpu_seconds`, measured here unless given), the most memory it held (`peak_mb`) and, for
        a failure, the text of `stack_text` in a stack file (`attemptlog.keep_stack`) its row names.
        """
        try:
            stack_file = attemptlog.keep_stack(self.row, self.number, stack_text)
            if cpu_seconds is None:
                cpu_seconds = round(max(0.0, usage.cpu_seconds() - self._cpu_at_start), 3)
            ledger.end_attempt(self.control, self.row, self.number, outcome, peak_mb=peak_mb, cpu_seconds=cpu_seconds, error_kind=error_kind, error_text=error_text, stack_file=stack_file)
            ledger.set_job_state(self.control, self.row, {"done": "done", "failed": "failed"}.get(outcome, "waiting"))
            self.control.commit()
        except Exception as exc:  # noqa: BLE001 -- the ledger's history never ends a unit that did its work
            log.swallowed("queue/attempts", "end", exc)
        finally:
            self.close()

    def close(self):
        try:
            self.control.close()
        except Exception as exc:  # noqa: BLE001
            log.swallowed("queue/attempts", "close", exc)


def begin(conn, run_id, row, args, pid=None):
    """
    Takes the attempt: raises the unit's attempt number and writes it to the unit's row in the galaxy. Returns the `Attempt`, or `None` for a unit that has no row (`row` is `None`). A unit whose
    arguments name no database has a ledger row and no galaxy row: its attempts are counted and not fenced.
    """
    if row is None:
        return None
    control = planrows.control_connection(conn, run_id)
    try:
        named = (args or {}).get("database")
        config = store.MySQLConfig.from_args(named) if named else None
        galaxy = store.get_connection(config) if config is not None else None
        try:
            number = ledger.begin_attempt(control, row, galaxy, pid=pid)
        finally:
            if galaxy is not None:
                galaxy.close()
    except Exception:
        control.close()
        raise
    return Attempt(control, row, number, config)


def failure(exc):
    """`(error kind, error text, outcome)` of an exception that ended an attempt: the kind is its class, the text its traceback, the outcome `timed_out` for the time limit and `failed` otherwise."""
    kind = type(exc).__name__
    text = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
    return kind, text, "timed_out" if kind == "JobTimeoutException" else "failed"


WAIT_NOTE_SECONDS = 15.0
"""float: How often a worker that goes on waiting at the gate writes the same reason again (a changed reason is written at once)."""

_NOTED = {}
"""dict: Job row to `(reason, time)` of the last note written, so a wait of an hour writes a few rows, not one at every poll."""


def note_wait(conn, run_id, row, reason, now=None):
    """
    Writes why a unit is not starting (OPS.224): the reason the memory gate gave the worker that holds it, onto the job's ledger row (`engine_jobs.waited`), where the Queue page reads it. Called at
    every poll of a wait at the gate; written when the reason changes and then at most every `WAIT_NOTE_SECONDS`. Best effort: a failure to write is logged, never raised. `conn` is the run's Redis.
    """
    if row is None:
        return False
    now = time.time() if now is None else now
    before = _NOTED.get(row)
    if before is not None and before[0] == reason and now - before[1] < WAIT_NOTE_SECONDS:
        return False
    try:
        control = planrows.control_connection(conn, run_id)
    except Exception as exc:  # noqa: BLE001 -- a note on a wait is never a reason to stop the wait
        log.swallowed("queue/attempts", "note_wait", exc)
        return False
    try:
        taken = ledger.set_wait(control, row, reason)
        control.commit()
    except Exception as exc:  # noqa: BLE001
        log.swallowed("queue/attempts", "note_wait", exc)
        return False
    finally:
        control.close()
    _NOTED[row] = (reason, now)
    if len(_NOTED) > 1000:
        _NOTED.pop(next(iter(_NOTED)))
    return taken


def close_lost(conn, run_id, row, outcome, error_kind=None, error_text=None, stack_text=None):
    """
    Writes the end of an attempt that could not write it itself (OPS.223): a horse killed for memory (`killed_memory`) or by something else (`lost`), a unit the follower stopped as stuck
    (`stuck`, with the stack it had the horse write). Closes the job's attempts that are still `running` in the ledger and moves the job back to `waiting`, for the retry. Best effort:
    a failure to write is logged, never raised. `conn` is the run's Redis, which names the control database. Returns whether an attempt was closed.
    """
    if row is None:
        return False
    try:
        control = planrows.control_connection(conn, run_id)
    except Exception as exc:  # noqa: BLE001 -- the ledger's history never stops what saw the kill
        log.swallowed("queue/attempts", "close_lost", exc)
        return False
    try:
        newest = control.execute("SELECT attempt FROM engine_jobs WHERE id = ?", (row,)).fetchone()
        stack_file = attemptlog.keep_stack(row, newest["attempt"], stack_text) if newest else None
        numbers = ledger.end_running(control, row, outcome, error_kind=error_kind, error_text=error_text, stack_file=stack_file)
        if numbers:
            ledger.set_job_state(control, row, "waiting")
        control.commit()
        return bool(numbers)
    except Exception as exc:  # noqa: BLE001
        log.swallowed("queue/attempts", "close_lost", exc)
        return False
    finally:
        control.close()
