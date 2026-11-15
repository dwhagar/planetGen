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

from planetgen.db import store
from planetgen.queue import ledger, planrows
from planetgen.util import log

Replaced = ledger.Replaced
"""exception: Another attempt took the unit over: the caller rolls back and stops."""


class Attempt:
    """What a worker holds for one attempt: the unit's row, the attempt's number, the control connection that writes its history and where the galaxy is."""

    def __init__(self, control, row, number, config):
        self.control, self.row, self.number, self.config = control, row, number, config

    def commit(self, galaxy, cursor_key=None, done=False):
        """
        Fences the commit `galaxy` (a connection to the unit's database) is about to make, then makes it: the unit's row is updated `WHERE attempt = mine` (with the place reached and, for the last
        commit, the done time) in the transaction that holds the data.

        Raises:
            Replaced: The row is another attempt's; the transaction is rolled back.
        """
        try:
            ledger.commit_item(galaxy, self.row, self.number, cursor_key=cursor_key, done=done)
        except Replaced:
            galaxy.rollback()
            raise
        galaxy.commit()

    def finish(self):
        """The done mark of a unit whose commits were not fenced one by one: opens the unit's database, sets the row done where the attempt is still mine, commits."""
        if self.config is None:
            return
        galaxy = store.get_connection(self.config)
        try:
            self.commit(galaxy, done=True)
        finally:
            galaxy.close()

    def end(self, outcome, error_kind=None, error_text=None, peak_mb=None, cpu_seconds=None):
        """Writes how the attempt ended and the job's state to the ledger (a failure to write is logged: the unit's own outcome stands), and closes the control connection."""
        try:
            ledger.end_attempt(self.control, self.row, self.number, outcome, peak_mb=peak_mb, cpu_seconds=cpu_seconds, error_kind=error_kind, error_text=error_text)
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
