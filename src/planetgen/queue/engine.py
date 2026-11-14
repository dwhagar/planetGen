# planetgen/queue/engine.py

"""
The engine's loop over a unit's pieces (OPS.191; reports/job-engine-design-2026-11-13.md 3.4, reports/job-engine-items-stage-1-2-2026-11-13.md item P1).

A command that declares `pieces` and `do` (`commands.CommandSpec`) does not write its own loop. The engine asks `pieces(args, cursor)` for the list of pieces of the unit that come after
`cursor`, and calls `do(args, piece, ctx)` for each in order. What the engine does for every such command, so that no command writes it:

* **Counts.** Pieces done of pieces in all go to the unit's beat (`report(done, total)`) before the first piece and after every one, so a unit that works never goes quiet and the bar of the run
  shows real numbers. A command cannot forget it.
* **Stop and pause between two pieces.** Before a piece, at most every `CHECK_SECONDS`, the engine reads whether the run is still live and whether the queues are held. A run that is no longer
  live (canceled, failed, ended from outside) ends the unit with `UnitStopped`, which `units._do_unit` answers by leaving the unit undone; a pause holds the unit where it stands, shown in
  the beat as waiting on `pause`, until the queues are let go (a solo hold lets its own unit through, `pause.blocks`) or the run ends.
* **`ctx`**, the one object a piece is handed: `ctx.args` and `ctx.kind`, `ctx.cursor`, the unit's database connection (`ctx.db`, opened from `args["database"]` the first time it is asked
  for and closed by the engine), `ctx.commit()` and `ctx.sql(statement, label)` (`queue/longsql.py`). A converted command commits only through `ctx.commit()`: where the cursor is stored with the
  data (OPS.219) that call is the one place that writes it, and a direct `.commit()` would be a second path.

**The cursor is in the contract and not stored yet.** Every unit starts at cursor 0. Storing it, so a retry or a resume starts after the last committed piece, waits for the ledger (OPS.219);
until then `ctx.commit()` is a plain commit.

**A command without `pieces` runs as it always has** (`commands.unit_function`) until the item of its family converts it. Nothing here changes what a command computes.
"""

import time

from planetgen.queue import commands, longsql, pause, runrecord, unitbeat
from planetgen.util import log

CHECK_SECONDS = 1.0
"""float: At most how often the engine reads the run's state and the pause flag between two pieces. Defined by the program, not a physical value: a unit of thousands of pieces a second
costs Redis one read of each a second, and a cancel or a pause is honoured within a second of the piece in hand ending."""

PAUSE_POLL_SECONDS = 2.0
"""float: How often a unit held by a pause looks at the flag and the run again. Defined by the program, not a physical value (the same as `pause.WAIT_SECONDS`, which a worker waits by)."""

PAUSE_TAG = "pause"
"""str: What the beat says a unit held between two pieces by a pause waits on (`unitbeat.waiting`)."""


class UnitStopped(Exception):
    """The run the unit belongs to is no longer live: the engine ends the unit between two pieces, with the pieces done so far kept."""


class Context:
    """
    What a piece is handed (`do(args, piece, ctx)`): the unit's arguments, where it stands, and the database and the long statements, so no command opens its own or writes its own commit.
    """

    def __init__(self, kind, args, report, *, redis=None, run_id=None, stage=None, unit=None, job=None, cursor=0, clock=time.monotonic):
        self.kind, self.args, self.report, self.cursor = kind, args, report, cursor
        self.redis, self.run_id, self.stage, self.unit, self.job = redis, run_id, stage, unit, job
        self.clock = clock
        self._db = None
        self._checked_at = None

    # --- the database ----------------------------------------------------------------------------------------

    @property
    def config(self):
        """The `MySQLConfig` the unit's arguments name (`args["database"]`, as `wrapper.database_args` writes it; the password is this process's settings')."""
        from planetgen.db import store
        named = self.args["database"]
        return store.MySQLConfig(host=named["host"], port=named["port"], user=named["user"], database=named["database"])

    @property
    def db(self):
        """The unit's connection, opened on first use and closed when the unit ends."""
        if self._db is None:
            from planetgen.db import store
            self._db = store.get_connection(self.config)
        return self._db

    def commit(self):
        """
        Commits what the piece wrote. In stage 1 of the engine this is a plain commit on `ctx.db`; with the ledger (OPS.219) it commits the cursor with the data. A command that has opened no
        connection has nothing to commit.
        """
        if self._db is not None:
            self._db.commit()

    def sql(self, statement, label, params=()):
        """One long statement with its label and start time in the beat (`longsql.run`); returns the cursor."""
        return longsql.run(self.db, statement, label, params, config=self.config)

    def moving(self):
        """The piece made progress its count does not show: the quiet clock starts again (`unitbeat.moving`)."""
        unitbeat.moving()

    def close(self):
        """Closes the unit's connection if it opened one; called by the engine, once, when the unit ends."""
        db, self._db = self._db, None
        if db is not None:
            db.close()

    # --- between two pieces ----------------------------------------------------------------------------------

    def between(self):
        """
        Called before every piece: when `CHECK_SECONDS` have passed since the last look, reads whether the run is still live and holds the unit while the queues are paused.

        Raises:
            UnitStopped: The run is no longer live.
        """
        if self.redis is None or self.run_id is None:
            return
        now = self.clock()
        if self._checked_at is not None and now - self._checked_at < CHECK_SECONDS:
            return
        self._check_live()
        if pause.blocks(self.redis, self.job):
            self._hold_while_paused()
        self._checked_at = self.clock()

    def _check_live(self):
        state = runrecord._s(self.redis.hget(runrecord.run_key(self.run_id), "state"))
        if state != "running":
            raise UnitStopped(f"the run is {state or 'gone'}")

    def _hold_while_paused(self):
        """Holds the unit between two pieces while the queues are paused; a run that ends meanwhile ends the hold (`UnitStopped`)."""
        log.debug(f"Unit {self.unit} of stage {self.stage} ({self.kind}) waits between two pieces: the queues are paused.", area="queue")
        with unitbeat.waiting(PAUSE_TAG):
            while pause.blocks(self.redis, self.job):
                self._check_live()
                time.sleep(PAUSE_POLL_SECONDS)


def run(kind, args, report, *, redis=None, run_id=None, stage=None, unit=None, job=None, clock=time.monotonic):
    """
    Runs one unit of command `kind` and returns its result: through the engine's loop for a command with `pieces` and `do`, as the command's own `unit` function for one without.

    Args:
        report: `report(done, estimate)` of the unit's beat (`unitbeat.UnitBeat`), or any callable taking the same two numbers.
        redis, run_id, stage, unit, job: Where the unit runs; with them the engine honours a cancel and a pause between two pieces. Without them (`commands.run_alone`) it does neither.

    Returns:
        The unit's JSON-able result: what the command's `gather` makes of the pieces' results, or, with none, the list of them in order.

    Raises:
        UnitStopped: The run ended between two pieces.
    """
    spec = commands.spec(kind)
    if spec.pieces is None:
        return commands.unit_function(kind)(args, report)
    pieces_of, do, gather = commands.pieces_function(kind), commands.do_function(kind), commands.gather_function(kind)
    cursor = 0              # OPS.219 stores it; until then a unit always starts at the beginning
    ctx = Context(kind, args, report, redis=redis, run_id=run_id, stage=stage, unit=unit, job=job, cursor=cursor, clock=clock)
    try:
        pieces = list(pieces_of(args, cursor))
        total = cursor + len(pieces)
        report(cursor, total)
        results = []
        for done, piece in enumerate(pieces, cursor + 1):
            ctx.between()
            results.append(do(args, piece, ctx))
            report(done, total)
    finally:
        ctx.close()
    return results if gather is None else gather(args, results)
