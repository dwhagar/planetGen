# planetgen/queue/engine.py

"""
The engine's loop over a unit's pieces (OPS.191; reports/job-engine-design-2026-11-13.md 3.4, reports/job-engine-items-stage-1-2-2026-11-13.md item P1).

A command that declares `pieces` and `do` (`commands.CommandSpec`) does not write its own loop. The engine asks `pieces(args, cursor)` for the list of pieces of the unit that come after
`cursor`, and calls `do(args, piece, ctx)` for each in order. What the engine does for every such command, so that no command writes it:

* **Counts.** Pieces done of pieces in all go to the unit's beat (`report(done, total)`) before the first piece and after every one, so a unit that works never goes quiet and the bar of the run
  shows real numbers. A command cannot forget it. A piece that is a dict may carry `weight`, how much of the unit's work it is (default 1; PERF.268): the beat counts weight done of weight in all, so
  a unit of pieces that cost unequal amounts (a dense sector against a sparse one) has a bar, a judge's expected rate and a profile that follow the work and not the number of pieces. The cursor,
  the fence and the done mark count pieces, never weights.
* **Stop and pause between two pieces.** Before a piece, at most every `CHECK_SECONDS`, the engine reads whether the run is still live and whether the queues are held. A run that is no longer
  live (canceled, failed, ended from outside) ends the unit with `UnitStopped`, which `units._do_unit` answers by leaving the unit undone; a pause holds the unit where it stands, shown in
  the beat as waiting on `pause`, until the queues are let go (a solo hold lets its own unit through, `pause.blocks`) or the run ends.
* **`ctx`**, the one object a piece is handed: `ctx.args` and `ctx.kind`, `ctx.cursor`, the unit's database connection (`ctx.db`, opened from `args["database"]` the first time it is asked
  for and closed by the engine), `ctx.commit()`, `ctx.sql(statement, label)` (`queue/longsql.py`) and the second-level count, `ctx.step(name, total)` and `ctx.advance(n)` (OPS.192: a bar under the unit's own on the pages). A converted command commits only through `ctx.commit()`: where the cursor is stored with the
  data (OPS.219) that call is the one place that writes it, and a direct `.commit()` would be a second path.

**Checkpoints (OPS.219).** A command that `resumes` (`commands.CommandSpec`) starts after the pieces an earlier attempt of the unit finished. `ctx.commit()` is the one place the cursor is written: with a unit that
has a ledger row it is the fence (OPS.217: `attempts.Attempt.commit`) and the same statement stores the number of pieces finished and, in the same transaction, what the newly finished ones returned
(`ledger.commit_item`), so the data, the cursor and the results land together or not at all. After every piece the engine checkpoints (`Context.checkpoint`), so the cursor is never behind data a piece left for the engine to commit. A unit run again reads its cursor and results
(`attempts.Attempt.resume_point`), asks `pieces(args, cursor)` for the rest and gathers the stored results with the new ones, so its result is the one an undisturbed unit gives. A piece's result goes through
JSON before the unit keeps it, so a resumed unit and an undisturbed one hand `gather` the same kinds of value. A command that does not `resume` always starts at 0 and stores no cursor. A unit with no
ledger row (a command run alone, a test's own stage) starts at 0 and commits plainly.

**A command without `pieces` runs as it always has** (`commands.unit_function`) until the item of its family converts it. Nothing here changes what a command computes.
"""

import json
import math
import time

from planetgen.queue import commands, longsql, pause, runrecord, unitbeat
from planetgen.util import log, ticker

CHECK_SECONDS = 1.0
"""float: At most how often the engine reads the run's state and the pause flag between two pieces. Defined by the program, not a physical value: a unit of thousands of pieces a second
costs Redis one read of each a second, and a cancel or a pause is honoured within a second of the piece in hand ending."""

PAUSE_POLL_SECONDS = 2.0
"""float: How often a unit held by a pause looks at the flag and the run again. Defined by the program, not a physical value (the same as `pause.WAIT_SECONDS`, which a worker waits by)."""

PAUSE_TAG = "pause"
"""str: What the beat says a unit held between two pieces by a pause waits on (`unitbeat.waiting`)."""


def weight_of(piece):
    """
    How much of its unit's work a piece is (PERF.268): the `weight` of a piece that is a dict and names a positive finite number, else 1.

    Args:
        piece: One piece, as the command's `pieces` made it.

    Returns:
        float: The weight, over 0.
    """
    weight = piece.get("weight") if isinstance(piece, dict) else None
    if isinstance(weight, bool) or not isinstance(weight, (int, float)) or not math.isfinite(weight) or weight <= 0:
        return 1.0
    return float(weight)


class UnitStopped(Exception):
    """The run the unit belongs to is no longer live: the engine ends the unit between two pieces, with the pieces done so far kept."""


class Context:
    """
    What a piece is handed (`do(args, piece, ctx)`): the unit's arguments, where it stands, and the database and the long statements, so no command opens its own or writes its own commit.
    """

    def __init__(self, kind, args, report, *, redis=None, run_id=None, stage=None, unit=None, job=None, cursor=0, clock=time.monotonic, attempt=None, halt=None):
        self.kind, self.args, self.report, self.cursor = kind, args, report, cursor
        self.redis, self.run_id, self.stage, self.unit, self.job = redis, run_id, stage, unit, job
        self.clock = clock
        self.attempt = attempt
        self.done = cursor                  # the pieces finished, those an earlier attempt did included (OPS.219)
        self.results = []                   # what every finished piece returned, in order; `None` for a piece that returned nothing
        self._stored = cursor               # how many of them the ledger holds
        self._resumes = False
        self.halt = halt
        self._closers = []
        self.state = {}                 # what the pieces of this unit hand to each other: they run in one process, in order, so a piece may leave what the next reads
        self.began = time.monotonic()
        self._db = None
        self._checked_at = None

    # --- the database ----------------------------------------------------------------------------------------

    @property
    def config(self):
        """The `MySQLConfig` the unit's arguments name (`args["database"]`, as `wrapper.database_args` writes it; the password is this process's settings')."""
        from planetgen.db import store
        return store.MySQLConfig.from_args(self.args["database"])

    @property
    def db(self):
        """The unit's connection, opened on first use and closed when the unit ends."""
        if self._db is None:
            from planetgen.db import store
            self._db = store.get_connection(self.config)
        return self._db

    def kept(self, name, make):
        """
        What the unit works out once and every piece reads (the galaxy's shape, a table of components): `make()` the first time `name` is asked for, the same value after, until the unit ends
        (`state` is the unit's own, so a value never outlives it and a later unit reads the database afresh).
        """
        if name not in self.state:
            self.state[name] = make()
        return self.state[name]

    def commit(self):
        """
        Commits what the piece wrote, fenced (OPS.217): given the unit's attempt (`attempts.Attempt`) the unit's row in the galaxy is updated `WHERE attempt = mine` in the same transaction, and a
        unit another attempt has taken over rolls back and raises `attempts.Replaced`. With the cursor in the ledger (OPS.219) the same statement stores it. Without an attempt (a unit with no ledger
        row, a command run alone) it is a plain commit. A command that has opened no connection has nothing to commit.
        """
        if self._db is None:
            return
        if self.attempt is not None and self.attempt.config is not None:
            if self._resumes:
                pending = [(index, result) for index, result in enumerate(self.results[self._stored:], self._stored) if result is not None]
                self.attempt.commit(self._db, cursor_key=str(self.done), results=pending)
                self._stored = self.done
            else:
                self.attempt.commit(self._db)
        else:
            self._db.commit()

    def commit_piece(self, result=None):
        """
        Commits the piece in hand as finished (OPS.222), for a piece whose own commit must be the one that counts it: its data, the cursor past it and `result` land in one fenced transaction, so a
        unit stopped after it starts at the piece that follows and one stopped before it does this piece again. A piece that would write twice if it ran twice (rows added to a total, rows inserted
        with no key to refuse them) commits so. `result` is what the piece returns (the engine stores it here and not again). Without a database, an attempt or `resumes` it is a plain `commit()`.
        """
        if self._db is None:
            return
        if self._resumes and self.attempt is not None and self.attempt.config is not None:
            number = len(self.results)          # the place of the piece in hand: the pieces before it, an earlier attempt's included, are `results`
            pending = [(index, stored) for index, stored in enumerate(self.results[self._stored:], self._stored) if stored is not None]
            stored = json.loads(json.dumps(result))
            if stored is not None:
                pending.append((number, stored))
            self.attempt.commit(self._db, cursor_key=str(number + 1), results=pending)
            self._stored = number + 1
        else:
            self.commit()

    def checkpoint(self):
        """
        After every piece (OPS.219): stores the number of pieces finished and what the new ones returned, with a commit of anything the piece left uncommitted, fenced like any `commit()`. A unit that does
        not `resumes`, or has no ledger row, stores nothing. A piece that wrote through a connection of its own has committed that before this runs, so a unit stopped between the two does the one piece
        again: a resumable command's pieces are safe to run again from their start (`commands.CommandSpec`).
        """
        if not self._resumes or self.attempt is None or self.attempt.config is None:
            return
        self.db         # the unit's connection, opened here for a piece that wrote through its own
        self.commit()

    def sql(self, statement, label, params=()):
        """One long statement with its label and start time in the beat (`longsql.run`); returns the cursor."""
        return longsql.run(self.db, statement, label, params, config=self.config)

    def step(self, name, total=None):
        """The piece begins a second-level count (OPS.192): `name`, `total` of what it counts (`None` where that is not known). The Queue page draws it as a bar under the unit's own, with its numbers and rate."""
        unitbeat.step(name, total)

    def advance(self, amount=1):
        """The second-level count has done `amount` more."""
        unitbeat.advance(amount)

    def moving(self):
        """The piece made progress its count does not show: the quiet clock starts again (`unitbeat.moving`)."""
        unitbeat.moving()

    def phase(self, text, keep=None, loud=True):
        """
        The piece says what it is doing (the standard's rule 6, OPS.170): a log line with the time since the unit began, the unit's phase for the heartbeat and, where `keep` is a list,
        the line kept in it for the unit's result. `loud` false makes the line a debug record, for one unit of many.
        """
        line = ticker.stamp(text, self.began)
        (log.normal if loud else log.debug)(line)
        unitbeat.phase(text[:120])
        if keep is not None:
            keep.append(line)

    def close_db(self):
        """Closes the unit's connection now, where the rest of the unit needs none held (`ctx.db` opens another on the next use)."""
        db, self._db = self._db, None
        if db is not None:
            db.close()

    def on_close(self, closer):
        """Registers `closer()` to be called, last in first out, when the unit ends, however it ends: a resource the pieces share (a pool of processes) that no single piece can close."""
        self._closers.append(closer)

    def close(self):
        """Closes what the unit opened, called by the engine, once, when the unit ends: what `on_close` registered, then the unit's connection if it opened one."""
        while self._closers:
            self._closers.pop()()
        self.close_db()

    # --- between two pieces ----------------------------------------------------------------------------------

    def between(self):
        """
        Called before every piece: first what `halt` says (a local read, every piece), then, when `CHECK_SECONDS` have passed since the last look, reads whether the run is still live and holds the unit while the queues are paused.

        Raises:
            UnitStopped: The run is no longer live.
        """
        if self.halt is not None:
            reason = self.halt()
            if reason:
                raise UnitStopped(reason)
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


def run(kind, args, report, *, redis=None, run_id=None, stage=None, unit=None, job=None, clock=time.monotonic, attempt=None, halt=None):
    """
    Runs one unit of command `kind` and returns its result: through the engine's loop for a command with `pieces` and `do`, as the command's own `unit` function for one without.

    Args:
        report: `report(done, estimate)` of the unit's beat (`unitbeat.UnitBeat`), or any callable taking the same two numbers.
        attempt (attempts.Attempt, optional): The unit's attempt in the ledger; with it every `ctx.commit()` is fenced.
        halt (callable, optional): Asked before every piece (a local read, no throttle): returns why the unit is to stop there, or `None`. The job process the template forked (OPS.226) passes what the supervisor asks and whether it is still there.
        redis, run_id, stage, unit, job: Where the unit runs; with them the engine honours a cancel and a pause between two pieces. Without them (`commands.run_alone`) it does neither.

    Returns:
        The unit's JSON-able result: what the command's `gather(args, results, ctx)` makes of the pieces' results (and `ctx.state`), or, with none, the list of them in order.

    Raises:
        UnitStopped: The run ended between two pieces.
    """
    spec = commands.spec(kind)
    if spec.pieces is None:
        return commands.unit_function(kind)(args, report)
    pieces_of, do, gather = commands.pieces_function(kind), commands.do_function(kind), commands.gather_function(kind)
    resumes = spec.resumes and attempt is not None
    cursor, before = attempt.resume_point() if resumes else (0, [])
    ctx = Context(kind, args, report, redis=redis, run_id=run_id, stage=stage, unit=unit, job=job, cursor=cursor, clock=clock, attempt=attempt, halt=halt)
    ctx._resumes = resumes
    ctx.results = list(before)
    try:
        pieces = list(pieces_of(args, cursor))
        weights = [weight_of(piece) for piece in pieces]
        weight_done = cursor * (sum(weights) / len(weights)) if weights else float(cursor)      # the cursor counts pieces: the ones before it weigh what the rest weigh on average (PERF.268)
        total = weight_done + sum(weights)
        report(weight_done, total)
        for done, (piece, weight) in enumerate(zip(pieces, weights), cursor + 1):
            ctx.between()
            result = do(args, piece, ctx)
            ctx.results.append(json.loads(json.dumps(result)) if spec.resumes else result)
            ctx.done = done
            ctx.checkpoint()
            weight_done += weight
            report(weight_done, total)
        return ctx.results if gather is None else gather(args, ctx.results, ctx)
    finally:
        ctx.close()
