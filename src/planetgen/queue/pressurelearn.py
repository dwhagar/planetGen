# planetgen/queue/pressurelearn.py

"""
The edges of the memory controller, learned from recorded runs (OPS.249, H4 of the job engine; answers OPS.151's question: tune the edges from a real run's record).

`pressure.Controller` (OPS.231) holds the bulk and upkeep lanes when the host reads *pressed* and asks the youngest bulk job to stop when it reads *starved*. The edges that decide the
reading (`pressure.HOLD_SOME`, `YIELD_SOME`, `YIELD_FULL`, and without PSI `PRESSURE_SHARE` and `STARVED_SHARE`) were typed. This module records what the host read, a minute at a time, and
works out from the record how far to move them.

**Record.** Once a minute the controller's look is summed up in one row of `pressure_minutes` (control schema v25): the highest `some` and `full` of the minute, the lowest available memory, the
reserve, the highest level, the seconds the controller held admission, the number of jobs it asked to stop and whether the host read by PSI. A row is the whole of a minute whatever the number of looks. The record never
blocks the supervisor: a failed write drops the minute, and the writer leaves the database alone for `BACKOFF_SECONDS`. `PLANETGEN_PRESSURE_LEARN=0` turns the record and the learning off.

**One number.** The five edges move together by one *sensitivity*, `k`, between `MIN_SENSITIVITY` and `MAX_SENSITIVITY` (so the typed edges are the middle and the clamp is an eighth to a
double of them): the `some` and `full` edges are divided by `k`, and the shares of the reserve are multiplied by it. The ratios between the edges stay the typed stand-ins; only where the
whole ladder sits is learned. `k` is worked out from the record every time, so it has no state of its own to drift or to clear: forget the rows and the edges are the typed ones.

**What the record says.**

* *A kill the hold did not see coming* (a `killed_memory` attempt of the ledger with no pressed reading in the `LEAD_SECONDS` before it): the edge was too high. The strength of the reading before
  each such kill (the highest `some` over the edge, or the reserve's share the lowest available memory fell to) is taken at its lower quartile, and `k` is set so that this reading would have
  been pressed with `MARGIN` to spare. It needs `MIN_KILLS` of them; fewer says nothing.
* *Holds that were not needed.* A hold is one run of minutes with seconds held. It was *needed* when the host reached starved in it, or a job was asked to stop, or a job was killed in it or in the
  `LEAD_SECONDS` after. With `MIN_EPISODES` holds, no kill at all in the window and at least `IDLE_SHARE` of the holds not needed, the edge moves up to where the upper quartile of those holds would
  not have been pressed, never more than `MAX_RAISE`.

Both rules look at the last `WINDOW_DAYS` days. The numbers here are design stand-ins; they are what a record is read with, not what the host is held to.
"""

import os
import time
from typing import NamedTuple

from planetgen.queue import pressure
from planetgen.util import log

ENV_VAR = "PLANETGEN_PRESSURE_LEARN"
"""str: Set to `0`, `off`, `no` or `false` and nothing is recorded and the typed edges stand."""

WINDOW_DAYS = 14
"""int: How far back the record is read. Stand-in."""

LEAD_SECONDS = 180.0
"""float: How long before a kill (and after a hold) the record is looked at. Stand-in."""

MIN_KILLS = 3
"""int: The kills the hold did not see coming that the edge needs before it moves down. Stand-in."""

MIN_EPISODES = 8
"""int: The holds the edge needs before it moves up. Stand-in."""

IDLE_SHARE = 0.8
"""float: The share of holds that were not needed before the edge moves up. Stand-in."""

MARGIN = 0.8
"""float: The edge is put this share of the way to the reading that came before a kill, so such a reading would have been pressed. Stand-in."""

MAX_RAISE = 1.5
"""float: The most the edge moves up in one reading of the record, as a multiple of where it is. Stand-in."""

MIN_SENSITIVITY, MAX_SENSITIVITY = 0.5, 2.0
"""float: The clamp of `k`; at the first the edges are twice the typed ones, at the second half. Stand-in."""

BACKOFF_SECONDS = 60.0
"""float: After a failed write the recorder leaves the database alone for this long."""

REFRESH_SECONDS = 600.0
"""float: How often the controller asks for its edges again."""

_MINUTE = 60.0


class Edges(NamedTuple):
    """The edges `pressure.read` judges a reading by, with the sensitivity they were worked out at (1 is the typed ones) and a note of why."""

    hold_some: float
    yield_some: float
    yield_full: float
    pressure_share: float
    starved_share: float
    sensitivity: float = 1.0
    note: str = "the typed edges"

    def words(self):
        return self.note


def typed():
    """The typed edges of `pressure`."""
    return scaled(1.0, "the typed edges")


def scaled(sensitivity, note):
    """The typed edges at `sensitivity`: the PSI edges divided by it, the shares of the reserve multiplied by it (a share cannot pass 1)."""
    return Edges(pressure.HOLD_SOME / sensitivity, pressure.YIELD_SOME / sensitivity, pressure.YIELD_FULL / sensitivity,
                 min(1.0, pressure.PRESSURE_SHARE * sensitivity), min(1.0, pressure.STARVED_SHARE * sensitivity), sensitivity, note)


def enabled():
    return os.environ.get(ENV_VAR, "1").strip().lower() not in ("0", "off", "no", "false")


class Minute(NamedTuple):
    """One recorded minute: `start` (seconds since the epoch, the minute's own), the highest `some` and `full`, the lowest available MB, the reserve MB, the highest level, the seconds held, the jobs asked to stop, and whether the host read by PSI."""

    start: float
    some: float
    full: float
    available_mb: float
    reserve_mb: float
    level: int
    held_seconds: float
    yields: int
    psi: bool = True


def strength(minute):
    """
    How far over the typed pressed edge a minute read: 1 is just pressed, under 1 calm. With PSI it is the highest `some` over `pressure.HOLD_SOME`; without it the reserve's pressed share over the share of
    it the host had left (the lowest available memory over the reserve, against `pressure.PRESSURE_SHARE`).
    """
    if minute.psi:
        return minute.some / pressure.HOLD_SOME
    if minute.available_mb <= 0.0:
        return float("inf")
    return pressure.PRESSURE_SHARE * minute.reserve_mb / minute.available_mb


def _quantile(values, share):
    ordered = sorted(values)
    return ordered[min(len(ordered) - 1, max(0, round(share * (len(ordered) - 1))))]


def _episodes(minutes):
    """The holds of a record: lists of consecutive minutes (a minute apart) in which the controller held admission."""
    found, current = [], []
    for minute in sorted(minutes, key=lambda item: item.start):
        if minute.held_seconds > 0.0:
            if current and minute.start - current[-1].start > _MINUTE + 1.0:
                found.append(current)
                current = []
            current.append(minute)
        elif current:
            found.append(current)
            current = []
    if current:
        found.append(current)
    return found


def learn(minutes, kills):
    """
    The `Edges` the record gives. `minutes` are `Minute`s over the window, `kills` the times (seconds since the epoch) of the memory kills in it. With too little record, or a record that
    says nothing, the typed edges.
    """
    minutes = sorted(minutes, key=lambda item: item.start)
    if not minutes:
        return typed()
    by_start = {minute.start: minute for minute in minutes}
    uncaught = []
    for killed in kills:
        lead = [minute for start, minute in by_start.items() if killed - LEAD_SECONDS - _MINUTE <= start <= killed]
        if not lead:
            continue
        top = max(strength(minute) for minute in lead)
        if 0.0 < top < 1.0:
            uncaught.append(top)
    if len(uncaught) >= MIN_KILLS:
        reading = _quantile(uncaught, 0.25)
        sensitivity = min(MAX_SENSITIVITY, max(1.0, 1.0 / (MARGIN * reading)))
        return scaled(sensitivity, f"{len(uncaught)} memory kills came with no pressed reading before them (the lower quartile read {reading:.0%} of the edge): the edges are {sensitivity:.2f} times as sensitive")
    if kills:
        return typed()
    episodes = _episodes(minutes)
    if len(episodes) < MIN_EPISODES:
        return typed()
    idle = [episode for episode in episodes if not any(minute.level >= pressure.STARVED or minute.yields for minute in episode)]
    if len(idle) / len(episodes) < IDLE_SHARE:
        return typed()
    peak = _quantile([max(strength(minute) for minute in episode) for episode in idle], 0.75)
    raised = min(MAX_RAISE, max(1.0, peak))
    sensitivity = max(MIN_SENSITIVITY, 1.0 / raised)
    if sensitivity >= 1.0:
        return typed()
    return scaled(sensitivity, f"{len(idle)} of {len(episodes)} holds ended with no starved reading, no job asked to stop and no kill: the edges are {raised:.2f} times as high")


class Recorder:
    """
    Sums the controller's looks up a minute at a time and writes the rows (`pressure_minutes`), and hands back the edges the record gives (`edges`, from `learn`). `connect` makes a control connection. Never raises.
    """

    def __init__(self, connect, clock=time.time, monotonic=time.monotonic):
        self.connect, self.clock, self.monotonic = connect, clock, monotonic
        self.current = None         # the open minute: [start, some, full, available_mb, reserve_mb, level, held_seconds, yields, last look, psi]
        self.retry_after = 0.0
        self.cached = typed()
        self.cached_at = None

    def note(self, reading, held, asked=False):
        """One look at `reading`: `held` says whether the controller held admission, `asked` whether it asked a job to stop on this look."""
        if not enabled():
            return
        now = self.clock()
        start = now - now % _MINUTE
        if self.current is not None and self.current[0] != start:
            self._write(self.current)
            self.current = None
        if self.current is None:
            self.current = [start, 0.0, 0.0, reading.available_mb, reading.reserve_mb, 0, 0.0, 0, now, reading.source == "psi"]
        row = self.current
        row[1], row[2] = max(row[1], reading.some), max(row[2], reading.full)
        row[3], row[4], row[5] = min(row[3], reading.available_mb), reading.reserve_mb, max(row[5], reading.level)
        if held:
            row[6] = min(_MINUTE, row[6] + max(0.0, now - row[8]))
        row[7] += 1 if asked else 0
        row[8] = now

    def _write(self, row):
        if self.monotonic() < self.retry_after:
            return
        try:
            conn = self.connect()
            try:
                conn.execute("INSERT INTO pressure_minutes (minute_start, some_max, full_max, available_min_mb, reserve_mb, level_max, held_seconds, yields, psi)"
                             " VALUES (FROM_UNIXTIME(?), ?, ?, ?, ?, ?, ?, ?, ?) ON DUPLICATE KEY UPDATE some_max = GREATEST(some_max, VALUES(some_max)),"
                             " full_max = GREATEST(full_max, VALUES(full_max)), available_min_mb = LEAST(available_min_mb, VALUES(available_min_mb)),"
                             " level_max = GREATEST(level_max, VALUES(level_max)), held_seconds = LEAST(60, held_seconds + VALUES(held_seconds)), yields = yields + VALUES(yields)",
                             (row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], 1 if row[9] else 0))
                conn.commit()
                if int(row[0] // _MINUTE) % 1440 == 0:
                    conn.execute("DELETE FROM pressure_minutes WHERE minute_start < NOW() - INTERVAL ? DAY", (2 * WINDOW_DAYS,))
                    conn.commit()
            finally:
                conn.close()
        except Exception as exc:  # noqa: BLE001 -- a record is never worth the supervisor's turn
            self.retry_after = self.monotonic() + BACKOFF_SECONDS
            log.swallowed("queue/pressurelearn", "_write", exc)

    def edges(self):
        """The edges the record gives, read again every `REFRESH_SECONDS`; the typed ones when learning is off, there is no record, or the database does not answer."""
        if not enabled():
            return typed()
        now = self.monotonic()
        if self.cached_at is not None and now - self.cached_at < REFRESH_SECONDS:
            return self.cached
        self.cached_at = now
        try:
            conn = self.connect()
            try:
                minutes, kills = read_record(conn)
            finally:
                conn.close()
            edges = learn(minutes, kills)
        except Exception as exc:  # noqa: BLE001 -- the typed edges stand
            log.swallowed("queue/pressurelearn", "edges", exc)
            edges = typed()
        self.cached = edges
        return edges


def read_record(conn):
    """`(minutes, kills)` of the last `WINDOW_DAYS` days: the `Minute`s of `pressure_minutes` and the times (seconds since the epoch) the ledger's memory kills ended."""
    minutes = [Minute(float(row["start"]), row["some_max"], row["full_max"], row["available_min_mb"], row["reserve_mb"], row["level_max"], row["held_seconds"], row["yields"], bool(row["psi"]))
               for row in conn.execute("SELECT UNIX_TIMESTAMP(minute_start) AS start, some_max, full_max, available_min_mb, reserve_mb, level_max, held_seconds, yields, psi FROM pressure_minutes"
                                       " WHERE minute_start >= NOW() - INTERVAL ? DAY ORDER BY minute_start", (WINDOW_DAYS,)).fetchall()]
    kills = [float(row["ended"]) for row in conn.execute("SELECT UNIX_TIMESTAMP(ended_at) AS ended FROM engine_attempts WHERE outcome = 'killed_memory' AND ended_at IS NOT NULL"
                                                         " AND ended_at >= NOW() - INTERVAL ? DAY ORDER BY ended_at", (WINDOW_DAYS,)).fetchall()]
    return minutes, kills
