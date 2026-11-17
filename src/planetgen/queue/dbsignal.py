# planetgen/queue/dbsignal.py

"""
MySQL as a signal of the job engine (OPS.248, H3): admission slows when the latency of the engine's own statements says the database is the limit.

The memory controller (`pressure.Controller`, OPS.231) holds the bulk and upkeep lanes when the host is short of memory. This is that same controller with a second `look`: how long one of the
engine's own statements takes. Every few seconds a thread of its own sends the statement the supervisor sends most (the newest row of the ledger's `engine_jobs`, by its key) and times it; the
supervisor's loop never waits on MySQL, it reads the last figure (`Prober.latest`). A probe that has not come back counts as the time it has been out, so a database that has stopped answering reads
as slow at once, not when it finally answers.

**The reading.** The calm figure (`Baseline`) is learned: a running mean of the probes that were calm, so a fast database and a slow one each judge themselves. A probe is *pressed* at
`PRESSED_FACTOR` times the calm figure and no less than `PRESSED_FLOOR_MS`, *starved* at `STARVED_FACTOR` times and no less than `STARVED_FLOOR_MS`; until `BASELINE_PROBES` calm probes
are in, nothing is pressed (there is nothing to compare with yet), though a probe over the starved floor still counts. The floors keep a database that answers in a millisecond from reading as pressed at four.

**What it does.** Exactly what the memory controller does with a pressed host: after `pressure.HOLD_AFTER_SECONDS` in a row nothing of the bulk or maintenance lanes is admitted (the site's
jobs still go), and after `pressure.RELEASE_AFTER_SECONDS` of calm it admits again. It asks no job to stop (a slow database is not helped by killing a job), so it is made without a launcher.
Every change is said in the supervisor's list, as the memory controller's are. The factors and floors are design stand-ins; H4 (OPS.249) has them learned from recorded runs.
Nothing here changes what a run generates.
"""

import threading
import time
from typing import NamedTuple

from planetgen.queue import pressure
from planetgen.util import log

PROBE_SECONDS = 2.0
"""float: How often the prober sends its statement. Stand-in."""

PRESSED_FACTOR = 4.0
"""float: How many times the calm figure a probe takes before the database reads as pressed. Stand-in."""

STARVED_FACTOR = 10.0
"""float: How many times the calm figure before it reads as starved. Stand-in."""

PRESSED_FLOOR_MS = 50.0
"""float: The least a probe takes to read as pressed, however fast the database is when calm. Stand-in."""

STARVED_FLOOR_MS = 250.0
"""float: The least a probe takes to read as starved. Stand-in."""

BASELINE_PROBES = 5
"""int: How many calm probes the calm figure needs before it judges. Stand-in."""

BASELINE_WEIGHT = 0.1
"""float: The weight of each new calm probe in the calm figure (a running mean). Stand-in."""

PROBE_SQL = "SELECT id FROM engine_jobs ORDER BY id DESC LIMIT 1"
"""str: The statement timed: the newest job row by its key, the shape of what the supervisor's ledger reads send."""

WHY = "the database is the limit"


class Reading(NamedTuple):
    """One look: the probe's time and the calm figure it was judged against (ms), and the `level` (`pressure.CALM`, `PRESSED` or `STARVED`)."""

    latency_ms: float
    baseline_ms: float
    level: int

    def words(self):
        return f"a ledger read took {self.latency_ms:.0f} ms where it takes {self.baseline_ms:.0f} ms when calm"


class Baseline:
    """The calm figure: a running mean of the calm probes, learned over the first `BASELINE_PROBES`."""

    def __init__(self):
        self.mean_ms = None
        self.count = 0

    def ready(self):
        return self.count >= BASELINE_PROBES

    def level(self, latency_ms, learn=True):
        """The level of a probe against what has been learned; a calm probe adds to what is learned (`learn` is false for a probe that was already counted)."""
        if self.ready():
            if latency_ms >= max(STARVED_FLOOR_MS, STARVED_FACTOR * self.mean_ms):
                return pressure.STARVED
            if latency_ms >= max(PRESSED_FLOOR_MS, PRESSED_FACTOR * self.mean_ms):
                return pressure.PRESSED
        elif latency_ms >= STARVED_FLOOR_MS:
            return pressure.STARVED
        if not learn:
            return pressure.CALM
        self.mean_ms = latency_ms if self.mean_ms is None else (1.0 - BASELINE_WEIGHT) * self.mean_ms + BASELINE_WEIGHT * latency_ms
        self.count += 1
        return pressure.CALM


class Prober:
    """
    The thread that times the probe statement. `connect` opens the connection it uses (its own, so the supervisor's never waits); `clock` gives the time. `latest()` is the last figure in ms,
    or the time the probe in flight has been out when that is longer, or `None` before the first probe.
    """

    def __init__(self, connect, clock=time.monotonic, sql=PROBE_SQL):
        self.connect, self.clock, self.sql = connect, clock, sql
        self.last_ms = None
        self.sent_at = None
        self.completed = 0      # how many probes have come back
        self.stop = threading.Event()
        self.thread = None

    def probe(self, conn):
        """Sends the statement once and records how long it took."""
        self.sent_at = self.clock()
        conn.execute(self.sql).fetchone()
        self.last_ms = (self.clock() - self.sent_at) * 1000.0
        self.sent_at = None
        self.completed += 1

    def latest(self):
        out_ms = None if self.sent_at is None else (self.clock() - self.sent_at) * 1000.0
        if out_ms is not None and (self.last_ms is None or out_ms > self.last_ms):
            return out_ms
        return self.last_ms

    def run(self):
        conn = None
        while not self.stop.is_set():
            try:
                conn = conn or self.connect()
                self.probe(conn)
            except Exception as exc:  # noqa: BLE001 -- a database that cannot answer is what is measured: the probe stays out (it reads as slow by the time it has been out) until a new connection answers
                log.swallowed("queue/dbsignal", "probe", exc)
                conn = None
                self.sent_at = self.sent_at if self.sent_at is not None else self.clock()
            self.stop.wait(PROBE_SECONDS)

    def start(self):
        self.thread = threading.Thread(target=self.run, name="dbsignal-probe", daemon=True)
        self.thread.start()
        return self

    def close(self):
        self.stop.set()


class Look:
    """What the controller calls each turn: the prober's latest figure judged against the learned calm figure. No probe yet reads as calm and teaches nothing."""

    def __init__(self, prober, baseline=None):
        self.prober, self.baseline = prober, baseline or Baseline()
        self.counted = 0        # the probes already taught to the calm figure: the controller looks more often than the prober probes

    def __call__(self):
        latency = self.prober.latest()
        if latency is None:
            return Reading(0.0, 0.0, pressure.CALM)
        level = self.baseline.level(latency, learn=self.prober.completed != self.counted)
        self.counted = self.prober.completed
        return Reading(latency, self.baseline.mean_ms or 0.0, level)


def controller(prober, clock=time.time, baseline=None):
    """The database's controller: the memory controller's hold and release with this signal for its `look`, and no job to ask to stop."""
    return pressure.Controller(launcher=None, look=Look(prober, baseline), clock=clock, why=WHY)
