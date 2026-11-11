# planetgen/queue/unitbeat.py

"""
The beat of a running unit (OPS.136; reports/heartbeat-research-2026-11-09.md, docs/design/operations/redis-and-the-work-queue.md 5.16): what a horse writes about the unit it is running, so
the follower can judge it (`queue/beatjudge.py`, OPS.137) and a person can see it.

* **The work loop only sets numbers.** The `report(done, estimate)` a unit calls (as often as it likes, every 0.25 s or every item) records the counter in this process and writes nothing to Redis;
  `phase(text)` and `waiting(tag)` do the same for what it is doing and what it waits for (`gate`, `sky lock`, `mysql`).
* **A thread writes.** A thread of the horse writes the record into the stage's progress hash (the field of the unit, no new Redis structure) every `beatjudge.BEAT_SECONDS` (10 s), and sooner,
  after `COUNTER_WRITE_SECONDS`, when the counter has moved since the last write, so the bar of a run keeps its pace. It beats while the work loop is frozen, which is what tells a stuck unit
  (fresh beat, counter unchanged) from a dead one (no beat).
* **The clock is Redis's.** `beat`, `started` and `moved` are seconds on the Redis server clock (`TIME`), so a lane machine and the server cannot disagree; `moved` is worked out from this process's
  monotonic clock when the counter changed, set against the server time of the write.
* **The record** is one JSON object in the unit's field: `done`, `estimate`, `beat`, `started`, `moved` (absent until `done` first changes), `phase`, `waiting_on`, `restarts` (how many times the
  current phase label has begun again) and `resident_mb` (optional). `read` gives a `beatjudge.Beat` for each unit of a stage.
"""

import json
import threading
import time
from contextlib import contextmanager

from planetgen.queue import beatjudge, runrecord

COUNTER_WRITE_SECONDS = 2.0
"""float: A counter that moved is written at most this often, between two beats."""

TICK_SECONDS = 0.5
"""float: How often the thread looks at whether a write is due (nothing is sent to Redis on a tick that writes nothing)."""


def server_time(conn):
    """Seconds on the Redis server clock."""
    seconds, microseconds = conn.time()
    return float(seconds) + float(microseconds) / 1_000_000.0


class UnitBeat:
    """
    The record of one running unit and the thread that writes it. Used as the `report` of a unit (`report(done, estimate)`); a context manager that starts the thread on entry and stops it,
    after a last write, on exit.
    """

    def __init__(self, conn, run_id, stage, unit, clock=time.monotonic, beat_seconds=None, counter_write_seconds=None, resident=None, tick_seconds=TICK_SECONDS):
        self.conn, self.key, self.unit = conn, runrecord.progress_key(run_id, stage), str(unit)
        self.clock = clock
        self.beat_seconds = beatjudge.BEAT_SECONDS if beat_seconds is None else beat_seconds
        self.counter_write_seconds = COUNTER_WRITE_SECONDS if counter_write_seconds is None else counter_write_seconds
        self.resident = resident
        self.tick_seconds = tick_seconds
        self._lock = threading.Lock()
        self._stop = threading.Event()
        self._thread = None
        self._done = self._estimate = 0.0
        self._moved_at = None                       # this process's monotonic clock when `done` last changed
        self._phase, self._waiting, self._phases = "", "", {}
        self._started_at = self.clock()
        self._started_server = None
        self._last_write = None
        self._written_done = None

    # --- what the work loop sets: numbers in this process, nothing sent ---------------------------------------

    def __call__(self, done, estimate):
        with self._lock:
            if float(done) != self._done:
                self._moved_at = self.clock()
            self._done, self._estimate = float(done), float(estimate)

    def phase(self, text):
        """The unit begins `text`; the same label beginning again counts as a restart of the phase."""
        with self._lock:
            self._phases[text] = self._phases.get(text, 0) + 1
            self._phase = text

    @contextmanager
    def waiting(self, tag):
        """The unit waits for `tag` (`gate`, `sky lock` or `mysql`) inside the block."""
        with self._lock:
            before, self._waiting = self._waiting, tag
        try:
            yield
        finally:
            with self._lock:
                self._waiting = before

    # --- the record and its thread ---------------------------------------------------------------------------

    def record(self, now_server=None):
        """The record as it stands, a dict ready to be written; `now_server` is the Redis server time of the write."""
        now_server = server_time(self.conn) if now_server is None else now_server
        if self._started_server is None:
            self._started_server = now_server - (self.clock() - self._started_at)
        with self._lock:
            made = {"done": self._done, "estimate": self._estimate, "beat": now_server, "started": self._started_server,
                    "phase": self._phase, "waiting_on": self._waiting, "restarts": max(0, self._phases.get(self._phase, 1) - 1)}
            if self._moved_at is not None:
                made["moved"] = now_server - (self.clock() - self._moved_at)
            written = self._done
        if self.resident is not None:
            made["resident_mb"] = round(float(self.resident()), 1)
        self._written_done = written
        return made

    def write(self):
        """One write of the record to the stage's progress hash."""
        self.conn.hset(self.key, self.unit, json.dumps(self.record()))
        self._last_write = self.clock()

    def _due(self):
        if self._last_write is None:
            return True
        since = self.clock() - self._last_write
        if since >= self.beat_seconds:
            return True
        with self._lock:
            moved = self._done != self._written_done
        return moved and since >= self.counter_write_seconds

    def _run(self):
        while not self._stop.is_set():
            if self._due():
                try:
                    self.write()
                except Exception as exc:        # a Redis that cannot be reached must not end good work (the follower judges only fresh reads)
                    from planetgen.util import log
                    log.trace(f"unit {self.unit}: the beat could not be written: {type(exc).__name__}: {exc}", area="queue")
            self._stop.wait(self.tick_seconds)

    def __enter__(self):
        global _current
        self._outer, _current = _current, self
        self._thread = threading.Thread(target=self._run, name=f"beat-{self.unit}", daemon=True)
        self._thread.start()
        return self

    def __exit__(self, *exc_info):
        global _current
        _current = self._outer
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=5.0)
        return False


_current = None         # the beat of the unit this process is running (one unit per horse), for code too deep to be handed the reporter


@contextmanager
def waiting(tag):
    """The running unit waits for `tag` inside the block; nothing happens outside a unit (a command run by hand)."""
    if _current is None:
        yield
        return
    with _current.waiting(tag):
        yield


def phase(text):
    """The running unit begins `text`; nothing happens outside a unit."""
    if _current is not None:
        _current.phase(text)


def parse(raw):
    """A stored field as a dict, `None` for one that is not a record."""
    try:
        found = json.loads(raw if isinstance(raw, str) else raw.decode())
    except (ValueError, AttributeError, TypeError):
        return None
    return found if isinstance(found, dict) else None


def beat_of(found):
    """The `beatjudge.Beat` of a stored record."""
    return beatjudge.Beat(beat=float(found["beat"]), started=float(found.get("started", found["beat"])), done=float(found.get("done", 0.0)),
                          estimate=float(found.get("estimate") or 0.0), moved=None if found.get("moved") is None else float(found["moved"]),
                          phase=str(found.get("phase", "")), waiting_on=str(found.get("waiting_on", "")), restarts=int(found.get("restarts", 0)),
                          resident_mb=None if found.get("resident_mb") is None else float(found["resident_mb"]))


def read(conn, run_id, stage):
    """`{unit: Beat}` of the units of a stage that have a record."""
    found = {}
    for unit, raw in conn.hgetall(runrecord.progress_key(run_id, stage)).items():
        record = parse(raw)
        if record is None or "beat" not in record:
            continue
        try:
            found[int(unit)] = beat_of(record)
        except (ValueError, TypeError, KeyError):
            continue
    return found
