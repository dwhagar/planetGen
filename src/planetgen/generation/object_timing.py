# planetgen/generation/object_timing.py

"""
How long each generated object takes, per object and per object type
(PERF.92), on top of the per-stage and per-job records (PERF.67).

Every generator reports the same way: a generated class is registered with
`track(cls, kind)`, which wraps its constructor, so no generator carries
timing code of its own. While a *unit of work* (a sector, a run of scatter shells)
runs inside `unit(...)`, each object built is timed; the time an object
spends building the objects inside it (a planet's moons, a system's stars)
is taken out of its own, so the per-kind totals add up to the unit's time
instead of counting nested objects twice.

At billions of objects a row each is too much (Boss, 2026-10-14), so the
unit keeps only, per kind: count, total, minimum and maximum seconds, saved as
one `object_timing_stats` row per kind and sector. (PERF.197 took out the
histogram and the repeatable sample of single objects: nothing read them.)

Times here are wall-clock seconds in a worker, like the rest of the stats.
Nothing in this module can fail a run: saving is skipped without a run id
or when the stats are switched off (`PLANETGEN_GENERATION_STATS=0`).
"""

import contextlib
import functools
import math
import os
import threading
import time

KINDS = (
    "system", "star", "planet", "moon", "belt", "comet", "binary-star", "asteroid-field", "black-hole",
    "neutron-star", "nebula", "quasar", "rogue-planet", "interstellar-comet", "supernova-remnant",
)
"""tuple: Every object kind the generators create. A test checks that
each one reports and that every registered class names one of these."""

_KIND_OF = {}
"""dict: Class -> kind, filled by `track`."""

_local = threading.local()


def kind_of(cls):
    """The kind registered for `cls` or the nearest base class that has
    one, or `None`."""
    for base in cls.__mro__:
        if base in _KIND_OF:
            return _KIND_OF[base]
    return None


def stats_enabled():
    """Whether timings may be saved (`PLANETGEN_GENERATION_STATS` is not off)."""
    return os.environ.get("PLANETGEN_GENERATION_STATS", "1").strip().lower() not in ("0", "off", "no", "false")


class Aggregate:
    """One kind's count, total and extremes, in seconds."""

    __slots__ = ("count", "total", "minimum", "maximum")

    def __init__(self):
        self.count = 0
        self.total = 0.0
        self.minimum = math.inf
        self.maximum = 0.0

    def add(self, seconds):
        self.count += 1
        self.total += seconds
        if seconds < self.minimum:
            self.minimum = seconds
        if seconds > self.maximum:
            self.maximum = seconds

    @property
    def mean(self):
        return self.total / self.count if self.count else 0.0


class Collector:
    """The timings of one unit of work: per-kind aggregates."""

    def __init__(self, label):
        self.label = label
        self.aggregates = {}
        self._stack = []
        self._open = set()

    def begin(self, obj, kind):
        """Starts timing `obj`; returns the token `end` takes, or `None`
        when `obj` is already being timed (a subclass constructor calling
        its base's)."""
        key = id(obj)
        if key in self._open:
            return None
        self._open.add(key)
        frame = [kind, time.perf_counter(), 0.0, key]
        self._stack.append(frame)
        return frame

    def end(self, frame, kind=None):
        """Stops timing: takes the object's time off its parent's own time
        and records it (under `kind` when the finished object says what it
        turned out to be)."""
        total = time.perf_counter() - frame[1]
        self._stack.pop()
        self._open.discard(frame[3])
        own = max(total - frame[2], 0.0)
        if self._stack:
            self._stack[-1][2] += total
        kind = kind or frame[0]
        aggregate = self.aggregates.get(kind)
        if aggregate is None:
            aggregate = self.aggregates[kind] = Aggregate()
        aggregate.add(own)


def active():
    """The current unit's collector, or `None` outside a unit."""
    return getattr(_local, "collector", None)


@contextlib.contextmanager
def unit(label):
    """Times the objects built inside the block as one unit of work named
    `label` (a sector's address, a run of scatter shells); yields its `Collector`.
    Units don't nest: a unit inside another shares the outer's collector."""
    outer = active()
    if outer is not None:
        yield outer
        return
    collector = Collector(label)
    _local.collector = collector
    try:
        yield collector
    finally:
        _local.collector = None


def track(cls, kind, kind_after=None):
    """Registers `cls` as generating objects of `kind` and wraps its
    constructor so each one built inside a `unit` is timed. `kind_after`,
    when given, is called with the finished object and returns its kind (a
    `Planet` is a moon once it knows it is one). Returns `cls`."""
    if kind not in KINDS:
        raise ValueError(f"unknown object kind {kind!r}")
    _KIND_OF[cls] = kind
    original = cls.__init__

    @functools.wraps(original)
    def timed_init(self, *args, **kwargs):
        collector = getattr(_local, "collector", None)
        if collector is None:
            return original(self, *args, **kwargs)
        frame = collector.begin(self, kind_of(type(self)) or kind)
        if frame is None:
            return original(self, *args, **kwargs)
        try:
            return original(self, *args, **kwargs)
        finally:
            collector.end(frame, kind_after(self) if kind_after is not None and type(self) is cls else None)

    cls.__init__ = timed_init
    return cls


def save(conn, run_id, collector, sector_id=None):
    """Writes a finished unit's aggregates (the caller commits,
    so they land with the sector). Does nothing without a run id or when the
    stats are off. Returns the number of rows written."""
    if run_id is None or not stats_enabled():
        return 0
    rows = 0
    for kind, aggregate in sorted(collector.aggregates.items()):
        conn.execute(
            "INSERT INTO object_timing_stats (run_id, sector_id, kind, object_count, total_us, min_us, max_us)"
            " VALUES (?, ?, ?, ?, ?, ?, ?)",
            (run_id, sector_id, kind, aggregate.count, round(aggregate.total * 1e6),
             round(aggregate.minimum * 1e6), round(aggregate.maximum * 1e6)))
        rows += 1
    return rows


def kind_summaries(conn):
    """
    What each kind of object cost over every saved run (PERF.196): one `{"kind", "objects", "total_seconds", "mean_seconds",
    "fastest_seconds", "slowest_seconds"}` per kind, the kind that took the most time first. The admin Stats page shows it (the Seconds per
    object table) through `GET /api/admin/generation-stats`.
    """
    cursor = conn.execute(
        "SELECT kind, SUM(object_count) AS objects, SUM(total_us) AS micros, MIN(min_us) AS fastest, MAX(max_us) AS slowest"
        " FROM object_timing_stats GROUP BY kind")
    rows = [{"kind": row["kind"], "objects": int(row["objects"]), "total_seconds": float(row["micros"]) / 1e6,
             "mean_seconds": float(row["micros"]) / 1e6 / float(row["objects"]),
             "fastest_seconds": float(row["fastest"]) / 1e6, "slowest_seconds": float(row["slowest"]) / 1e6}
            for row in cursor.fetchall() if row["objects"]]
    return sorted(rows, key=lambda row: -row["total_seconds"])
