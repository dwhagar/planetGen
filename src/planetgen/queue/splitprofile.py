# planetgen/queue/splitprofile.py

"""
What the queue engine knows of a kind of job before it cuts it (OPS.161): seconds and memory per item, from the stored profile.

**Reader of the OPS.155 rows** (`jobprofiles.unit_rows`, keyed `unit:<command>`). The engine uses the unit's own seconds per item and a fit of its peak memory against its items. Rows from this machine's
shape and generator epoch (`profileclaim._shape_of`, the one rule of what the same machine is) count in full when there are at least `MIN_SAMPLES` of them; otherwise rows from elsewhere count at half
trust (`OTHER_TRUST`, the design note's half weight); with fewer than `MIN_SAMPLES` of either the reading is `None` and the job keeps its suggested split. Any failure to read gives `None`.
A test (or a lane before the rows exist) puts readings in `STUB`, which wins over the table.
"""

from dataclasses import dataclass

from planetgen.util import log

MIN_SAMPLES = 3
"""int: The fewest samples a reading stands on before the engine trusts it (the design note's fallback chain: a profile with at least 3 samples on this machine shape)."""


@dataclass(frozen=True)
class Reading:
    """
    What the profile holds for one kind of job (`planetgen/queue/splitprofile.py`, docs/design/operations/self-tuning-jobs.md, "One stored profile per job kind"):
    `seconds_per_item` the unit's own worked seconds for one item; `mb_base` the memory a unit holds before its first item (the worker and its tables) and `mb_per_item` the
    slope of the peak memory per item (with samples from two or more size buckets the line through the lowest and highest bucket's medians; with one bucket the raw peak over the items, both the bucket's medians, and no base, which overstates the slope and so never cuts a unit too large); `samples` how many finished units this stands on; `trust` 1.0 for this machine's shape
    and this epoch, 0.5 for a profile measured elsewhere (the half trust of the design note, `OTHER_TRUST`).
    """

    seconds_per_item: float
    mb_base: float = 0.0
    mb_per_item: float = 0.0
    samples: int = 0
    trust: float = 1.0

    def usable(self):
        """Whether the engine may choose a cut from this: enough samples, a positive time per item and some trust."""
        return self.samples >= MIN_SAMPLES and self.seconds_per_item > 0 and self.trust > 0


OTHER_TRUST = 0.5
"""float: The weight of a profile measured on another machine shape or epoch."""

STUB = {}
"""dict: Readings by job kind, for a test; empty in production."""


def _median(values):
    ordered = sorted(values)
    middle = len(ordered) // 2
    return ordered[middle] if len(ordered) % 2 else 0.5 * (ordered[middle - 1] + ordered[middle])


def _memory(rows):
    """`(mb_base, mb_per_item)` from the rows that hold a peak and items: the line through the lowest and highest size bucket's medians, or the raw peak over the items for one bucket."""
    buckets = {}
    for row in rows:
        if row["peak_mb"] and row["items"] and row["items"] > 0:
            buckets.setdefault(row["size_bucket"], []).append((row["items"], row["peak_mb"]))
    if not buckets:
        return 0.0, 0.0
    points = {bucket: (_median([i for i, _ in pairs]), _median([m for _, m in pairs])) for bucket, pairs in buckets.items()}
    low, high = points[min(points)], points[max(points)]
    if len(points) >= 2 and high[0] > low[0]:
        slope = (high[1] - low[1]) / (high[0] - low[0])
        if slope > 0:
            return max(low[1] - slope * low[0], 0.0), slope
    return 0.0, max(m / i for i, m in points.values())


def from_rows(rows, host_memory, host_cores, epoch):
    """The `Reading` of one kind from its `jobprofiles.unit_rows`, or `None` when neither this machine's rows nor the rest reach `MIN_SAMPLES`."""
    from planetgen.queue import profileclaim
    rows = [r for r in rows if r["seconds_per_item"] and r["seconds_per_item"] > 0]
    mine = [r for r in rows if profileclaim._shape_of(r, host_memory, host_cores, epoch)]
    used, trust = (mine, 1.0) if len(mine) >= MIN_SAMPLES else (rows, OTHER_TRUST)
    if len(used) < MIN_SAMPLES:
        return None
    base, slope = _memory(used)
    return Reading(_median([r["seconds_per_item"] for r in used]), base, slope, len(used), trust)


def lookup(kind):
    """The `Reading` of the units of command `kind`, or `None` when the profile holds nothing usable for it (a `STUB` entry first, then the `unit:<kind>` rows of `jobprofiles`)."""
    if kind in STUB:
        return STUB[kind]
    try:
        from planetgen.queue import jobprofiles, memorygate, profileclaim
        from planetgen.util import generator_epoch
        if not (jobprofiles.enabled() and profileclaim.enabled()):       # the one switch of the self-tuning readers (OPS.156) and the recording switch
            return None
        host = memorygate.read_host()
        return from_rows(jobprofiles.unit_rows("unit:" + kind), host.total_mb, host.cores, generator_epoch.generator_epoch())
    except Exception as exc:  # noqa: BLE001 -- learning never blocks a job; the suggested split stands
        log.swallowed("queue/splitprofile", "lookup", exc)
        return None
