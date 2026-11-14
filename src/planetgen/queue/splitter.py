# planetgen/queue/splitter.py

"""
The queue engine's choice of a job's split (OPS.161, docs/design/operations/dynamic-queue-chain.md).

A stage's cut is made once, when the stage opens (the planner's moment, `wrapper._open`), and the decision is written into the stage's plan in the run record, so a restarted
unit uses the same cut and a run can be read back from its record. The choice is a pure function, `choose`:

* **A job flagged `SUGGESTED`, or one that declares nothing, keeps its suggested split exactly** (`breakpoints`). That is every job today.
* **A job flagged `BYTE_IDENTICAL` gets a cut from the profile and the free resources** (`splitprofile.lookup`, `memorygate.read_host`): about the number of items that take the
  target seconds of a unit (`TARGET_SHARE` of `queue.unit_ceiling_seconds`), fewer if the unit's memory (`mb_base + mb_per_item x items`) would not fit its share of the free
  memory (`MEMORY_SHARE`, over the workers the stage may have), between the job's floor and ceiling (the typed split stays the default and the bound), cut only at legal
  positions (a span that must stay whole is never cut inside) and into at most `units.MAX_UNITS` units. With no usable profile the suggested split is used and the decision says so.

Every decision is a `Decision` that says what was chosen and why, with the suggested figures beside the new ones; `announce` writes it to the run's log (Boss, 3:03 PM EST on
2026-11-12: "ALL of this is logged including decisions when the queue changes values").
"""

import bisect
import math
from dataclasses import dataclass, field

from planetgen.queue import breakpoints, splitprofile
from planetgen.util import log

TARGET_SHARE = 0.2
"""float: The share of `queue.unit_ceiling_seconds` a unit is aimed at (300 s ceiling, 60 s unit)."""

MEMORY_SHARE = 0.5
"""float: The share of the host's free memory the units of a stage may hold at once; the rest stays for the other lanes and the site (the gate's `STAGE_MEMORY_SHARE`)."""

BY_SUGGESTED = "suggested"
BY_PROFILE = "profile"


@dataclass(frozen=True)
class Resources:
    """What the stage may use now: the host's free memory (MB), its cores, and the workers the stage may have at once."""

    available_mb: float
    cores: int
    workers: int


@dataclass(frozen=True)
class Decision:
    """
    The engine's choice for one job: `cuts` (the positions the job is cut at), `by` (`BY_SUGGESTED` or `BY_PROFILE`), the `flag` the job declared and the `reason` in words.
    `old_*` are the job's suggested split, `new_*` the one chosen (equal when `by` is suggested); `samples` is how many finished units the profile stood on (0 for none).
    """

    kind: str
    flag: str
    by: str
    reason: str
    items: int
    cuts: tuple
    old_units: int
    old_items: int
    new_units: int
    new_items: int
    samples: int = 0
    seconds_per_item: float = None
    limit: str = None
    resources: tuple = field(default=None)

    def as_dict(self, cuts=None):
        """The decision as plain JSON for the run record: the cuts only when they were chosen (`by` profile) or `cuts` is true; a suggested split is the job's own."""
        keep = (self.by == BY_PROFILE) if cuts is None else cuts
        record = {"kind": self.kind, "flag": self.flag, "by": self.by, "reason": self.reason, "items": self.items,
                  "old_units": self.old_units, "old_items": self.old_items, "new_units": self.new_units, "new_items": self.new_items,
                  "samples": self.samples, "seconds_per_item": self.seconds_per_item, "limit": self.limit,
                  "resources": None if self.resources is None else list(self.resources)}
        if keep:
            record["cuts"] = list(self.cuts)
        return record

    def message(self):
        """The one line of the run's log for this decision."""
        if self.by == BY_SUGGESTED:
            return f"Split of {self.kind}: {self.new_units:,} unit(s) of about {self.new_items:,} item(s), the job's suggested split: {self.reason}."
        return (f"Split of {self.kind} changed from {self.old_units:,} unit(s) of about {self.old_items:,} item(s) to {self.new_units:,} of about "
                f"{self.new_items:,}: {self.reason} ({self.samples} sample(s)).")


def resources_now(kind=None):
    """
    The `Resources` of the host now: `memorygate.read_host` for the free memory and the cores, and for the workers the stage may have the smaller of the cores and
    `workers.stage_worker_cap(kind)`. A reading that fails gives a single worker on 1 core with no memory to spare, so the suggested figures stand.
    """
    from planetgen.queue import memorygate, workers
    try:
        host = memorygate.read_host()
        cap = workers.stage_worker_cap(kind) if kind else None
        return Resources(host.available_mb, host.cores, max(1, min(cap or host.cores, host.cores)))
    except Exception as exc:  # noqa: BLE001 -- learning never blocks a job; the typed figures stand
        log.swallowed("queue/splitter", "resources_now", exc)
        return Resources(0.0, 1, 1)


def target_seconds():
    """The unit length the engine aims at, seconds: `TARGET_SHARE` of the unit ceiling setting."""
    from planetgen.util.settings import get_settings
    return TARGET_SHARE * get_settings().queue.unit_ceiling_seconds


def _suggested(kind, decl, reason, reading=None, resources=None):
    sizes = decl.sizes()
    usual = decl.suggested_unit_items()
    return Decision(kind, decl.split, BY_SUGGESTED, reason, decl.items, decl.suggested, len(sizes), usual, len(sizes), usual,
                    samples=0 if reading is None else reading.samples, seconds_per_item=None if reading is None else reading.seconds_per_item,
                    resources=None if resources is None else (round(resources.available_mb), resources.cores, resources.workers))


def _snap(decl, per_unit):
    """The cuts that put about `per_unit` items in a unit, at legal positions only: the largest legal cut at or before the aim (so the unit stays under it), else the next
    one after it (a span that stays whole is longer than the aim)."""
    legal = decl.legal_cuts
    cuts, start = [], 0
    while start + per_unit < decl.items:
        aim = start + per_unit
        index = bisect.bisect_right(legal, aim) - 1
        cut = legal[index] if index >= 0 and legal[index] > start else None
        if cut is None:
            after = bisect.bisect_right(legal, start)
            if after >= len(legal):
                break
            cut = legal[after]
        cuts.append(cut)
        start = cut
    return tuple(cuts)


def choose(kind, decl, reading=None, resources=None, seconds=None, max_units=None):
    """
    The split of job `kind` that declares `decl` (`breakpoints.Declaration`), given the profile's `reading` (`splitprofile.Reading`, or `None`), the `resources` and the unit
    length `seconds` to aim at. Pure: it reads nothing and changes nothing.

    Returns:
        Decision: The cut and why. Always legal: every cut is one of `decl.legal_cuts`.
    """
    if max_units is None:
        from planetgen.queue import units
        max_units = units.MAX_UNITS
    if decl.items == 0 or not decl.legal_cuts:
        return _suggested(kind, decl, "the job is one unit by nature" if decl.items else "the job has no item", reading, resources)
    if not decl.dynamic:
        return _suggested(kind, decl, "its output depends on the cut, so its suggested split is used", reading, resources)
    if reading is None or not reading.usable():
        return _suggested(kind, decl, "the profile has nothing usable for the kind yet", reading, resources)
    if resources is None:
        return _suggested(kind, decl, "the machine's resources were not given", reading)
    seconds = seconds if seconds else target_seconds()
    by_time = seconds / reading.seconds_per_item
    budget = resources.available_mb * MEMORY_SHARE / max(resources.workers, 1)
    by_memory = math.inf
    if reading.mb_per_item > 0:
        by_memory = max(budget - reading.mb_base, 0.0) / reading.mb_per_item
    want = min(by_time, by_memory)
    limit = "memory" if by_memory < by_time else "time"
    floor, ceiling = decl.clamp()
    per_unit = int(max(1, want))
    if per_unit < floor:
        per_unit, limit = floor, "floor"
    elif per_unit > ceiling:
        per_unit, limit = ceiling, "ceiling"
    spread = -(-decl.items // max(resources.workers, 1))
    if per_unit > spread and spread >= floor:
        per_unit, limit = spread, "workers"
    per_unit = max(per_unit, -(-decl.items // max_units))
    cuts = _snap(decl, per_unit)
    while len(cuts) + 1 > max_units:
        per_unit = int(per_unit * 1.25) + 1
        cuts = _snap(decl, per_unit)
    sizes = decl.sizes(cuts)
    usual = int(sorted(sizes)[len(sizes) // 2]) if sizes else 0
    old_sizes = decl.sizes()
    shown = (round(resources.available_mb), resources.cores, resources.workers)
    reason = (f"about {seconds:,.0f} s a unit at {reading.seconds_per_item * 1000:,.2f} ms an item on {resources.workers} worker(s) with {resources.available_mb:,.0f} MB free; "
              f"limited by {limit}")
    return Decision(kind, decl.split, BY_PROFILE, reason, decl.items, cuts, len(old_sizes), decl.suggested_unit_items(), len(sizes), usual,
                    samples=reading.samples, seconds_per_item=reading.seconds_per_item, limit=limit, resources=shown)


def decide(kind, decl, resources=None):
    """
    `choose` with the stub or stored profile (`splitprofile.lookup`) and the host now (`resources_now`): the call a command's `expand` makes to cut its stage (PERF.261).
    Any failure to read either leaves the suggested split, and says so.
    """
    try:
        reading = splitprofile.lookup(kind)
    except Exception as exc:  # noqa: BLE001 -- learning never blocks a job
        log.swallowed("queue/splitter", "decide", exc)
        reading = None
    if not decl.dynamic:
        return choose(kind, decl, reading)
    return choose(kind, decl, reading, resources or resources_now(kind))


def announce(conn, run_id, stage, decision):
    """Writes `decision` to the run's log (`runrecord.log`): the line of `Decision.message` with the old and new figures, the samples and the reason as fields."""
    from planetgen.queue import runrecord
    level = "info"
    runrecord.log(conn, run_id, decision.message(), level=level, kind="split.decided", stage=stage, job=decision.kind, split_flag=decision.flag, by=decision.by,
                  old_units=decision.old_units, new_units=decision.new_units, old_items=decision.old_items, new_items=decision.new_items,
                  samples=decision.samples, limit=decision.limit or "")


def spans(items, decision):
    """The `(start, end)` item range of each unit of `decision` over `items` (what `Declaration.spans` gives for the chosen cuts)."""
    edges = [0, *decision.cuts, items]
    return list(zip(edges, edges[1:])) if items else []
