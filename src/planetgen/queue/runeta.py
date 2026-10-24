# planetgen/queue/runeta.py

"""
The time left in a run (PERF.96, PERF.33), as a pure function of the counters `runrecord.read_run` returns: no
Redis, no clock of its own (`now` is passed), so the terminal bar, the Generate and Queue pages, the banner
(UX.3) and a migration bar all show the same number and a test can drive it with a fake clock.

For each stage the estimate is, in the order of what is known (section 2.3 of
docs/design/operations/performance-eta-queue-and-caching.md):

1. a **rate**: the stage's live rate (cost finished per second, the decayed sums the units update) blended with
   the rate recorded for this kind of work, `progress_rate.blended_rate`; the time left is the cost left over the
   rate, counting down between updates, never below the time the rest takes once the next unit is overdue.
   Source `"live"`, or `"recorded"` while the recorded rate carries it.
2. **no unit has reported yet**, or the stage reports nothing at all (Boss, 2026-10-15 05:30Z: a phase with no
   updates still shows a time remaining): the time earlier runs took this stage (`expected_seconds`) less what
   has passed, held at a tenth of it while the stage runs over; source `"earlier runs"`.
3. no history either: the **time per cost unit so far in this run** applied to the cost left; source
   `"this run so far"`.
4. nothing to go on: unknown (`None`), shown as "estimating".

A stage that has not started contributes by 2 (or its cost at the recorded rate, or at the run's rate so far).
The run's time left is the sum; `partial` says a stage added nothing for want of an estimate.
"""

from dataclasses import dataclass

from planetgen.queue import progress_rate

OVERRUN_FRACTION = 0.1
"""float: A stage that outlasts what earlier runs took keeps this share of that time as its estimate."""

EARLIER_RUNS = "earlier runs"
THIS_RUN = "this run so far"


@dataclass
class Eta:
    """The estimate: `seconds` (`None` when unknown), the `low` and `high` ends of the likely range, `source`
    of the current stage's figure (`live`, `recorded`, `earlier runs`, `this run so far` or `None`),
    `from_earlier_runs` (the figure leans on earlier runs, not this one), `partial` (a stage had no estimate and
    was left out) and `stalled` (nothing finished for longer than a stall: the figure is held, not counted down).
    """

    seconds: float = None
    low: float = None
    high: float = None
    source: str = None
    from_earlier_runs: bool = False
    partial: bool = False
    stalled: bool = False


def _run_rate(run, now):
    """Cost finished per second over the whole run so far, or `None` before anything has finished."""
    started = run.get("started")
    done = sum(stage.get("cost_done") or 0.0 for stage in run["stages"])
    if started is None or not done or now <= started:
        return None
    return done / (now - started)


def _finished(stage):
    return stage["ended"] is not None or (stage["units"] and stage["done"] >= stage["units"])


def _stalled(stage, now):
    last = stage.get("last_end") or stage.get("started")
    if last is None or stage["done"] == 0:
        return False
    return now - last > max(progress_rate.STALL_FLOOR_SECONDS, 3.0 * (stage.get("longest_gap") or 0.0))


def _stage_remaining(stage, now, run_rate):
    """`(seconds or None, source or None, stalled)` for one stage."""
    if _finished(stage):
        return 0.0, None, False
    cost_left = max((stage["cost_total"] or 0.0) - stage["cost_done"], 0.0)
    prior = stage.get("prior_rate")
    expected = stage.get("expected_seconds")
    if stage["started"] is None:                                  # not begun
        if expected:
            return float(expected), EARLIER_RUNS, False
        if prior and cost_left:
            return cost_left / prior, "recorded", False
        if run_rate and cost_left:
            return cost_left / run_rate, THIS_RUN, False
        return None, None, False
    elapsed = max(now - stage["started"], 0.0)
    live = stage["u"] / stage["t"] if stage.get("t") else None
    rate = progress_rate.blended_rate(live, stage["done"], prior, elapsed)
    if rate:
        stalled = _stalled(stage, now)
        since = 0.0 if stalled or stage.get("last_end") is None else max(now - stage["last_end"], 0.0)
        mean_cost = (stage["cost_done"] / stage["done"]) if stage["done"] else (cost_left / max(stage["units"], 1))
        seconds = max(0.0, cost_left / rate - since, (cost_left - mean_cost) / rate)
        recorded = bool(prior) and (live is None or stage["done"] < progress_rate.HOLD_UNITS
                                    or elapsed < progress_rate.HOLD_SECONDS)
        return seconds, "recorded" if recorded else "live", stalled
    if expected:
        return max(expected - elapsed, OVERRUN_FRACTION * expected), EARLIER_RUNS, False
    if run_rate and cost_left:
        return cost_left / run_rate, THIS_RUN, False
    return None, None, False


def estimate(run, now):
    """
    The time left in `run` (`runrecord.read_run`'s dict) at `now`.

    Returns:
        Eta: `seconds` is 0 for a run that is over, `None` while nothing can be said.
    """
    if run["state"] not in ("running", "cancelling"):
        return Eta(seconds=0.0, low=0.0, high=0.0)
    run_rate = _run_rate(run, now)
    total, source, partial, stalled, current_seen = 0.0, None, False, False, False
    for index, stage in enumerate(run["stages"]):
        seconds, stage_source, stage_stalled = _stage_remaining(stage, now, run_rate)
        if seconds is None:
            partial = True
            continue
        total += seconds
        if not _finished(stage) and not current_seen:
            current_seen, source, stalled = True, stage_source, stage_stalled
    if partial and total == 0.0 and not all(_finished(stage) for stage in run["stages"]):
        return Eta(partial=True)
    low, high = progress_rate.eta_range(total) or (None, None)
    return Eta(seconds=total, low=low, high=high, source=source,
               from_earlier_runs=source in ("recorded", EARLIER_RUNS), partial=partial, stalled=stalled)


def public_view(run, now):
    """
    What a public read needs and nothing private (UX.3's banner, built on this): whether a run is going,
    which stage, how far, and the time left. No user, database or arguments.
    """
    if run is None or run["state"] not in ("running", "cancelling"):
        return {"active": False}
    eta = estimate(run, now)
    stages = run["stages"]
    current = min(max(run["current"], 0), len(stages) - 1)
    cost_total = sum(stage["cost_total"] or 0.0 for stage in stages)
    cost_done = sum(stage["cost_done"] for stage in stages)
    return {
        "active": True, "state": run["state"], "title": run["title"],
        "stage": {"index": current, "count": len(stages), "label": stages[current]["label"],
                  "done": stages[current]["done"], "units": stages[current]["units"]},
        "percent": round(100.0 * cost_done / cost_total, 1) if cost_total else None,
        "eta_seconds": eta.seconds, "eta_low": eta.low, "eta_high": eta.high, "eta_source": eta.source,
        "from_earlier_runs": eta.from_earlier_runs, "partial": eta.partial, "stalled": eta.stalled,
        "updated": run["updated"],
    }
