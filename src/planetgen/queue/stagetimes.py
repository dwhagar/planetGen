# planetgen/queue/stagetimes.py

"""
What this server has learned about how long each kind of stage takes (PERF.184): when a stage of a queue run ends, its
wall seconds and its cost per worker second go into a Redis hash, one field per stage key, as a running average that
leans on the newest run. The next run of the stage starts from them (`runrecord.Stage.expected_seconds`, `worker_rate`),
so a stage that has no row in the control database's stage history (`generation.stats.stage_seconds`: the warm-up, the
fills) still has an estimate from its second run on, and the figure tightens as runs accumulate. A stage that has
never run has none, and `runeta` then falls back to `roughtime`.

Kept in Redis beside the run records (the same prefix, so a test's isolated prefix isolates it too): losing it only
returns the estimates to the untimed estimate (`roughtime`).
"""

import json

from planetgen.queue import runrecord

KEEP_WEIGHT = 0.5
"""float: The share of the stored average that stays when a new run is folded in (the newest run weighs the rest)."""

MIN_SECONDS = 0.5
"""float: A stage that ended faster than this tells nothing about a typical run and is not recorded."""


def times_key():
    return runrecord._k("stage-times")


def _fold(old, new):
    return new if old is None else KEEP_WEIGHT * old + (1.0 - KEEP_WEIGHT) * new


def record(conn, key, wall_seconds, cost_done, busy_seconds):
    """
    Folds one finished stage into `key`'s average: `wall_seconds` from begin to end, and (when units reported their
    time) the cost per worker second, `cost_done` over `busy_seconds`, the sum of the units' own seconds. A stage that
    took less than `MIN_SECONDS` is ignored.
    """
    if wall_seconds is None or wall_seconds < MIN_SECONDS:
        return
    raw = conn.hget(times_key(), key)
    old = json.loads(runrecord._s(raw)) if raw else {}
    entry = {"seconds": _fold(old.get("seconds"), float(wall_seconds)), "runs": int(old.get("runs") or 0) + 1}
    rate = (cost_done / busy_seconds) if busy_seconds and busy_seconds > 0 and cost_done and cost_done > 0 else None
    entry["worker_rate"] = _fold(old.get("worker_rate"), rate) if rate else old.get("worker_rate")
    conn.hset(times_key(), key, json.dumps(entry))


def lookup(conn, keys):
    """`{key: {"seconds", "worker_rate", "runs"}}` for the stage keys in `keys` that have been recorded."""
    found = {}
    for key in keys:
        raw = conn.hget(times_key(), key)
        if raw:
            found[key] = json.loads(runrecord._s(raw))
    return found
