# planetgen/queue/runview.py

"""
A run as the pages show it (PERF.99): `runrecord.read_run`'s dict turned into plain JSON-able numbers and words, in one
place, so the Queue page's API, the Generate page and the job log all say the same thing about the same run. Pure
functions of the run and the clock; the ETA is `runeta.estimate`'s.
"""

from planetgen.queue import runeta, runrecord

WAITING, RUNNING, DONE = "waiting", "running", "done"
"""str: A stage's state: not begun, begun and not ended, ended."""


def stage_state(stage):
    if stage["started"] is None:
        return WAITING
    return DONE if stage["ended"] is not None or (stage["units"] and stage["done"] >= stage["units"]) else RUNNING


def stage_view(stage, run_over, now):
    """One stage: its words and counters. `seconds` is how long it ran (so far, for a running one)."""
    state = stage_state(stage)
    seconds = None
    if stage["started"] is not None:
        end = stage["ended"] if stage["ended"] is not None else (now if not run_over else stage["last_end"] or stage["started"])
        seconds = max(end - stage["started"], 0.0)
    return {"key": stage["key"], "label": stage["label"], "state": state, "units": stage["units"], "done": stage["done"],
            "cost": stage["cost_total"], "cost_done": stage["cost_done"], "started": stage["started"],
            "ended": stage["ended"], "seconds": seconds, "longest_unit": stage.get("longest_unit"),
            "peak_mb": stage.get("peak_mb")}


def eta_view(run, now):
    """The time left in a live run, as the ETA reports it, or `None` for a run that is over."""
    if run["state"] not in runrecord.LIVE_STATES:
        return None
    eta = runeta.estimate(run, now)
    return {"seconds": eta.seconds, "low": eta.low, "high": eta.high, "source": eta.source,
            "from_earlier_runs": eta.from_earlier_runs, "partial": eta.partial, "stalled": eta.stalled,
            "rough": eta.rough}


def summary(run, now):
    """
    The run for a list and for the head of its page: identity, state, times, how far it is, the time left, whether it
    has stopped answering (`stale`), and the stage it is in. No stage detail (see `detail`).
    """
    live = run["state"] in runrecord.LIVE_STATES
    stages = run["stages"]
    cost_total = sum(stage["cost_total"] or 0.0 for stage in stages)
    cost_done = sum(stage["cost_done"] or 0.0 for stage in stages)
    current = stages[min(max(run["current"], 0), len(stages) - 1)] if stages else None
    started = run["started"] or run["created"]
    end = run["finished"] if run["finished"] is not None else (now if live else run["updated"])
    return {
        "id": run["id"], "title": run["title"], "kind": run["kind"], "who": run["who"] or None,
        "database": run["database"] or None, "job": run["job"], "argv": run["argv"], "state": run["state"],
        "error": run["error"],
        "live": live, "stale": runrecord.is_stale(run, now), "created": run["created"], "started": started,
        "updated": run["updated"], "finished": run["finished"],
        "seconds": None if started is None else max((end or now) - started, 0.0),
        "stage": None if current is None else {"index": stages.index(current), "count": len(stages),
                                               "label": current["label"], "state": stage_state(current)},
        "units": sum(stage["units"] for stage in stages), "units_done": sum(stage["done"] for stage in stages),
        "percent": round(100.0 * cost_done / cost_total, 1) if cost_total else (100.0 if run["state"] == "done" else None),
        "eta": eta_view(run, now),
    }


def detail(run, now):
    """`summary` and every stage."""
    over = run["state"] not in runrecord.LIVE_STATES
    return {**summary(run, now), "stages": [stage_view(stage, over, now) for stage in run["stages"]]}


def read(conn, run_id, now):
    """`summary` of the run `run_id`, or `None` when it is unknown."""
    run = runrecord.read_run(conn, run_id)
    return None if run is None else summary(run, now)


__all__ = ["DONE", "RUNNING", "WAITING", "detail", "eta_view", "read", "stage_state", "stage_view", "summary"]
