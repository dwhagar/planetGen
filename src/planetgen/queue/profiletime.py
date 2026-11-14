# planetgen/queue/profiletime.py

"""
The expected seconds of a unit, from the job profile (OPS.159, step 6 of `docs/design/operations/self-tuning-jobs.md`).

Three readers want to know how long a unit of a command takes before this run has told them: the heartbeat judge (`commandgaps.history_of`: `beatjudge.History.expected_seconds`, with the
unit's estimate it gives the rate the unit should show), the ETA of a stage nothing has timed (`runeta._rough_remaining` through `roughtime.stage_seconds`) and the stage's own expected time.
Until now they had the Redis running averages (`stagetimes`, `unitrates`) and, for the untimed, `roughtime.UNIT_SECONDS`, a typed 20 s. This reads the median worked seconds of the
kind's last `profileclaim.WINDOW` finished units on this machine (`jobprofiles` rows), by the same rules as the memory claim (`profileclaim`):

* at least `MIN_SAMPLES` samples from this machine shape and generator epoch; samples from another shape or epoch count at half weight (the mean of their median and the typed figure, or
  their median alone where there is none); fewer than that, no figure and the readers keep what they had (`None`);
* worked out once per run (the run's id is the token), moving at most 25 percent from the last run's figure;
* every change is logged with the kind, old and new value, samples and reason, and kept for the Queue page (`profileclaim.log_change`).

A figure here never changes what a run generates: it only sets expectations and the judge's threshold, and `PLANETGEN_SELF_TUNING=0` turns it off.
"""

import json
import statistics
import time

from planetgen.queue import memorygate, profileclaim
from planetgen.util import log


def _key():
    return memorygate._key("time-profile")


COLUMNS = ("seconds", "seconds_per_item")
"""tuple: The figures kept per kind: the worked seconds of a unit, and the worked seconds an item of it takes (the judge's rate, which does not depend on the size of the unit)."""


def seconds_rows(kind, limit=profileclaim.WINDOW, conn=None, column="seconds"):
    """The newest `limit` finished rows of `kind` that hold `column`, as `seconds`, newest first, with the machine shape and epoch that measured them."""
    if column not in COLUMNS:
        raise ValueError(column)
    own = conn is None
    if own:
        from planetgen.db import store
        conn = store.get_control_connection()
    try:
        return list(conn.execute(f"SELECT {column} AS seconds, outcome, host_memory_mb, host_cores, generator_epoch FROM job_profiles WHERE kind = ? AND outcome = 'done' AND {column} IS NOT NULL"
                                 " ORDER BY id DESC LIMIT ?", (kind, int(limit))).fetchall())
    finally:
        if own:
            conn.close()


def target(rows, typed_seconds, host_memory, host_cores, epoch):
    """`(seconds or None, details)` from the rows of one kind, newest first (the rules in the module docstring)."""
    own = [r["seconds"] for r in rows if profileclaim._shape_of(r, host_memory, host_cores, epoch)]
    other = [r["seconds"] for r in rows if not profileclaim._shape_of(r, host_memory, host_cores, epoch)]
    if len(own) >= profileclaim.MIN_SAMPLES:
        return statistics.median(own), {"source": "own", "samples": len(own)}
    if len(other) >= profileclaim.MIN_SAMPLES:
        theirs = statistics.median(other)
        return (theirs if typed_seconds is None else 0.5 * theirs + 0.5 * typed_seconds), {"source": "other", "samples": len(other)}
    return None, {"source": "none", "samples": len(own) + len(other)}


def item_seconds(conn, command, token=None, read_rows=None, now=None):
    """The expected worked seconds of one item of a unit of `command` (the judge's expected rate is its inverse), by the rules of `unit_seconds`."""
    return unit_seconds(conn, command, token=token, read_rows=read_rows, now=now, column="seconds_per_item")


def unit_seconds(conn, command, token=None, typed_seconds=None, read_rows=None, now=None, column="seconds"):
    """
    The expected worked seconds of a unit of `command`, or `None` for 'use what you had'. `token` is the run's id (any other reader passes none and is held to the hour); `read_rows(kind)`
    gives the kind's rows (`seconds_rows` by default). Never raises.
    """
    if not profileclaim.enabled():
        return None
    kind = f"unit:{command}"
    try:
        now = time.time() if now is None else now
        token = str(token) if token is not None else f"hour-{int(now // 3600)}"
        field = kind if column == "seconds" else f"{kind}|{column}"
        stored = conn.hget(_key(), field)
        row = json.loads(memorygate._decode(stored)) if stored is not None else None
        if row is not None and row["token"] == token:
            return row["figure"]
        host = memorygate.read_host()
        from planetgen.util import generator_epoch
        wanted, details = target((read_rows or (lambda k: seconds_rows(k, column=column)))(kind), typed_seconds, host.total_mb, host.cores, generator_epoch.generator_epoch())
        if wanted is None:
            return None if row is None else row["figure"]
        previous = row["figure"] if row else typed_seconds
        figure = wanted if previous is None else min(max(wanted, previous * (1.0 - profileclaim.MAX_MOVE)), previous * (1.0 + profileclaim.MAX_MOVE))
        conn.hset(_key(), field, json.dumps({"token": token, "figure": figure}))
        if previous is None or abs(figure - previous) >= (0.05 if column == "seconds" else 1e-6):
            reason = f"median of {details['samples']} samples ({details['source']})" + ("" if abs(figure - wanted) <= 0.01 * wanted else f", limited to 25 percent a run from {previous:.4g} s")
            profileclaim.log_change(conn, "expected unit seconds" if column == "seconds" else "expected seconds per item", kind, previous, figure, reason, now=now, **details)
        return figure
    except Exception as exc:        # noqa: BLE001 -- an expectation is never worth a run
        log.debug(f"Self-tuning: the expected seconds of {kind} were not read from the profile ({type(exc).__name__}: {exc}).", area="queue")
        return None
