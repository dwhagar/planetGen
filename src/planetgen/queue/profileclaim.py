# planetgen/queue/profileclaim.py

"""
The memory claim worked out from the job profile (OPS.156, step 2 of `docs/design/operations/self-tuning-jobs.md`).

A kind with enough history on this machine claims the 95th percentile of its recent peaks (`jobprofiles` rows, newest `WINDOW`) times a margin that falls from `MARGIN_AT_MIN` (1.25, at
`MIN_SAMPLES` samples) to `MARGIN_AT_FULL` (1.10, at `WINDOW`). The typed figure is then no floor, so a constant that was too generous on this machine is corrected downward, and the largest
peak ever seen (`memorygate.learned`) no longer holds a claim up for good.

In order, per kind:

1. **Own history**: at least `MIN_SAMPLES` finished samples from this machine shape (memory within `SHAPE_TOLERANCE`, same cores) and the current generator epoch: the figure above.
2. **Other history**: at least `MIN_SAMPLES` samples from another machine shape or epoch count at half weight: the mean of their figure (at the full `MARGIN_AT_MIN`) and the typed figure (their figure alone where there is no typed figure).
3. **No history**: nothing here, and `memorygate.claim_mb` gives what it always gave. A fresh install, a new kind, a wiped control database or `PLANETGEN_SELF_TUNING=0` therefore behave as before.

**A memory kill wins.** A kill of the kind in the window (`jobprofiles` `killed_memory`, its row holds the size the retry was given less the margin) sets a floor at that size times `LEARN_MARGIN`,
which the slow-move rule below does not hold back.

**Slow and frozen.** The figure is worked out once per run (the run's id is the token; any other job uses the hour) and kept in Redis: the jobs of the run claim that figure, and a new run moves it by at most
`MAX_MOVE` (25 percent) from the last. The first move away from a typed figure is limited the same way.

**Logged** (Boss, 3:03 PM EST 2026-11-12): every change of a figure logs the kind, old and new, the samples behind it, the source and the reason, and is kept in the Redis list `changes` (newest `KEEP_CHANGES`) for the Queue page.
**Never blocks a job**: any failure gives `None`, the typed path.
"""

import json
import math
import os
import time

from planetgen.queue import memorygate
from planetgen.util import log

ENV_VAR = "PLANETGEN_SELF_TUNING"
"""str: Set to `0`, `off`, `no` or `false` and the profile is not read: every claim is what it was before OPS.156."""

WINDOW = 30
MIN_SAMPLES = 3
MARGIN_AT_MIN = 1.25
MARGIN_AT_FULL = 1.10
MAX_MOVE = 0.25
SHAPE_TOLERANCE = 0.10
KEEP_CHANGES = 100


def enabled():
    return os.environ.get(ENV_VAR, "1").strip().lower() not in ("0", "off", "no", "false")


def margin(samples):
    """1.25 at 3 samples, falling in a line to 1.10 at 30 or more."""
    share = min(max((samples - MIN_SAMPLES) / (WINDOW - MIN_SAMPLES), 0.0), 1.0)
    return MARGIN_AT_MIN + (MARGIN_AT_FULL - MARGIN_AT_MIN) * share


def percentile(values, share=0.95):
    """The `share` quantile of `values` by linear interpolation between the nearest ranks."""
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    at = share * (len(ordered) - 1)
    low = int(math.floor(at))
    return ordered[low] + (ordered[min(low + 1, len(ordered) - 1)] - ordered[low]) * (at - low)


def _shape_of(row, host_memory, host_cores, epoch):
    near = abs(row["host_memory_mb"] - host_memory) <= SHAPE_TOLERANCE * host_memory
    return near and row["host_cores"] == host_cores and row["generator_epoch"] == epoch


def target(rows, typed_mb, host_memory, host_cores, epoch):
    """
    `(figure MB or None, kill_floor MB or None, details)` from the profile rows of one kind (`jobprofiles.peaks`), newest first. `figure` is `None` when there is no usable history. `details`
    names the source (`own`, `other`) and the samples behind the figure.
    """
    own = [r for r in rows if r["outcome"] == "done" and _shape_of(r, host_memory, host_cores, epoch)]
    other = [r for r in rows if r["outcome"] == "done" and not _shape_of(r, host_memory, host_cores, epoch)]
    kills = [r["peak_mb"] for r in rows if r["outcome"] == "killed_memory" and _shape_of(r, host_memory, host_cores, epoch)]
    floor = max(kills) * memorygate.LEARN_MARGIN if kills else None
    if len(own) >= MIN_SAMPLES:
        peaks = [r["peak_mb"] for r in own]
        k = margin(len(peaks))
        return percentile(peaks) * k, floor, {"source": "own", "samples": len(peaks), "p95_mb": round(percentile(peaks)), "margin": round(k, 3)}
    if len(other) >= MIN_SAMPLES:
        peaks = [r["peak_mb"] for r in other]
        theirs = percentile(peaks) * MARGIN_AT_MIN
        figure = theirs if typed_mb is None else 0.5 * theirs + 0.5 * typed_mb
        return figure, floor, {"source": "other", "samples": len(peaks), "p95_mb": round(percentile(peaks)), "margin": MARGIN_AT_MIN}
    return None, floor, {"source": "none", "samples": len(own) + len(other)}


def settle(previous, wanted, floor):
    """The figure for a new run: `wanted` moved from `previous` by at most `MAX_MOVE`, and never under the kill `floor`."""
    figure = wanted if previous is None else min(max(wanted, previous * (1.0 - MAX_MOVE)), previous * (1.0 + MAX_MOVE))
    return max(figure, floor) if floor else figure


def _changes_key():
    return memorygate._key("changes")


def claim(conn, kind, typed_mb, token=None, read_rows=None, now=None):
    """
    The claim of `kind` from the profile, MB, or `None` for 'use the typed path'. `token` is the run's id (any other job passes none and is held to the hour). `read_rows(kind)` gives the
    kind's profile rows (`jobprofiles.peaks` by default; a test gives its own). Never raises.
    """
    pinned = pinned_mb(conn, kind)
    if pinned is not None:
        return pinned           # OPS.158: an admin's figure wins, and stays until it is cleared
    if not enabled():
        return None
    try:
        now = time.time() if now is None else now
        token = str(token) if token is not None else f"hour-{int(now // 3600)}"
        stored = conn.hget(memorygate._key("profile"), kind)
        row = json.loads(memorygate._decode(stored)) if stored is not None else None
        if row is not None and row["token"] == token:
            return row["figure"]
        from planetgen.queue import jobprofiles
        rows = (read_rows or jobprofiles.peaks)(kind)
        host = memorygate.read_host()
        from planetgen.util import generator_epoch
        wanted, floor, details = target(rows, typed_mb, host.total_mb, host.cores, generator_epoch.generator_epoch())
        if wanted is None:
            return None if floor is None else max(floor, row["figure"] if row else 0.0, typed_mb or 0.0)
        previous = row["figure"] if row else typed_mb
        figure = settle(previous, wanted, floor)
        conn.hset(memorygate._key("profile"), kind, json.dumps({"token": token, "figure": figure}))
        if previous is None or abs(figure - previous) >= 1.0:
            _log_change(conn, kind, previous, figure, wanted, floor, details, typed_mb, now)
        return figure
    except Exception as exc:        # noqa: BLE001 -- the profile is never worth a job: the typed claim stands
        log.debug(f"Self-tuning: the claim of {kind} was not read from the profile ({type(exc).__name__}: {exc}).", area="queue")
        return None


def _log_change(conn, kind, old, new, wanted, floor, details, typed_mb, now):
    held = "a memory kill floor" if floor and new == floor else ("limited to 25 percent a run" if abs(new - wanted) >= 1.0 else "the profile's figure")
    reason = f"{held}; source {details['source']}, {details['samples']} samples" + (f", p95 {details['p95_mb']} MB x margin {details['margin']}" if "p95_mb" in details else "")
    log.normal(f"Self-tuning: claim of {kind} {'was ' + format(old, ',.0f') + ' MB' if old is not None else 'had no earlier figure'}, now {new:,.0f} MB (typed {typed_mb if typed_mb is not None else 'none'}; {reason}).")
    conn.lpush(_changes_key(), json.dumps({"at": now, "what": "claim", "kind": kind, "old": old, "new": new, "typed": typed_mb, "reason": reason, **details}))
    conn.ltrim(_changes_key(), 0, KEEP_CHANGES - 1)


def _pins_key():
    return memorygate._key("pins")


def pinned_mb(conn, kind):
    """The figure an admin fixed for `kind`, MB, or `None` (OPS.158)."""
    try:
        raw = conn.hget(_pins_key(), kind)
        return None if raw is None else float(memorygate._decode(raw))
    except Exception as exc:        # noqa: BLE001 -- the pins are not worth a job
        log.swallowed("queue", f"reading the pin of {kind}", exc)
        return None


def pin(conn, kind, mb, who="an admin"):
    """Fixes the claim of `kind` at `mb` until `unpin` or `clear`; logged. Returns the figure."""
    old = pinned_mb(conn, kind)
    conn.hset(_pins_key(), kind, f"{float(mb):.3f}")
    log_change(conn, "pin", kind, old, float(mb), f"fixed by {who}")
    return float(mb)


def unpin(conn, kind, who="an admin"):
    """Lets the profile decide the claim of `kind` again; logged. Returns whether a pin was there."""
    old = pinned_mb(conn, kind)
    if old is None:
        return False
    conn.hdel(_pins_key(), kind)
    log_change(conn, "unpin", kind, old, 0.0, f"released by {who}: the profile decides again")
    return True


def clear(conn, kind, who="an admin"):
    """
    Forgets everything learned about `kind`: the settled claim, the largest peak seen (`memorygate.learned`), the stored starting limit and its profile rows in the control database
    (a pin stays, `unpin` lifts it). The next run claims from the typed table and measures again. Logged. Returns the profile rows deleted.
    """
    old = None
    stored = conn.hget(memorygate._key("profile"), kind)
    if stored is not None:
        old = json.loads(memorygate._decode(stored))["figure"]
    conn.hdel(memorygate._key("profile"), kind)
    conn.hdel(memorygate._key("learned"), kind)
    conn.hdel(memorygate._key("adaptive:settled"), kind)
    removed = 0
    try:
        from planetgen.db import store
        db = store.get_control_connection()
        try:
            removed = db.execute("DELETE FROM job_profiles WHERE kind = ?", (kind,)).rowcount
            db.commit()
        finally:
            db.close()
    except Exception as exc:        # noqa: BLE001 -- the Redis side is cleared whatever the database says
        log.debug(f"Self-tuning: the profile rows of {kind} were not deleted ({exc}).", area="queue")
    log_change(conn, "clear", kind, old, 0.0, f"forgotten by {who}: {removed} profile rows deleted, the next run claims from the typed table")
    return removed


def overview(conn):
    """
    What the page shows (OPS.158): `{"self_tuning", "kinds", "changes"}`. A kind is every one that has a settled figure, a pin or profile rows: `kind`, `typed_mb`, `figure_mb` (what a run
    claims from the profile now), `pinned_mb`, `samples`, `last_sample` (UTC text) and `changed` (the newest change made to it). The kinds come out sorted by name.
    """
    from planetgen.queue import claims
    kinds = {}
    for name, raw in conn.hgetall(memorygate._key("profile")).items():
        kinds.setdefault(memorygate._decode(name), {})["figure_mb"] = json.loads(memorygate._decode(raw))["figure"]
    for name, raw in conn.hgetall(_pins_key()).items():
        kinds.setdefault(memorygate._decode(name), {})["pinned_mb"] = float(memorygate._decode(raw))
    try:
        from planetgen.db import store
        db = store.get_control_connection()
        try:
            for row in db.execute("SELECT kind, COUNT(*) AS n, MAX(recorded_at) AS last FROM job_profiles GROUP BY kind").fetchall():
                entry = kinds.setdefault(row["kind"], {})
                entry["samples"], entry["last_sample"] = int(row["n"]), str(row["last"])[:19] if row["last"] else None
        finally:
            db.close()
    except Exception as exc:        # noqa: BLE001 -- a page without the control database still shows the Redis side
        log.debug(f"Self-tuning: the profile rows were not counted ({exc}).", area="queue")
    recent = changes(conn, KEEP_CHANGES)
    rows = []
    for name in sorted(kinds):
        entry = kinds[name]
        newest = next((c for c in recent if c["kind"] == name), None)
        rows.append({"kind": name, "typed_mb": claims.typed_mb(name), "figure_mb": entry.get("figure_mb"), "pinned_mb": entry.get("pinned_mb"),
                     "samples": entry.get("samples", 0), "last_sample": entry.get("last_sample"), "changed": newest})
    return {"self_tuning": enabled(), "kinds": rows, "changes": recent[:20]}


def log_change(conn, what, kind, old, new, reason, now=None, **extra):
    """
    Logs a change the engine made (`what`: `claim`, `workers` or `start limit`; later `split` and `time`) with the kind, the old and new value and the reason, and keeps it in the list `changes`
    reads (Boss, 3:03 PM EST 2026-11-12). Never raises.
    """
    try:
        now = time.time() if now is None else now
        log.normal(f"Self-tuning: {what} of {kind} {'was ' + format(old, ',.0f') if old is not None else 'had no earlier figure'}, now {new:,.0f} ({reason}).")
        conn.lpush(_changes_key(), json.dumps({"at": now, "what": what, "kind": kind, "old": old, "new": new, "reason": reason, **extra}))
        conn.ltrim(_changes_key(), 0, KEEP_CHANGES - 1)
    except Exception as exc:        # noqa: BLE001
        log.debug(f"Self-tuning: a change of {kind} was not logged ({exc}).", area="queue")


def changes(conn, limit=20):
    """The newest changes the engine made, newest first (for the Queue page)."""
    return [json.loads(memorygate._decode(item)) for item in conn.lrange(_changes_key(), 0, limit - 1)]
