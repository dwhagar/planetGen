# planetgen/queue/jobprofiles.py

"""
What each kind of job costs in memory and time, kept as samples in the control database (OPS.155, the first piece of the self-tuning jobs in
`docs/design/operations/self-tuning-jobs.md`).

**Record only.** Nothing reads these rows but the readers named below, so no run changes until they are built: the gate's claim, the worker counts, the splits and the ETA still come from the typed tables and
the Redis hashes. The rows are there so the data gathers while the readers (OPS.156 to OPS.160) are built.

* **Peak and average.** Memory, time and processor each have both. A unit's come from a `usage.Usage` the unit's beat thread feeds while it runs; a gated job's from the gate's heartbeat; a stage's
  from its units (the longest and the mean unit, the largest peak, the time-weighted mean of the means).
* **One row per finished unit or job**, in `job_profiles` (`control_schema.sql`). The writers are the paths that already measured it (REUSE RULE): a unit
  (`units._do_unit`, with the items and worked seconds `unitrates.note` takes and the peak `units.peak_resident_mb` takes), any other job at the gate's release
  (`memorygate.release`: its peak and the seconds it held its claim), a stage (`stagetimes.record`: its wall seconds) and a memory kill (`units.horse_killed`).
  A unit's row is written by the unit, not by the gate's release, so a unit is not counted twice.
* **Per unit and per item.** A row has the items the unit did, the seconds it worked and its peak MB, and from those the seconds per item and the MB per item. MB per
  item is the raw peak over the items, kept beside the size bucket (a power of two of the items) so a reader can fit `peak = base + slope x items` from rows of two
  or more buckets: a unit holds a fixed base before its first item, so the raw figure alone overstates the slope.
* **The machine.** Each row carries the host's memory and cores, the MySQL buffer pool and the generator epoch it was measured under, so a reader can trust a row
  from another machine shape or epoch at half weight.
* **Bounded.** A kind and bucket keep their newest `KEEP_PER_BUCKET` rows (`prune`, run now and then by the writer that crosses the limit).
* **Never blocks a job.** Every call fails open: on any error the row is dropped, and after a failure the writer leaves the database alone for `BACKOFF_SECONDS`.
  A reading that is missing or implausible (no seconds, a peak below `memorygate.MIN_CLAIM_MB`) is stored as empty, and a row with nothing left is not stored.
  `PLANETGEN_JOB_PROFILES=0` turns the writers off.
"""

import math
import os
import time

from planetgen.util import log

ENV_VAR = "PLANETGEN_JOB_PROFILES"
"""str: Set to `0`, `off`, `no` or `false` and nothing is recorded (tests set it; `tests/test_job_profiles.py` turns it on)."""

KEEP_PER_BUCKET = 200
"""int: The rows a kind keeps in one size bucket; the oldest beyond twice this go at the next prune."""

BACKOFF_SECONDS = 60.0
"""float: After a failed write, the writer does not try again for this long."""

OUTCOMES = ("done", "killed_memory", "timed_out", "stuck")
"""tuple: How a unit or job ended (the column's CHECK)."""

MIN_SECONDS = 0.05
"""float: A job that worked for less than this says nothing about its time."""

_state = {"retry_after": 0.0, "machine": None}


def enabled():
    return os.environ.get(ENV_VAR, "1").strip().lower() not in ("0", "off", "no", "false")


def size_bucket(items):
    """The power of two of `items` (`floor(log2(items))`), or -1 for a job that has no items."""
    if items is None or items < 1:
        return -1
    return int(math.floor(math.log2(items)))


def _machine(conn):
    """`(memory MB, cores, MySQL buffer pool MB or None)` of this machine, worked out once per process."""
    if _state["machine"] is None:
        from planetgen.queue import memorygate
        host = memorygate.read_host()
        pool = None
        try:
            pool = int(conn.execute("SELECT @@innodb_buffer_pool_size AS b").fetchone()["b"] / (1024 * 1024))
        except Exception:       # noqa: BLE001 -- an unknown pool size is stored empty
            pass
        _state["machine"] = (int(host.total_mb), int(host.cores), pool)
    return _state["machine"]


def _plain(value, floor=0.0):
    return float(value) if value is not None and value >= floor and value > 0 else None


def sample(kind, items=None, seconds=None, peak_mb=None, outcome="done", avg_mb=None, cpu_peak=None, cpu_avg=None, seconds_peak=None, seconds_avg=None, parts=1):
    """
    The row `record` would store for these readings, as a dict, or `None` when it holds nothing useful. Implausible readings are emptied: seconds under
    `MIN_SECONDS`, a memory under `memorygate.MIN_CLAIM_MB`, items under 1 and a processor reading of zero. A row for a job that did not end `done` is kept even when
    it has no reading, because the kill is the data.
    """
    from planetgen.queue import memorygate
    items = float(items) if items is not None and items >= 1 else None
    seconds = _plain(seconds, MIN_SECONDS)
    peak_mb, avg_mb = _plain(peak_mb, memorygate.MIN_CLAIM_MB), _plain(avg_mb, memorygate.MIN_CLAIM_MB)
    cpu_peak, cpu_avg = _plain(cpu_peak), _plain(cpu_avg)
    seconds_peak, seconds_avg = _plain(seconds_peak, MIN_SECONDS), _plain(seconds_avg, MIN_SECONDS)
    if outcome not in OUTCOMES:
        return None
    if outcome == "done" and all(value is None for value in (seconds, peak_mb, avg_mb, cpu_peak, cpu_avg)):
        return None
    return {"kind": kind, "size_bucket": size_bucket(items), "items": items, "parts": max(int(parts or 1), 1), "seconds": seconds, "seconds_peak": seconds_peak,
            "seconds_avg": seconds_avg, "peak_mb": peak_mb, "avg_mb": avg_mb, "cpu_peak": cpu_peak, "cpu_avg": cpu_avg, "outcome": outcome,
            "seconds_per_item": seconds / items if seconds and items else None, "mb_per_item": peak_mb / items if peak_mb and items else None}


_COLUMNS = ("kind", "size_bucket", "items", "parts", "seconds", "seconds_peak", "seconds_avg", "peak_mb", "avg_mb", "cpu_peak", "cpu_avg", "seconds_per_item",
            "mb_per_item", "outcome")


def _said(row):
    """The sample as one line of fields for the log: the kind, then every reading it has."""
    figures = " ".join(f"{name}={row[name]:.4g}" if isinstance(row[name], float) else f"{name}={row[name]}"
                       for name in _COLUMNS[1:] if row[name] is not None and name != "size_bucket")
    return f"{row['kind']} bucket={row['size_bucket']} {figures}"


def record(kind, items=None, seconds=None, peak_mb=None, outcome="done", conn=None, **more):
    """
    Stores one sample for the job kind `kind` (the gate's strings: `unit:<command>`, `cut:<command>`, `api:<ref>`, `fn:<name>`, and `stage:<key>` for a stage). `more`
    is any of `avg_mb`, `cpu_peak`, `cpu_avg`, `seconds_peak`, `seconds_avg` and `parts` (see `sample`). `conn` is an open control connection (a test gives one);
    otherwise one is opened from the pool and closed. Every sample stored is logged (Boss, 3:03 PM EST 2026-11-12). Returns whether a row was stored; never raises.
    """
    if not enabled() or time.monotonic() < _state["retry_after"]:
        return False
    row = sample(kind, items, seconds, peak_mb, outcome, **more)
    if row is None:
        return False
    own = conn is None
    try:
        if own:
            from planetgen.db import store
            conn = store.get_control_connection()
        memory, cores, pool = _machine(conn)
        from planetgen.util import generator_epoch
        cursor = conn.execute(
            f"INSERT INTO job_profiles ({', '.join(_COLUMNS)}, host_memory_mb, host_cores, buffer_pool_mb, generator_epoch, recorded_at)"
            f" VALUES ({', '.join('?' * (len(_COLUMNS) + 4))}, NOW(6))",
            tuple(row[name] for name in _COLUMNS) + (memory, cores, pool, generator_epoch.generator_epoch()))
        conn.commit()
        log.debug(f"Job profile sample: {_said(row)} host={memory}MB/{cores}cores.", area="queue")
        if cursor.lastrowid % 25 == 0:
            prune(conn, row["kind"], row["size_bucket"])
        return True
    except Exception as exc:        # noqa: BLE001 -- a profile is never worth a job
        _state["retry_after"] = time.monotonic() + BACKOFF_SECONDS
        log.debug(f"Job profiles: {kind} was not recorded ({exc}).", area="queue")
        return False
    finally:
        if own and conn is not None:
            try:
                conn.close()
            except Exception:       # noqa: BLE001
                pass


def peaks(kind, limit=30, conn=None):
    """
    The newest `limit` rows of `kind` that hold a peak memory, newest first, as dicts (`peak_mb`, `outcome` and the machine shape and epoch that measured them): the reading of
    `profileclaim` (OPS.156). Raises on a database that cannot answer, and the caller falls back.
    """
    own = conn is None
    if own:
        from planetgen.db import store
        conn = store.get_control_connection()
    try:
        return list(conn.execute("SELECT peak_mb, outcome, host_memory_mb, host_cores, generator_epoch FROM job_profiles WHERE kind = ? AND peak_mb IS NOT NULL"
                                 " ORDER BY id DESC LIMIT ?", (kind, int(limit))).fetchall())
    finally:
        if own:
            conn.close()


def unit_rows(kind, limit=60, conn=None):
    """
    The newest `limit` finished rows of `kind` that hold seconds per item (newest first), as dicts with the items, size bucket, peak MB and the machine shape and epoch that measured
    them: the reading of `splitprofile` (OPS.161). Raises on a database that cannot answer, and the caller falls back.
    """
    own = conn is None
    if own:
        from planetgen.db import store
        conn = store.get_control_connection()
    try:
        return list(conn.execute("SELECT items, size_bucket, seconds_per_item, peak_mb, host_memory_mb, host_cores, generator_epoch FROM job_profiles"
                                 " WHERE kind = ? AND outcome = 'done' AND seconds_per_item IS NOT NULL AND items IS NOT NULL ORDER BY id DESC LIMIT ?",
                                 (kind, int(limit))).fetchall())
    finally:
        if own:
            conn.close()


def stage_figures(entries):
    """
    The figures of a stage from what its units reported (`{"seconds", "wall", "peak_mb", "avg_mb", "cpu_peak", "cpu_seconds"}` each): the units counted, the longest and mean unit,
    the largest memory peak, the mean memory weighted by each unit's seconds, the busiest unit's processor and the processor seconds over the units' seconds.
    """
    entries = [e for e in entries if e]
    if not entries:
        return {"parts": 1}
    worked = [e["seconds"] for e in entries if e.get("seconds")]
    weighted = [(e["avg_mb"], e.get("wall") or e.get("seconds") or 0.0) for e in entries if e.get("avg_mb")]
    walls = sum(e.get("wall") or 0.0 for e in entries)
    return {"parts": len(entries), "seconds_peak": max(worked) if worked else None, "seconds_avg": sum(worked) / len(worked) if worked else None,
            "peak_mb": max((e["peak_mb"] for e in entries if e.get("peak_mb")), default=None),
            "avg_mb": sum(a * w for a, w in weighted) / sum(w for _, w in weighted) if weighted and sum(w for _, w in weighted) > 0 else None,
            "cpu_peak": max((e["cpu_peak"] for e in entries if e.get("cpu_peak")), default=None),
            "cpu_avg": sum(e.get("cpu_seconds") or 0.0 for e in entries) / walls if walls > 0 else None}


def prune(conn, kind, bucket, keep=KEEP_PER_BUCKET):
    """Keeps the newest `keep` rows of the kind and bucket once there are more than twice that; returns the rows deleted."""
    count = conn.execute("SELECT COUNT(*) AS n FROM job_profiles WHERE kind = ? AND size_bucket = ?", (kind, bucket)).fetchone()["n"]
    if count <= 2 * keep:
        return 0
    edge = conn.execute("SELECT id FROM job_profiles WHERE kind = ? AND size_bucket = ? ORDER BY id DESC LIMIT 1 OFFSET ?",
                        (kind, bucket, keep - 1)).fetchone()
    if edge is None:
        return 0
    removed = conn.execute("DELETE FROM job_profiles WHERE kind = ? AND size_bucket = ? AND id < ?", (kind, bucket, edge["id"])).rowcount
    conn.commit()
    return removed
