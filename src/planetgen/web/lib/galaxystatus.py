# planetgen/web/lib/galaxystatus.py

"""
Whether a background run is changing the galaxy, for the banner every visitor sees (UX.3).

The run that owns the machine lives in Redis (`queue/runrecord.py`: the `planetgen:run:active` key names it, its
record holds the stages and counters), whether it was started from the Generate page or the command line, so one
read answers for both. `runeta.public_view` is the public part (no user, no arguments); `hourdisplay.hours_left`
turns the time left into the one whole-hour figure Boss asked for. The figure is remembered in Redis
(`runrecord.banner_key`) so a visitor on one web process sees what a visitor on another does and the figure only
rises by the rule of `hourdisplay`.

A process keeps the answer for `CACHE_SECONDS`, so a page render and the 30 second poll of every open tab cost a
dictionary lookup. Every failure fails open: with Redis unreachable, or no run, or a run that shows no sign of
life (`runrecord.is_stale`: the Queue page offers to clear it), there is no banner.

A run that names another database than the site's is not shown: it is not changing what this site shows.
"""

import threading
import time

import redis

from planetgen.queue import hourdisplay, redisclient, runeta, runrecord

CACHE_SECONDS = 5.0
"""float: How long a process keeps the answer."""

SOCKET_SECONDS = 0.5
"""float: How long one Redis call may take before the page simply has no banner."""

BANNER_KEEP_SECONDS = 86400
"""int: How long the figure shown is remembered."""

INACTIVE = {"active": False}

_lock = threading.Lock()
_cache = {}


def _client():
    return redisclient.client(timeout=SOCKET_SECONDS, retries=1)


def clear_cache():
    """Forgets every kept answer (tests, and a run that has just changed state)."""
    with _lock:
        _cache.clear()


def message(hours, view):
    """The banner's sentence: that the galaxy is being changed and, when it can be said, when it should finish."""
    lead = "The galaxy is being changed by a background job, so some places may be missing or change while you look."
    if view.get("stalled"):
        return f"{lead} The job is waiting on progress; the finish time returns when it resumes."
    if hours is None:
        return f"{lead} Working out how long it will take."
    unit = "hour" if hours == 1 else "hours"
    guess = " (a guess from earlier runs)" if view.get("from_earlier_runs") else ""
    return f"{lead} It should finish in about {hours} {unit}{guess}."


def _remember(conn, run_id, eta_seconds):
    """The hours to show for `eta_seconds`, kept for the next visitor."""
    key = runrecord.banner_key(run_id)
    raw = conn.get(key)
    shown = int(raw) if raw else None
    hours = hourdisplay.hours_left(eta_seconds, shown)
    if hours != shown:
        conn.set(key, hours, ex=BANNER_KEEP_SECONDS)
    return hours


def compute(conn, database, now=None):
    """
    The public status of the run in progress, read from Redis `conn` (no cache).

    Returns:
        dict: `{"active": False}`, or `active`, `hours` (whole hours, `None` while there is no estimate),
        `percent` (or `None`), `stage` (label), and `text`, the banner's sentence.
    """
    now = time.time() if now is None else now
    run_id = runrecord.active_run(conn)
    run = runrecord.read_run(conn, run_id) if run_id else None
    if run is None or runrecord.is_stale(run, now):
        return dict(INACTIVE)
    if run["database"] and database and run["database"] != database:
        return dict(INACTIVE)
    view = runeta.public_view(run, now)
    if not view["active"]:
        return dict(INACTIVE)
    seconds = view["eta_seconds"]
    hours = None if seconds is None or view["stalled"] else _remember(conn, run_id, seconds)
    return {"active": True, "hours": hours, "percent": view["percent"], "stage": view["stage"]["label"],
            "text": message(hours, view)}


def status(database):
    """The same, kept for `CACHE_SECONDS`; no banner when Redis can't be read."""
    now = time.monotonic()
    with _lock:
        kept = _cache.get(database)
        if kept is not None and kept[0] > now:
            return dict(kept[1])
    try:
        answer = compute(_client(), database)
    except (redis.RedisError, OSError, ValueError, KeyError):
        answer = dict(INACTIVE)
    with _lock:
        _cache[database] = (now + CACHE_SECONDS, answer)
    return dict(answer)
