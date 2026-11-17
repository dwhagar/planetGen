# planetgen/queue/waiting.py

"""
Why jobs wait (OPS.233): what the Queue page shows of the engine's ledger. For each job class the jobs that have not started, in two groups: `waiting` (a job it needs has not finished) and `ready` (it may
start and the supervisor's gate has not let it), the age of the oldest of each, how many jobs gave each reason for not starting (the words the gate wrote into `engine_jobs.waited`, `Launcher._note`), when
the next retry may start, and how long the jobs that started in the last hour waited from their making to their start. One read of the control database; it writes nothing and changes no job.

The site's own jobs (an API job, a slow page's answer, a visitor's prefetch) are light jobs with a record in Redis and no row in the ledger (OPS.230, `queue/lightjobs.py`); `light_overview` lists the ones that
have not started, by lane (UX.320). The page reads both through `GET /api/admin/work/waiting`. What the controller last did is the supervisor panel's list of acts (OPS.214), which the memory controller of OPS.231 fills.
"""

import collections
import time

from planetgen.queue import lightjobs, memorygate
from planetgen.util import log

STARTED_WINDOW_SECONDS = 3600
"""int: How far back the starts count that `started` reports, in seconds."""

REASONS_PER_CLASS = 5
"""int: How many different reasons a class lists; the rest are counted under `more`."""


def _seconds(value):
    """A number of seconds from the database's figure (a `Decimal`, an int or None), 0.0 for none, never below 0."""
    return 0.0 if value is None else max(float(value), 0.0)


def overview(conn, with_gate=True):
    """
    The jobs that have not started, by class.

    Args:
        conn: A connection to the control database.
        with_gate (bool): Whether to add the gate's figures (`gate`, `memorygate.gate_words`); a failure to read them leaves `gate` empty and the rest whole.

    Returns:
        dict: `classes` (one dict for each class that has a job not started: `job_class`, `waiting`, `ready`, `oldest_waiting_seconds`, `oldest_ready_seconds`, `retry_in_seconds` or None, `reasons`
            as `{"reason", "count"}` biggest first and `more` the count of jobs under the reasons left out), `started` (for each class a job started from in the last hour: `job_class`, `count`,
            `median_seconds`, `longest_seconds` from the job being made to its first start) and `gate` (str).
    """
    classes = {}
    for row in conn.execute(
            "SELECT job_class, state, COUNT(*) AS n, TIMESTAMPDIFF(MICROSECOND, MIN(created_at), NOW(6)) / 1000000 AS oldest_made,"
            " TIMESTAMPDIFF(MICROSECOND, MIN(updated_at), NOW(6)) / 1000000 AS oldest_ready,"
            " TIMESTAMPDIFF(MICROSECOND, NOW(6), MIN(CASE WHEN not_before > NOW(6) THEN not_before END)) / 1000000 AS retry_in"
            " FROM engine_jobs WHERE state IN ('waiting', 'ready') GROUP BY job_class, state").fetchall():
        item = classes.setdefault(row["job_class"], {"job_class": row["job_class"], "waiting": 0, "ready": 0, "oldest_waiting_seconds": 0.0, "oldest_ready_seconds": 0.0,
                                                     "retry_in_seconds": None, "reasons": [], "more": 0})
        item[row["state"]] = int(row["n"])
        if row["state"] == "waiting":
            item["oldest_waiting_seconds"] = _seconds(row["oldest_made"])
        else:
            item["oldest_ready_seconds"] = _seconds(row["oldest_ready"])
        if row["retry_in"] is not None:
            soon = _seconds(row["retry_in"])
            item["retry_in_seconds"] = soon if item["retry_in_seconds"] is None else min(item["retry_in_seconds"], soon)
    for row in conn.execute("SELECT job_class, waited, COUNT(*) AS n FROM engine_jobs WHERE state = 'ready' AND waited IS NOT NULL GROUP BY job_class, waited ORDER BY n DESC, waited").fetchall():
        item = classes.get(row["job_class"])
        if item is None:
            continue
        if len(item["reasons"]) < REASONS_PER_CLASS:
            item["reasons"].append({"reason": row["waited"], "count": int(row["n"])})
        else:
            item["more"] += int(row["n"])
    started = {}
    for row in conn.execute(
            "SELECT j.job_class, TIMESTAMPDIFF(MICROSECOND, j.created_at, a.started_at) / 1000000 AS made_to_start FROM engine_attempts a JOIN engine_jobs j ON j.id = a.job_id"
            " WHERE a.attempt = 1 AND a.started_at >= DATE_SUB(NOW(6), INTERVAL ? SECOND) ORDER BY j.job_class, made_to_start", (STARTED_WINDOW_SECONDS,)).fetchall():
        started.setdefault(row["job_class"], []).append(_seconds(row["made_to_start"]))
    gate = ""
    if with_gate:
        try:
            gate = memorygate.gate_words(conn)
        except Exception as exc:  # noqa: BLE001 -- the page is for reading; a gate that cannot be read leaves its line out
            log.swallowed("queue/waiting", "gate_words", exc)
    return {"classes": sorted(classes.values(), key=lambda item: item["job_class"]),
            "started": [{"job_class": name, "count": len(times), "median_seconds": times[len(times) // 2], "longest_seconds": times[-1]} for name, times in sorted(started.items())],
            "gate": gate}


LIGHT_WINDOW = 500
"""int: How many of the oldest pending light jobs `light_overview` reads (a bound on the work of a poll; `total` says how many wait)."""

LIGHT_FUNCTIONS_PER_LANE = 4
"""int: How many different functions a lane lists, the most common first."""

NO_REASON = "no reason yet: the supervisor has not looked at it"
"""str: What a pending light job that the dispatch rule has not judged yet says."""


def light_overview(connection, now=None):
    """
    The light jobs that have not started (UX.320), by lane. Reads the pending list and each record once; writes nothing.

    Args:
        connection (redis.Redis): The server the supervisor reads.
        now (float, optional): The time now in Unix seconds (a test gives it).

    Returns:
        dict: `lanes` (one dict for each lane with a job waiting: `lane`, `count`, `oldest_seconds`, `deadline_seconds` the least time any job of the lane has left before it is dropped for want of anyone
            asking (negative: it is due, `None`: no job of the lane has a deadline), `reasons` as `{"reason", "count"}` biggest first, `functions` as `{"function", "count"}` biggest first), `total`
            (how many jobs wait) and `shown` (how many of them were read).
    """
    now = time.time() if now is None else now
    ids = [lightjobs._text(item) for item in connection.zrange(lightjobs.pending_key(), 0, LIGHT_WINDOW - 1)]
    total = int(connection.zcard(lightjobs.pending_key()))
    pipe = connection.pipeline(transaction=False)
    for job_id in ids:
        pipe.hmget(lightjobs.job_key(job_id), "state", "lane", "function", "created", "wanted_until", "waited")
    lanes = {}
    for found in pipe.execute():
        state, lane, function, created, wanted_until, waited = [None if value is None else lightjobs._text(value) for value in found]
        if state != "queued":
            continue
        item = lanes.setdefault(lane or lightjobs.LANE_INTERACTIVE, {"lane": lane or lightjobs.LANE_INTERACTIVE, "count": 0, "oldest_seconds": 0.0, "deadline_seconds": None,
                                                                     "reasons": collections.Counter(), "functions": collections.Counter()})
        item["count"] += 1
        if created:
            item["oldest_seconds"] = max(item["oldest_seconds"], now - float(created))
        if wanted_until:
            left = float(wanted_until) - now
            item["deadline_seconds"] = left if item["deadline_seconds"] is None else min(item["deadline_seconds"], left)
        item["reasons"][waited or NO_REASON] += 1
        item["functions"][(function or "?").rpartition(":")[2]] += 1
    for item in lanes.values():
        item["reasons"] = [{"reason": reason, "count": count} for reason, count in item["reasons"].most_common(REASONS_PER_CLASS)]
        item["functions"] = [{"function": name, "count": count} for name, count in item["functions"].most_common(LIGHT_FUNCTIONS_PER_LANE)]
    return {"lanes": sorted(lanes.values(), key=lambda item: item["lane"]), "total": total, "shown": len(ids)}
