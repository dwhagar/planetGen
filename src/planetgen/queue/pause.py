# planetgen/queue/pause.py

"""
The pause of the work queues (PERF.152): one flag in Redis, `planetgen:queue:paused`, written by the admin Queue
page's Pause and Resume and read by every worker before it takes a job. A paused queue starts no new unit; a unit
already running finishes. The flag holds who paused and when, as JSON, so the page can say so. It replaces the
`paused` column of `work_lease` (control schema v21 dropped that table, PERF.100), which only the old executor ever read.

Pause is for the bulk and maintenance lanes (the units of a run, the cache warm-up, the database check). The
site's light jobs (`queue/lightjobs.py`) are never held: an admin who pauses the generation must still be able to edit and resume.

(OPS.235: the solo hold that the stuck-unit follower made for a restarted RQ unit, OPS.145, is gone with RQ: the launcher's hand-outs never read it.)

Without a Redis server nothing is paused (there is nothing to pause).
"""

import json
import time

from planetgen.queue import redisqueue, runrecord
from planetgen.util import log

def key():
    """The pause flag, `planetgen:queue:paused` (the tests' prefix, `runrecord.PREFIX`, in a test): a string whose value is
    `{"by": name, "at": Unix seconds}` as JSON, with no expiry."""
    return runrecord._k("queue", "paused")


NEVER_HELD = ("planetgen-web-",)
"""tuple: The queues (by name prefix) whose workers never wait for the flag: the Generate page's job queues (a paused queue must not hold the page's own starter)."""

WAIT_SECONDS = 2.0
"""float: How often a worker whose queue is paused looks at the flag again."""


def holds(queue_name):
    """Whether a worker of the queue `queue_name` waits while the queues are paused: every queue but `NEVER_HELD`
    (the bulk and maintenance lanes, and a run's own queue)."""
    return not queue_name.startswith(NEVER_HELD)


def state(connection=None):
    """
    The pause as it stands.

    Returns:
        dict: `paused` (bool) and `by` and `at` (who paused and when, `None` when not paused or not recorded).
    """
    free = {"paused": False, "by": None, "at": None}
    try:
        connection = connection or redisqueue.connect()
        raw = connection.get(key())
    except Exception as exc:  # noqa: BLE001 -- no Redis: nothing runs, nothing is paused
        log.swallowed("queue/pause", "state", exc)
        return free
    if raw is None:
        return free
    try:
        body = json.loads(raw)
        by, at = body.get("by"), body.get("at")
    except (ValueError, AttributeError):
        by = at = None
    return {"paused": True, "by": by, "at": at}


def is_paused(connection=None):
    """Whether the queues are paused (`state`)."""
    return state(connection)["paused"]


def pause(who, connection=None):
    """
    Pauses the queues; the first pause is the one kept.

    Returns:
        bool: Whether it was not paused already.
    """
    connection = connection or redisqueue.connect()
    body = {"by": str(who or "")[:64] or None, "at": time.time()}
    return bool(connection.set(key(), json.dumps(body), nx=True))


def resume(connection=None):
    """
    Lets the workers take jobs again.

    Returns:
        bool: Whether it was paused.
    """
    return bool((connection or redisqueue.connect()).delete(key()))
