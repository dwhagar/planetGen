# planetgen/queue/pause.py

"""
The pause of the work queues (PERF.152): one flag in Redis, `planetgen:queue:paused`, written by the admin Queue
page's Pause and Resume and read by every worker before it takes a job. A paused queue starts no new unit; a unit
already running finishes. The flag holds who paused and when, as JSON, so the page can say so. It replaces the
`paused` column of `work_lease` (control schema v21 dropped that table, PERF.100), which only the old executor ever read.

Pause is for the bulk and maintenance lanes (the units of a run, the cache warm-up, the database check). The
interactive and long lanes (`redisqueue.INTERACTIVE_QUEUE`, `LONG_QUEUE`) are never held: an admin who pauses the
generation must still be able to edit and resume.

**A solo hold (OPS.145, Boss 5:29 PM EST on 2026-11-09)** is the same flag with two more fields, `reason` and `solo` (the unit whose job is let through). After the
heartbeat's stuck action (`unitstuck`, OPS.140) stops a unit and queues it again, the follower holds every other task so the restarted unit runs alone: the workers take
no other job (`blocks`), the gate gives the unit the whole memory budget of the bulk lanes (`memorygate.Claim.solo`) and the scheduler does not charge its time to its
kind (`memorygate.release`), so the jobs held are not starved after it. Units already running finish or stop at their next batch; nothing new starts. The hold lifts when
the unit completes or fails (the follower, `lift_solo`), when its job has ended whoever noticed (`state` drops it), or when an admin presses Resume. A pause someone made
by hand is never replaced by a solo hold: the Resume it needs is theirs.

Without a Redis server nothing is paused (there is nothing to pause).
"""

import json
import time

from planetgen.queue import redisqueue, runrecord

def key():
    """The pause flag, `planetgen:queue:paused` (the tests' prefix, `runrecord.PREFIX`, in a test): a string whose value is
    `{"by": name, "at": Unix seconds}` as JSON, with no expiry."""
    return runrecord._k("queue", "paused")


NEVER_HELD = (redisqueue.INTERACTIVE_QUEUE, redisqueue.LONG_QUEUE, "planetgen-web-", "planetgen-api-")
"""tuple: The queues (by name prefix) whose workers never wait for the flag: the two API lanes, the Generate page's
job queues (a paused queue must not hold the page's own starter) and the API queues the tests name."""

WAIT_SECONDS = 2.0
"""float: How often a worker whose queue is paused looks at the flag again."""


def holds(queue_name):
    """Whether a worker of the queue `queue_name` waits while the queues are paused: every queue but `NEVER_HELD`
    (the bulk and maintenance lanes, and a run's own queue)."""
    return not queue_name.startswith(NEVER_HELD)


def state(connection=None):
    """
    The pause as it stands.

    A solo hold (`hold_solo`) whose job has ended is dropped here and reads as not paused.

    Returns:
        dict: `paused` (bool), `by` and `at` (who paused and when, `None` when not paused or not recorded), `reason` (why, `None` for a pause by hand) and `solo` (the
            unit let through, `{"job", "run", "stage", "unit"}`, or `None`).
    """
    free = {"paused": False, "by": None, "at": None, "reason": None, "solo": None}
    try:
        connection = connection or redisqueue.connect()
        raw = connection.get(key())
    except Exception:  # noqa: BLE001 -- no Redis: nothing runs, nothing is paused
        return free
    if raw is None:
        return free
    try:
        body = json.loads(raw)
        by, at, reason, solo = body.get("by"), body.get("at"), body.get("reason"), body.get("solo")
    except (ValueError, AttributeError):
        by = at = reason = solo = None
    if not isinstance(solo, dict) or not solo.get("job"):
        solo = None
    if solo is not None and _solo_is_over(connection, solo):
        _drop(connection, raw)
        return free
    return {"paused": True, "by": by, "at": at, "reason": reason, "solo": solo}


ENDED_STATUSES = ("finished", "failed", "canceled", "stopped")
"""tuple: The statuses of an RQ job that has ended, with which a solo hold on it has nothing left to hold for."""


def _solo_is_over(connection, solo):
    """Whether the job of a solo hold has ended or is gone (its follower died, or an admin cleared the run), so the hold has nothing left to wait for."""
    try:
        job = redisqueue.fetch_job(solo["job"], connection)
        if job is None:
            return True
        status = job.get_status(refresh=True)
        return str(getattr(status, "value", status)) in ENDED_STATUSES
    except Exception:  # noqa: BLE001 -- the job cannot be read: the hold stays, an admin can lift it
        return False


def _drop(connection, raw):
    """Deletes the flag if it is still what was read (a hold that was replaced meanwhile stays)."""
    try:
        if connection.get(key()) == raw:
            connection.delete(key())
    except Exception:  # noqa: BLE001 -- it expires by being read again
        pass


def solo_of(connection=None):
    """The unit a solo hold lets through (`{"job", "run", "stage", "unit"}`), or `None` when the queues are not held for one."""
    return state(connection)["solo"]


def blocks(connection=None, job_id=None):
    """
    Whether a worker must hold off now. `job_id` is the job it has taken, `None` before it takes one. A pause by hand blocks every worker that `holds`; a solo hold
    blocks every job but the unit it is held for (and no worker before it has taken one: it must be able to reach that job in the queue).
    """
    now = state(connection)
    if not now["paused"]:
        return False
    if now["solo"] is None:
        return True
    return job_id is not None and job_id != now["solo"]["job"]


def is_paused(connection=None):
    """Whether the queues are paused (`state`)."""
    return state(connection)["paused"]


def pause(who, connection=None, reason=None, solo=None):
    """
    Pauses the queues; the first pause is the one kept. `reason` and `solo` are for `hold_solo`.

    Returns:
        bool: Whether it was not paused already.
    """
    connection = connection or redisqueue.connect()
    body = {"by": str(who or "")[:64] or None, "at": time.time()}
    if reason:
        body["reason"] = str(reason)[:300]
    if solo:
        body["solo"] = solo
    return bool(connection.set(key(), json.dumps(body), nx=True))


def hold_solo(run_id, stage, unit, job_id, connection=None):
    """
    Holds the queues so that the restarted unit `unit` of stage `stage` of run `run_id` (its job `job_id`) runs alone (OPS.145). The reason names it on the Queue page, and the
    hold is by the run's follower. Not made over a pause that is on already (a pause by hand is kept, and so is a hold for another unit).

    Returns:
        bool: Whether the hold was made.
    """
    return pause(f"run {run_id}", connection, reason=f"solo run of unit {unit} (stage {stage}) after a restart",
                 solo={"job": job_id, "run": run_id, "stage": stage, "unit": unit})


def lift_solo(run_id, stage, unit, connection=None):
    """Lifts the solo hold held for that unit, if it is the one on (`False` when the queues are paused by hand or held for another unit, which are left as they are)."""
    connection = connection or redisqueue.connect()
    now = state(connection)
    solo = now["solo"]
    if solo is None or (solo.get("run"), solo.get("stage"), solo.get("unit")) != (run_id, stage, unit):
        return False
    return resume(connection)


def resume(connection=None):
    """
    Lets the workers take jobs again.

    Returns:
        bool: Whether it was paused.
    """
    return bool((connection or redisqueue.connect()).delete(key()))
