# planetgen.queue.warmphases

"""
What the warm-up has done, phase by phase, for the admin page (PERF.141, docs/design/operations/cache-warm-up.md).

`planetgen warm` writes, per database, the Redis hash `planetgen:warm:phases:<database>`: one field per phase (`opening`,
`built-levels`, ...) holding a JSON object, and `_run`, the state of the run as a whole. The admin queue page reads it
through `GET /api/admin/work/warm` (`read`), together with what the automatic start (`warmstart`) is waiting for.

A phase is `waiting` (the run has it ahead), `running`, then one of `built` (it built something), `unchanged` (nothing
under it changed since the last run: nothing was built), `stopped` (the cache budget or a cap left something out),
`nothing` (no family of it had anything to warm) or `failed`. A phase the last run did not ask for keeps what an earlier
run wrote, so `planetgen warm lists` leaves the others as they were. Nothing here raises for Redis: the warm-up works
without its record.
"""

import json
import os
import time

import redis

from planetgen.admin import updateguard
from planetgen.db import objectindex
from planetgen.queue import redisclient, warmstart

WAITING = "waiting"
RUNNING = "running"
BUILT = "built"
UNCHANGED = "unchanged"
STOPPED = "stopped"
NOTHING = "nothing"
FAILED = "failed"
INTERRUPTED = "interrupted"
DONE = "done"
RUN_FIELD = "_run"

STATE_TEXT = {
    WAITING: "Waiting", RUNNING: "Running", BUILT: "Built", UNCHANGED: "Unchanged, nothing built",
    STOPPED: "Stopped at the cache budget or a cap", NOTHING: "Nothing to warm", FAILED: "Failed", DONE: "Done",
    INTERRUPTED: "Interrupted",
    None: "Not run yet",
}
"""dict: What the admin page says for each state (`None`: no run has recorded the phase)."""

KEEP_SECONDS = 30 * 86400
"""int: How long the record of a database nothing warmed for is kept."""


def key(database):
    """The hash of `database`'s phases."""
    return f"planetgen:warm:phases:{database}"


def _client():
    return redisclient.client(timeout=0.5, retries=1, decode_responses=True)


def _write(database, **fields):
    """Writes the JSON of each field of the hash; says nothing when Redis does not answer."""
    try:
        client = _client()
        pipe = client.pipeline(transaction=True)
        for name, value in fields.items():
            pipe.hset(key(database), name, json.dumps(value))
        pipe.expire(key(database), KEEP_SECONDS)
        pipe.execute()
    except (redis.RedisError, OSError):
        pass


def _read_field(database, name):
    try:
        raw = _client().hget(key(database), name)
        return json.loads(raw) if raw else None
    except (redis.RedisError, OSError, ValueError):
        return None


def begin(database, phases, reason=None):
    """A run starts: `phases` (names) wait, and the run is `running`."""
    now = time.time()
    fields = {name: {"state": WAITING, "queued_at": now} for name in phases}
    fields[RUN_FIELD] = {"state": RUNNING, "started_at": now, "phases": list(phases), "reason": reason,
                           "pid": os.getpid()}
    _write(database, **fields)


def start_phase(database, name):
    """Phase `name` has its first command running."""
    current = _read_field(database, name) or {}
    if current.get("state") == RUNNING:
        return
    _write(database, **{name: {"state": RUNNING, "started_at": time.time()}})


def finish_phase(database, name, result):
    """Phase `name` ended with `result` (`warmup.PhaseTally.result`: state, built, cached, left_out, seconds, line)."""
    current = _read_field(database, name) or {}
    _write(database, **{name: {**result, "started_at": current.get("started_at"), "finished_at": time.time()}})


def end(database, error=None):
    """The run ended. With `error` (a message) the phase that was running is `failed` and the phases still waiting are
    left waiting no more: they go back to what they were, `not run`."""
    run = _read_field(database, RUN_FIELD) or {}
    fields = {RUN_FIELD: {**run, "state": FAILED if error else DONE, "finished_at": time.time(), "error": error}}
    if error:
        for name in run.get("phases", []):
            current = _read_field(database, name) or {}
            if current.get("state") == RUNNING:
                fields[name] = {**current, "state": FAILED, "finished_at": time.time(), "line": error}
            elif current.get("state") == WAITING:
                fields[name] = {**current, "state": None, "line": "Not reached: the run stopped before it."}
    _write(database, **fields)


def read(database, phases):
    """
    What the admin page shows.

    Args:
        database (str): The site's database.
        phases (list[tuple]): `warmup.PHASES`: `(name, what it warms, families)`.

    Returns:
        dict: `{"run": {...} or None, "phases": [{"name", "text", "families", "state", "state_text", "line",
        "started_at", "finished_at"}], "auto": {"enabled", "wanted": [kinds], "holder": pid or None},
        "index": {"size_bytes", "kinds": [{"kind", "rows", "mode", "built_at", "seconds"}]}}` (the object index as its
        last runs recorded it, PERF.140); `available` is
        false when Redis does not answer (the phases are then all `Not run yet`).
    """
    stored, available = {}, True
    try:
        stored = _client().hgetall(key(database))
    except (redis.RedisError, OSError):
        available = False
    rows = []
    for name, text, families in phases:
        try:
            data = json.loads(stored.get(name) or "null") or {}
        except ValueError:
            data = {}
        state = data.get("state")
        rows.append({"name": name, "text": text, "families": list(families), "state": state,
                     "state_text": STATE_TEXT.get(state, str(state)), "line": data.get("line") or "",
                     "built": data.get("built"), "cached": data.get("cached"), "left_out": data.get("left_out"),
                     "seconds": data.get("seconds"), "started_at": data.get("started_at"),
                     "finished_at": data.get("finished_at")})
    try:
        run = json.loads(stored.get(RUN_FIELD) or "null")
    except ValueError:
        run = None
    if run and run.get("state") == RUNNING and not updateguard.pid_alive(run.get("pid")):
        run["state"] = INTERRUPTED       # its process is gone (killed, the machine restarted): nothing is running
        for row in rows:
            if row["state"] in (RUNNING, WAITING):
                row["state"], row["state_text"] = INTERRUPTED, STATE_TEXT[INTERRUPTED]
    auto = {"enabled": False, "wanted": [], "holder": None}
    try:
        auto["enabled"] = warmstart.enabled()
        if available:
            client = _client()
            auto["wanted"] = warmstart.pending(database, client)
            holder = client.get(warmstart.lock_key(database))
            auto["holder"] = int(holder) if holder and holder.isdigit() else None
    except (redis.RedisError, OSError):
        available = False
    recorded = objectindex.recorded(database) if available else {}
    index = {"size_bytes": objectindex.recorded_size(database) if available else 0,
             "kinds": [{"kind": kind, "rows": recorded[kind].get("rows"), "mode": recorded[kind].get("mode"),
                        "built_at": recorded[kind].get("built_at"), "seconds": recorded[kind].get("seconds")}
                       for kind in objectindex.KINDS if kind in recorded]}
    return {"available": available, "run": run, "phases": rows, "auto": auto, "index": index}
