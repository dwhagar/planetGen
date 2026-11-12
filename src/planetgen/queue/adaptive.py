# planetgen/queue/adaptive.py

"""
The number of bulk workers follows the host's memory (adaptive worker control).

The memory gate (`memorygate`) holds a job at the gate until the host has room for it, and lets it in after `memorygate.MAX_WAIT_SECONDS` by the
table cap. A generation run on a host with more workers than its memory carries therefore has workers that sit at the gate, each holding the star
tables (about 1.4 GB), for minutes at a time: the run reports "waited 300 s" every few systems and goes slowly, and the waiting workers are
part of what keeps the room from coming free. The cap on the bulk workers (`workers.wanted`, `workers.ensure_workers`) was a number; this file
makes it follow the host.

**Lowering.** The limit goes down by one worker when

* a bulk job has waited `LOWER_AFTER_SECONDS` (30 s) at the gate for *memory* (`memorygate.Decision.memory`; a wait for a core, or for a higher
  ticket to go first, is not a memory wait), or
* the host's available memory is under `PRESSURE_SHARE` of the memory the gate keeps back (`memorygate.reserve_mb`) and a bulk job is waiting
  or running: the kernel is about to choose a victim, before any job has waited long enough to say so.

A lowering asks the bulk workers over the new limit to finish and leave (`workers.shed_bulk`), the ones that hold a job back at the gate first:
such a worker puts its job back at the front of its queue and exits, so nothing is lost and the memory of its tables comes free. Two lowerings are
`LOWER_COOLDOWN_SECONDS` apart so the first has had its effect before the next is judged. The limit is never below `MIN_WORKERS` (1).

**Raising.** The limit goes up by one (`RAISE_STEP_SECONDS` apart) when no memory wait or lowering has happened for `RAISE_QUIET_SECONDS`, bulk jobs are
queued that no worker has, and the memory the gate would hand out (`memorygate.decide`'s own free figure) holds one more worker of the largest
claim held with `RAISE_HEADROOM` on top. When the limit reaches the configured cap it is dropped, and a limit not touched for `LIMIT_SECONDS` is forgotten
(a run that ended leaves nothing behind).

**Every change is logged** (a lowering as a warning line with the job that waited and what it waited for, a raise as a normal line) and kept in a short list
in Redis (`recent`) for a status page. Nothing here changes what a run generates: the work seed does not depend on how many workers share the work
(`work.task_seed`). The control is off with the gate (`memorygate.DISABLE_ENV`).
"""

import json
import time
import uuid

from planetgen.queue import memorygate, workers
from planetgen.util import log

MIN_WORKERS = 1
"""int: The limit never goes below this."""

LOWER_AFTER_SECONDS = 30.0
"""float: How long a bulk job waits for memory before the limit goes down (Boss, 2026-11-11)."""

LOWER_COOLDOWN_SECONDS = 20.0
"""float: The least time between two lowerings: long enough for the worker asked to leave to give its memory back."""

PRESSURE_SHARE = 0.75
"""float: The host is under pressure when its available memory is under this share of `memorygate.reserve_mb`. Stand-in: a design edge, not a measurement."""

RAISE_QUIET_SECONDS = 120.0
"""float: How long without a memory wait or a lowering before the limit may go up."""

RAISE_STEP_SECONDS = 120.0
"""float: The least time between two raises."""

RAISE_HEADROOM = 1.25
"""float: A raise needs free memory for one more worker of the largest bulk claim held, times this."""

LIMIT_SECONDS = 3600
"""int: A limit nobody touches for this long is forgotten."""

RECENT_KEEP = 50
"""int: How many changes the Redis list keeps."""

LOCK_SECONDS = 5


def _key(name):
    return memorygate._key(f"adaptive:{name}")


def _now(now):
    return time.time() if now is None else now


def _state(conn):
    return {memorygate._decode(k): memorygate._decode(v) for k, v in conn.hgetall(_key("state")).items()}


def limit_of(conn):
    """The limit now, or `None` when the control holds the workers to nothing but the configured cap."""
    value = conn.hget(_key("state"), "limit")
    return None if value is None else max(MIN_WORKERS, int(memorygate._decode(value)))


def recent(conn, count=RECENT_KEEP):
    """The latest changes of the limit, newest first (`{"t", "from", "to", "why"}`), for a status page."""
    return [json.loads(memorygate._decode(raw)) for raw in conn.lrange(_key("recent"), 0, count - 1)]


def _record(conn, before, after, why, now):
    conn.lpush(_key("recent"), json.dumps({"t": now, "from": before, "to": after, "why": why}))
    conn.ltrim(_key("recent"), 0, RECENT_KEEP - 1)
    conn.expire(_key("recent"), LIMIT_SECONDS * 24)


def _locked(conn):
    """Takes the control's lock (one change at a time); returns its token or `None`."""
    token = uuid.uuid4().hex
    return token if conn.set(_key("lock"), token, nx=True, ex=LOCK_SECONDS) else None


def _unlock(conn, token):
    if memorygate._decode(conn.get(_key("lock"))) == token:
        conn.delete(_key("lock"))


def _set(conn, limit, now, field):
    conn.hset(_key("state"), mapping={"limit": int(limit), field: now})
    conn.expire(_key("state"), LIMIT_SECONDS)


def lower(conn, why, first=(), now=None):
    """
    Lowers the limit by one worker below the bulk workers that stay (never under `MIN_WORKERS`), logs it, and asks the workers over the new limit to
    leave (`workers.shed_bulk`, `first` before any other). Returns the new limit, or `None` when nothing changed (too soon after the last
    lowering, only one worker, or another caller holds the lock).
    """
    now = _now(now)
    token = _locked(conn)
    if token is None:
        return None
    try:
        state = _state(conn)
        if now - float(state.get("lowered_at", 0.0)) < LOWER_COOLDOWN_SECONDS:
            return None
        staying = len(workers.staying_bulk_workers(conn))
        current = limit_of(conn)
        before = staying if current is None else min(current, staying)
        after = max(MIN_WORKERS, before - 1)
        if after >= before:
            return None
        _set(conn, after, now, "lowered_at")
        _record(conn, before, after, why, now)
        log.warning(f"Worker limit lowered from {before} to {after}: {why}.")
        workers.shed_bulk(conn, after, first=first)
        return after
    finally:
        _unlock(conn, token)


def on_memory_wait(conn, name, kind, memory_waited, now=None):
    """The gate's call: worker `name` has held a job of `kind` back for `memory_waited` seconds for memory. At `LOWER_AFTER_SECONDS` the limit goes down."""
    now = _now(now)
    if memory_waited < LOWER_AFTER_SECONDS:
        return None
    conn.hset(_key("state"), "memory_wait_at", now)
    from planetgen.queue import scheduler
    waiting = [ticket.name for ticket in scheduler.tickets(conn, now) if ticket.claim.lane == memorygate.LANE_BULK]
    return lower(conn, f"{kind} has waited {memory_waited:,.0f} s for memory", first=[name, *waiting], now=now)


def pressure(host):
    """Whether the host's available memory is under `PRESSURE_SHARE` of what the gate keeps back."""
    return host.available_mb < PRESSURE_SHARE * memorygate.reserve_mb(host)


def adapt(conn, cap, queued, host=None, held=None, now=None):
    """
    The limit of the bulk workers (called with the cap `workers.ensure_workers` worked out and the bulk jobs queued): `cap` less what the host's
    memory says, at least `MIN_WORKERS`. Lowers it when the host is under pressure, raises it when memory has been free for a while
    (see the module). Never raises an error: on any failure the cap stands.
    """
    if not memorygate.enabled():
        return cap
    try:
        return _adapt(conn, cap, queued, host, held, _now(now))
    except Exception as exc:  # noqa: BLE001 -- an optimisation must never keep a run from getting workers
        log.debug("adaptive worker control failed: %s", exc, area="queue")
        return cap


def _adapt(conn, cap, queued, host, held, now):
    host = host or memorygate.read_host()
    held = memorygate.holders(conn) if held is None else held
    bulk = [holder for holder in held if holder.claim.lane == memorygate.LANE_BULK]
    state = _state(conn)
    limit = limit_of(conn)
    if limit is not None and limit >= cap:
        conn.delete(_key("state"))              # the cap is the limit: nothing to hold back
        state, limit = {}, None
    if pressure(host) and (bulk or queued):
        lower(conn, f"the host has {host.available_mb:,.0f} MB available, under {PRESSURE_SHARE * memorygate.reserve_mb(host):,.0f} MB (the gate keeps "
                    f"{memorygate.reserve_mb(host):,.0f} MB back)", now=now)
        limit = limit_of(conn)
    elif limit is not None and limit < cap:
        live = len(workers.staying_bulk_workers(conn))
        last = max(float(state.get("lowered_at", 0.0)), float(state.get("memory_wait_at", 0.0)), float(state.get("raised_at", 0.0)))
        quiet = now - max(float(state.get("lowered_at", 0.0)), float(state.get("memory_wait_at", 0.0)))
        biggest = max([holder.claim.mb for holder in bulk] or [float(workers.BULK_WORKER_MEMORY_MB)])
        unrealised = sum(max(0.0, holder.claim.mb - holder.held_mb) for holder in held if holder.claim.lane not in memorygate.SHORT_LANES)
        free = host.available_mb - memorygate.reserve_mb(host) - unrealised
        if quiet >= RAISE_QUIET_SECONDS and now - last >= RAISE_STEP_SECONDS and queued > live and free >= RAISE_HEADROOM * biggest:
            token = _locked(conn)
            if token is not None:
                try:
                    after = limit + 1
                    _record(conn, limit, after, "memory has been free", now)
                    log.normal(f"Worker limit raised from {limit} to {after}: no job has waited for memory for {quiet:,.0f} s and "
                               f"{free:,.0f} MB are free for a worker of {biggest:,.0f} MB.")
                    if after >= cap:
                        conn.delete(_key("state"))
                    else:
                        _set(conn, after, now, "raised_at")
                    limit = after
                finally:
                    _unlock(conn, token)
    return cap if limit is None else max(MIN_WORKERS, min(cap, limit))
