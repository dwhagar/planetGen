# planetgen/web/lib/estimatecache.py

"""
The Generate page's size-and-time estimate, cached in front of the work (PERF.109).

`planetgen galaxy ... --estimate-only` is repeatable and costly (up to two minutes for a large shell), and the
same form is sent again each time an admin opens the estimate, reloads it or confirms it. The answer is kept in
Redis under the kind of job, a digest of its command line and the galaxy database, with the versions it was made
under:

- the generator epoch (`util/generator_epoch.py`),
- the galaxy's data version (`db/dataversion.py`: the incarnation and the count every data change bumps, `epoch`
  tag included, so a fill or an edit makes the next estimate a new one), and
- the speed statistics' version (`generation/stats.py` `fingerprint`: the recorded samples and sizes the time
  and bytes are worked out from).

A stored answer is served only while all of them still stand. It is also dropped after `KEEP_SECONDS`: the
free space of the disk is part of the answer and nothing bumps a version when it changes (the job itself
checks the disk again when it starts).

Two admins asking for the same estimate at once run it once: the first request claims the work with
`singleflight` and the other waits for the stored answer, building it after all only when the claimant never
delivers (`WAIT_SECONDS`). A request that cannot read the versions or reach Redis simply runs the estimate, as
the page did before.

On the web (ADM.141) the estimate is never worked out inside the request: a big galaxy takes minutes, and Apache
restarts the whole daemon of a request that outlives `request-timeout` (60 s). `estimate_in_background` answers from
the cache, or hands the work to the work queue (`web/lib/estimatejob.py`) and raises `Working` at once; the page
asks again every few seconds and the finished answer (or its failure, kept `FAILURE_SECONDS`) is read from the same
key. One claim for the question and one for the database keep a second copy of the command from starting while one
runs, so a form changed or a page reloaded mid-estimate never starts the command again.

The Redis key is `planetgen:estimate:<database>:<digest>`.
"""

import hashlib
import json

import redis

from planetgen.queue import redisclient
from planetgen.util import log
from planetgen.web.lib import singleflight

KEEP_SECONDS = 300
"""int: How long a stored estimate may be served: the disk's free space is in it."""

CLAIM_SECONDS = 660
"""int: How long a claim lives if its owner dies: longer than the estimate's own limit (`estimatejob.TIMEOUT_SECONDS`)."""

FAILURE_SECONDS = 30
"""int: How long the failure of an estimate is kept, so a page that asks again at once reads why instead of
starting the same failing command again."""

GATE_PREFIX = "gate:"
"""str: The claim that lets one estimate at a time run for a galaxy database (ADM.141): the command reads the whole
plan, and a form being typed into would otherwise start one for every value on the way."""

WAIT_SECONDS = 150.0
"""float: The longest a request waits for the claimant of the same estimate."""

SOCKET_SECONDS = 0.5
"""float: How long one Redis call may take before the estimate is simply run."""

_PREFIX = "planetgen:estimate:"


def cache_key(database, kind, argv):
    """The Redis key of the estimate of the job `kind` with command line `argv` in `database`."""
    digest = hashlib.sha256(json.dumps([kind, list(argv)], sort_keys=True).encode("utf-8")).hexdigest()[:32]
    return f"{_PREFIX}{database}:{digest}"


def _client():
    return redisclient.client(timeout=SOCKET_SECONDS, retries=1)


def _entry(key, versions):
    """What is stored under `key` while it was made under `versions` (`{"versions", "estimate"}` or `{"versions",
    "error"}`), else `None`."""
    try:
        raw = _client().get(key)
        entry = json.loads(raw) if raw else None
    except (redis.RedisError, OSError, ValueError) as exc:
        log.degraded("estimatecache", "reading an estimate from Redis failed, so it is made again", exc)
        return None
    if not isinstance(entry, dict) or entry.get("versions") != versions:
        return None
    return entry


def _read(key, versions):
    """The estimate stored under `key` while it was made under `versions`, else `None`."""
    entry = _entry(key, versions)
    return None if entry is None else entry.get("estimate")


def _write(key, versions, estimate):
    try:
        _client().set(key, json.dumps({"versions": versions, "estimate": estimate}, sort_keys=True),
                      ex=KEEP_SECONDS)
    except (redis.RedisError, OSError) as exc:
        log.degraded("estimatecache", "keeping an estimate in Redis failed", exc)  # the next request just runs the estimate again


def store(key, versions, estimate=None, error=None):
    """Keeps what the work queue worked out under `key` (the estimate, or `error`: why it could not be made, for
    `FAILURE_SECONDS`) for `estimate_in_background` to read, and gives up the claims of the work."""
    if error is None:
        _write(key, versions, estimate)
    else:
        try:
            _client().set(key, json.dumps({"versions": versions, "error": str(error)}, sort_keys=True),
                          ex=FAILURE_SECONDS)
        except (redis.RedisError, OSError) as exc:
            log.degraded("estimatecache", "keeping an estimate failure in Redis failed", exc)
    singleflight.release([key, _gate_of(key)])


def _gate_of(key):
    """The database-wide claim that goes with the estimate stored under `key`."""
    database = key[len(_PREFIX):].split(":", 1)[0]
    return f"{GATE_PREFIX}{database}"


def estimate(database, kind, argv, versions, compute):
    """
    The estimate for a job, from the cache when one made under the current versions is stored.

    Args:
        database (str): The galaxy database.
        kind (str): The kind of job (the form's `build_job` kind).
        argv (list): The job's command line, the estimate's arguments.
        versions (callable): `versions()` -> a JSON-able value that changes when the estimate's answer would
            (see the module note). Raising anything means "can't tell": `compute` runs and nothing is stored.
        compute (callable): `compute()` -> the estimate dict; it raises when the estimate can't be made, and
            then nothing is stored.

    Returns:
        dict: The estimate, a fresh copy.
    """
    try:
        current = json.loads(json.dumps(versions(), sort_keys=True))
    except Exception as exc:  # noqa: BLE001 -- fail open: no versions, no cache
        log.swallowed("web/lib/estimatecache", "estimate", exc)
        return compute()
    key = cache_key(database, kind, argv)
    stored = _read(key, current)
    if stored is not None:
        return stored
    owned, held = singleflight.claim([key], seconds=CLAIM_SECONDS)
    if held:
        found = singleflight.wait(held, lambda name: _read(name, current), seconds=WAIT_SECONDS)
        if key in found:
            return found[key]
    try:
        result = compute()
        _write(key, current, result)
        return result
    finally:
        if owned:
            singleflight.release(owned)


class Working(Exception):
    """The estimate is being worked out in the work queue; ask again in a few seconds."""


class Failed(Exception):
    """The estimate could not be made; the message says why."""


def estimate_in_background(database, kind, argv, versions, compute, start):
    """
    The estimate for a job without ever waiting for it (ADM.141): from the cache when one made under the current
    versions is stored, else handed to `start` (the work queue) and `Working` raised, so the request returns at
    once and the page asks again.

    Args:
        database (str): The galaxy database.
        kind (str): The kind of job.
        argv (list): The job's command line, the estimate's arguments.
        versions (callable): As in `estimate`.
        compute (callable): `compute()` -> the estimate dict, run in this request only when nothing can queue it
            (no Redis server: the page then waits, as it did before the queue); it raises when it cannot.
        start (callable): `start(key, versions)` -> whether the work was queued. The worker ends by `store`, which
            also gives up the claims.

    Returns:
        dict: The estimate, a fresh copy.

    Raises:
        Working: The estimate is running (queued by this call, or by an earlier one).
        Failed: The last try at it failed (the message says why).
    """
    try:
        current = json.loads(json.dumps(versions(), sort_keys=True))
    except Exception as exc:  # noqa: BLE001 -- fail open: no versions, no cache
        log.swallowed("web/lib/estimatecache", "estimate_in_background", exc)
        return compute()
    key = cache_key(database, kind, argv)
    entry = _entry(key, current)
    if entry is not None:
        if "estimate" in entry:
            return entry["estimate"]
        raise Failed(entry.get("error") or "the estimate could not be made")
    gate = _gate_of(key)
    owned, held = singleflight.claim([key, gate], seconds=CLAIM_SECONDS)
    if held:
        singleflight.release(owned)        # the same question is running, or another one is: wait for it
        raise Working()
    queued = False
    try:
        queued = bool(start(key, current))
    except Exception as exc:  # noqa: BLE001 -- no queue: the estimate is worked out here
        log.swallowed("web/lib/estimatecache", "estimate_in_background", exc)
        queued = False
    if queued:
        raise Working()
    try:
        result = compute()
        _write(key, current, result)
        return result
    finally:
        singleflight.release(owned)
