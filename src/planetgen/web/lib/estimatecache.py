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

The Redis key is `planetgen:estimate:<database>:<digest>`.
"""

import hashlib
import json

import redis

from planetgen.queue import redisclient
from planetgen.web.lib import singleflight

KEEP_SECONDS = 300
"""int: How long a stored estimate may be served: the disk's free space is in it."""

CLAIM_SECONDS = 180
"""int: How long a claim lives if its owner dies: longer than the estimate's own limit."""

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


def _read(key, versions):
    """The estimate stored under `key` while it was made under `versions`, else `None`."""
    try:
        raw = _client().get(key)
        entry = json.loads(raw) if raw else None
    except (redis.RedisError, OSError, ValueError):
        return None
    if not isinstance(entry, dict) or entry.get("versions") != versions:
        return None
    return entry.get("estimate")


def _write(key, versions, estimate):
    try:
        _client().set(key, json.dumps({"versions": versions, "estimate": estimate}, sort_keys=True),
                      ex=KEEP_SECONDS)
    except (redis.RedisError, OSError):
        pass  # the next request just runs the estimate again


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
    except Exception:  # noqa: BLE001 -- fail open: no versions, no cache
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
