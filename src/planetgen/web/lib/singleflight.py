# planetgen/web/lib/singleflight.py

"""
One build of each Galaxy Map tile or stage at a time (PERF.38).

When a map is opened or a fill makes tiles stale, several browsers and web
processes ask for the same missing tile at once, and each of them used to
build it: the API runs the same expensive query once per request, and the
database is already busy with the fill. Here a request first *claims* the
tiles it is about to fetch with a Redis `SET NX EX`; a tile somebody else
claimed is waited for instead (`wait`, polling the disk cache the claimant
writes to), and only a tile whose claimant never delivered is built after
all.

Everything fails open: with no Redis server (or one that doesn't answer in
half a second) every name counts as claimed by the caller, which is how the
cache behaved before.
"""

import time

import redis

from planetgen.queue import redisqueue

LOCK_SECONDS = 60
"""int: How long a claim lives if its owner dies without releasing it."""

WAIT_SECONDS = 20.0
"""float: The longest a request waits for the claimant of a tile before
building it itself."""

POLL_SECONDS = 0.1

SOCKET_SECONDS = 0.5
"""float: Connect and read timeout of the Redis connection: a lock service
that is slow is worth less than building the tile."""

_PREFIX = "planetgen:build:"

_clients = {}
"""dict: `url -> redis.Redis`."""


def _client():
    url = redisqueue.redis_url()
    client = _clients.get(url)
    if client is None:
        client = _clients[url] = redis.Redis.from_url(
            url, socket_timeout=SOCKET_SECONDS, socket_connect_timeout=SOCKET_SECONDS)
    return client


def claim(names):
    """
    Claims the builds `names` (strings naming a tile or stage of one
    generation).

    Returns:
        tuple: `(owned, held)`: the names this caller must build, and the
            ones another caller has claimed. All of them are `owned` when
            Redis can't be asked.
    """
    names = list(names)
    if not names:
        return [], []
    try:
        pipe = _client().pipeline(transaction=False)
        for name in names:
            pipe.set(_PREFIX + name, "1", nx=True, ex=LOCK_SECONDS)
        results = pipe.execute()
    except (redis.RedisError, OSError):
        return names, []
    owned = [name for name, got in zip(names, results) if got]
    held = [name for name, got in zip(names, results) if not got]
    return owned, held


def release(names):
    """Gives up claims made by `claim` (once the builds are on disk)."""
    names = list(names)
    if not names:
        return
    try:
        _client().delete(*[_PREFIX + name for name in names])
    except (redis.RedisError, OSError):
        pass  # the claim expires on its own


def wait(names, read):
    """
    Waits for other callers' builds of `names`.

    Args:
        names (list[str]): Claimed by others (`claim`'s `held`).
        read (callable): `read(name)` -> the finished build, or `None`
            when it isn't there yet.

    Returns:
        dict: `{name: build}` for each name that turned up within
            `WAIT_SECONDS`; the rest are the caller's to build.
    """
    found = {}
    pending = list(names)
    deadline = time.monotonic() + WAIT_SECONDS
    while pending:
        for name in list(pending):
            value = read(name)
            if value is not None:
                found[name] = value
                pending.remove(name)
        if not pending or time.monotonic() >= deadline:
            break
        time.sleep(POLL_SECONDS)
    return found
