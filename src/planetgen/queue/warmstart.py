# planetgen.queue.warmstart

"""
The warm-up that starts by itself (PERF.141, docs/design/operations/cache-warm-up.md).

Every change that moves a data version tag (`planetgen.db.dataversion.bump`: a generation run, `planetgen orbits`, an
edit through the API, a change the page cache found) calls `request`, which notes the change in Redis and starts
`planetgen warm --on-change DATABASE` in the background unless one already waits or runs for that database. That
process waits `cache.warm_settle_seconds`, so a burst of edits is one warm-up, then takes what was asked for and runs the
warm-up once (a second run after a run that was asked again meanwhile, at most `MAX_ROUNDS`), which rebuilds only what
the changes touched (PERF.133) and so costs almost nothing when little changed.

Redis keys, per database: `planetgen:warm:wanted:<database>`, a hash of the kinds of change asked for (kind to the time
of the last), and `planetgen:warm:auto:<database>`, held by the process that waits or runs (its process id, kept for
`LOCK_SECONDS`; one whose process is gone is taken over). A request made while one is held only adds to the hash: the
holder looks again when it finishes, after letting go, so a change is never left unwarmed.

Nothing starts while an update runs (the update starts its own warm-up at its end), with `cache.warm_on_change` off
(`PLANETGEN_WARM_ON_CHANGE=0`), or without Redis; a request that cannot start anything says so on stderr and never
raises, since a data change must not fail because its caches could not be warmed.
"""

import os
import subprocess
import sys
import time

import redis

from planetgen.admin import updateguard
from planetgen.queue import redisclient

MAX_ROUNDS = 5
"""int: The most warm-up runs one automatic start makes; changes that keep arriving wait for the next request."""

LOCK_SECONDS = 4 * 3600
"""int: How long the lock of a waiting or running automatic warm-up is kept (its process's death frees it sooner)."""

WANTED_SECONDS = 86400
"""int: How long a request is kept when nothing takes it."""


FULL_INDEX_TAGS = ("positions", "epoch")
"""tuple: The changes that move what a refresh cannot see by `modified_at` (an orbital update leaves it alone and moves
which sector holds an object; `everything` is a release or a plan change): the object index is built whole for them."""


def wanted_key(database):
    """The hash of the changes asked for on `database`."""
    return f"planetgen:warm:wanted:{database}"


def lock_key(database):
    """The key held, with its process id, by the automatic warm-up of `database` while it waits or runs."""
    return f"planetgen:warm:auto:{database}"


def _client():
    return redisclient.client(timeout=0.5, retries=1, decode_responses=True)


def enabled():
    """Whether changes start the warm-up by themselves (`cache.warm_on_change`)."""
    from planetgen.util.settings import get_settings
    return bool(get_settings().cache.warm_on_change)


def settle_seconds():
    """How long the automatic warm-up waits for more changes (`cache.warm_settle_seconds`)."""
    from planetgen.util.settings import get_settings
    return float(get_settings().cache.warm_settle_seconds)


def pending(database, client=None):
    """The kinds of change asked for and not yet taken, as a list."""
    return sorted((client or _client()).hkeys(wanted_key(database)))


def take(database, client=None):
    """Takes what was asked for on `database`: the kinds of change, and the request is gone."""
    client = client or _client()
    pipe = client.pipeline(transaction=True)
    pipe.hkeys(wanted_key(database))
    pipe.delete(wanted_key(database))
    return sorted(pipe.execute()[0])


def acquire(database, pid=None, client=None):
    """Takes the lock for process `pid` (default this one); a lock whose process is gone is taken over. Returns whether it is ours."""
    client = client or _client()
    pid = os.getpid() if pid is None else pid
    if client.set(lock_key(database), pid, nx=True, ex=LOCK_SECONDS):
        return True
    holder = client.get(lock_key(database))
    try:
        alive = updateguard.pid_alive(int(holder))
    except (TypeError, ValueError):
        alive = holder is not None
    if holder is not None and not alive:
        client.delete(lock_key(database))
        return bool(client.set(lock_key(database), pid, nx=True, ex=LOCK_SECONDS))
    return False


def release(database, pid=None, client=None):
    """Lets go of the lock when it is this process's (default) or `pid`'s."""
    client = client or _client()
    pid = os.getpid() if pid is None else pid
    if client.get(lock_key(database)) == str(pid):
        client.delete(lock_key(database))


def _cache_owner():
    """`(uid, gid)` of the tile cache's directory when this process is root and the directory belongs to another user
    (Apache's), else `None`: the warm-up runs as the user the web server reads the files as."""
    if os.geteuid() != 0:
        return None
    from planetgen.web.lib import tilecache
    try:
        stat = os.stat(tilecache.cache_dir())
    except (OSError, TypeError):
        return None
    return (stat.st_uid, stat.st_gid) if stat.st_uid != 0 else None


def spawn(database):
    """Starts `planetgen warm --on-change DATABASE` apart from this process, as the cache directory's owner when this
    process is root, and returns its process id."""
    argv = [sys.executable, "-m", "planetgen.cli.generate", "warm", "--on-change", database]
    options = {}
    owner = _cache_owner()
    if owner is not None:
        options = {"user": owner[0], "group": owner[1]}
    proc = subprocess.Popen(argv, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                            cwd="/", close_fds=True, start_new_session=True, **options)
    return proc.pid


def request(database, kind):
    """
    Notes that `kind` of change happened on `database` and starts the warm-up unless one already waits or runs.

    Returns:
        str: `off` (not enabled), `update` (an update runs and will warm), `running` (the change was added to the
        warm-up that waits or runs), `started`, or `failed` (Redis or the start failed; warned on stderr).
    """
    if database is None:
        return "off"
    try:
        if not enabled():
            return "off"
        if updateguard.is_updating():
            return "update"
        client = _client()
        pipe = client.pipeline(transaction=True)
        pipe.hset(wanted_key(database), kind, int(time.time()))
        pipe.expire(wanted_key(database), WANTED_SECONDS)
        pipe.execute()
        if not acquire(database, client=client):
            return "running"
        try:
            pid = spawn(database)
        except OSError:
            release(database, client=client)
            raise
        client.set(lock_key(database), pid, xx=True, keepttl=True)
        return "started"
    except (redis.RedisError, OSError) as exc:
        print(f"warning: the cache warm-up did not start after a {kind!r} change ({exc}); run `planetgen warm`.",
              file=sys.stderr)
        return "failed"


def drain(database, run, settle=None, max_rounds=MAX_ROUNDS, sleep=time.sleep, client=None):
    """
    The automatic warm-up's loop, in the process `request` started (which holds the lock): waits `settle` seconds,
    takes the changes asked for, calls `run(kinds)`, and repeats while more were asked for meanwhile. It lets go of
    the lock before it looks, so a change that arrives as it ends is run by it or by the request's own start.

    Args:
        database (str): The database the changes are on.
        run (callable): Warms once; called with the sorted kinds of change. May raise: the lock is let go first and
            the error goes on.
        settle (float, optional): Seconds to wait before each take (default `cache.warm_settle_seconds`).
        max_rounds (int): The most runs.
        sleep (callable): The wait, for tests.

    Returns:
        int: How many times `run` was called.
    """
    client = client or _client()
    settle = settle_seconds() if settle is None else settle
    rounds = 0
    try:
        while True:
            sleep(settle)
            if rounds >= max_rounds:
                break
            kinds = take(database, client)
            if kinds:
                rounds += 1
                run(kinds)
                continue
            release(database, client=client)
            if not pending(database, client) or not acquire(database, client=client):
                return rounds
        release(database, client=client)
        return rounds
    except BaseException:
        release(database, client=client)
        raise
