# planetgen/db/countcache.py

"""
Counts and filter-menu totals that are too big to take inside a web
request (PERF.64).

A table page used to count its whole table on every load (the rows
matching its filters, then each filter menu's options): about 3 s per
million systems per count, so a galaxy of ten million systems took half a
minute for the Systems list alone and tripped the 10 s statement limit
whenever the database was also busy with a fill. `cached` turns each of
those into a stored answer:

- An answer made under the galaxy's current content stamp is served as it
  is.
- When the stamp has moved on (a fill adds sectors every second), the
  older answer is served and one background thread, on a connection with
  no time limit, makes the next. At most one is made at a time per answer
  and per database server (a Redis `SET NX EX`), and not again sooner than
  `MIN_INTERVAL_SECONDS` or three times as long as the last one took, so a
  fill's constant changes cost the database one scan at a time, not one
  per page load.
- With no answer at all (the first load after a restart, or of a filter
  nobody has used) the request waits up to `WAIT_SECONDS` for the first
  one and otherwise returns the caller's `fallback`: a cheap estimate.

The answers are kept in this process and, when Redis answers, in Redis
(`PLANETGEN_COUNT_CACHE=off` turns the cache off; the tests use that, then
test it on). Everything fails open: a Redis or thread problem just means
the caller's fallback or an uncached count.
"""

import hashlib
import json
import threading
import time

import redis

from planetgen.queue import redisqueue
from planetgen.util.settings import get_settings

ENV_VAR = "PLANETGEN_COUNT_CACHE"
"""str: `off` makes `cached` just call `compute` (the tests, and anyone who prefers exact counts at any price);
the settings model's `page_cache.stored_counts`."""

WAIT_SECONDS = 3.0
"""float: The longest a request waits for the first answer."""

MIN_INTERVAL_SECONDS = 30.0
"""float: The least time between two background counts of one answer."""

MAX_AGE_SECONDS = 6 * 3600
"""int: An answer older than this is not served (it counts as missing)."""

MAX_ENTRIES = 400
"""int: Answers kept in this process."""

LOCK_SECONDS = 600
"""int: How long a count's Redis claim lives if its owner dies."""

SOCKET_SECONDS = 0.5

_memory = {}
"""dict: `(database, key) -> {"value", "stamp", "made", "started", "took"}`."""

_running = {}
"""dict: `(database, key) -> threading.Event`, set when that background count ends."""

_lock = threading.Lock()
_clients = {}


def enabled():
    """`page_cache.stored_counts` (`PLANETGEN_COUNT_CACHE`)."""
    return get_settings().page_cache.stored_counts


def _database(conn):
    config = getattr(conn, "_config", None)
    return None if config is None else f"{config.host}:{config.port}:{config.database}"


def _redis():
    url = redisqueue.redis_url()
    client = _clients.get(url)
    if client is None:
        client = _clients[url] = redis.Redis.from_url(
            url, socket_timeout=SOCKET_SECONDS, socket_connect_timeout=SOCKET_SECONDS)
    return client


def _redis_key(database, key, suffix=""):
    digest = hashlib.sha256(json.dumps(key, sort_keys=True, default=str).encode("utf-8")).hexdigest()[:24]
    return f"planetgen:counts:{database}:{digest}{suffix}"


def _recall(database, key):
    """The stored answer for `key`, from this process or Redis, or `None`."""
    with _lock:
        entry = _memory.get((database, key))
    if entry is not None:
        return entry
    try:
        raw = _redis().get(_redis_key(database, key))
    except (redis.RedisError, OSError):
        return None
    if raw is None:
        return None
    try:
        entry = json.loads(raw)
    except ValueError:
        return None
    with _lock:
        _memory[(database, key)] = entry
    return entry


def _store(database, key, entry):
    with _lock:
        _memory[(database, key)] = entry
        while len(_memory) > MAX_ENTRIES:
            _memory.pop(next(iter(_memory)))
    try:
        _redis().set(_redis_key(database, key), json.dumps(entry), ex=MAX_AGE_SECONDS)
    except (redis.RedisError, OSError, TypeError, ValueError):
        pass


def _claim(database, key):
    """Whether this caller may make the answer now (not another process)."""
    try:
        return bool(_redis().set(_redis_key(database, key, ":lock"), "1", nx=True, ex=LOCK_SECONDS))
    except (redis.RedisError, OSError):
        return True  # no lock service: the in-process claim alone


def _release(database, key):
    try:
        _redis().delete(_redis_key(database, key, ":lock"))
    except (redis.RedisError, OSError):
        pass


def _make(config, database, key, compute, stamp, done):
    from planetgen.db import query

    started = time.time()
    try:
        conn = query.open_readonly(config, statement_timeout_s=None)
        try:
            value = compute(conn)
        finally:
            conn.close()
        json.dumps(value)
        _store(database, key, {"value": value, "stamp": stamp, "made": time.time(), "started": started,
                               "took": time.time() - started})
    except Exception:  # noqa: BLE001 -- a count that fails is simply not stored; the next request tries again
        pass
    finally:
        _release(database, key)
        with _lock:
            _running.pop((database, key), None)
        done.set()


def _start(conn, database, key, compute, stamp, entry):
    """Starts the background count when none is running and the last one is old enough.

    Returns:
        threading.Event or None: Set when the count ends; `None` when none was started (a recent or
            running one stands).
    """
    now = time.time()
    if entry is not None:
        wait = max(MIN_INTERVAL_SECONDS, 3 * float(entry.get("took") or 0))
        if now - float(entry.get("started") or 0) < wait:
            return None
    with _lock:
        if (database, key) in _running:
            return _running[(database, key)]
        done = threading.Event()
        _running[(database, key)] = done
    if not _claim(database, key):
        with _lock:
            _running.pop((database, key), None)
        return None
    thread = threading.Thread(target=_make, args=(conn._config, database, key, compute, stamp, done),
                              name="planetgen-count", daemon=True)
    thread.start()
    return done


def cached(conn, key, compute, fallback, stamp):
    """
    One stored count or menu of counts (see this module's docstring).

    Args:
        conn (planetgen.db.store.Connection): The request's connection (its timeout applies to `fallback` only).
        key: A JSON-able name of the answer: what is counted and every filter.
        compute (callable): `compute(conn)` makes the exact answer, on a connection of its own.
        fallback (callable): `fallback(conn)` makes a cheap approximate answer on the request's connection.
        stamp (callable): `stamp(conn)` is the galaxy's content stamp.

    Returns:
        The answer: exact when made under the current stamp, else the last one made, else `fallback`'s.
    """
    database = _database(conn)
    if not enabled() or database is None:
        return compute(conn)
    key = json.dumps(key, sort_keys=True, default=str)
    current = stamp(conn)
    entry = _recall(database, key)
    if entry is not None and time.time() - float(entry.get("made") or 0) > MAX_AGE_SECONDS:
        entry = None
    if entry is not None and entry.get("stamp") == current:
        return entry["value"]
    done = _start(conn, database, key, compute, current, entry)
    if entry is not None:
        return entry["value"]
    if done is not None:
        done.wait(WAIT_SECONDS)
        fresh = _recall(database, key)
        if fresh is not None:
            return fresh["value"]
    return fallback(conn)


def forget_all():
    """Drops every answer kept in this process (tests)."""
    with _lock:
        _memory.clear()
