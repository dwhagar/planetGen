# planetgen/queue/redisclient.py

"""
The one place a Redis client is made (PERF.94, finding F9 of
docs/design/redis-and-the-work-queue.md).

Every caller asks `client` for its connection and gets a shared, cached one
for that server and those settings, so the web process, the API jobs, the
single-flight locks, the lockout store and the work queue no longer each
build their own pool with their own settings (or none):

* a connect and a read timeout, so a hung Redis fails a call quickly instead
  of hanging the web thread that made it (`SOCKET_SECONDS`; a caller that
  would rather build the tile itself asks for less);
* a retry with a short backoff on a dropped connection or a timeout
  (`RETRIES`), which also covers a Redis that restarts between two calls;
* a health check on a connection idle for `HEALTH_CHECK_SECONDS`, so a stale
  one is replaced before use;
* TCP keepalive;
* one pool per distinct setting, kept for the life of the process
  (redis-py's pools are thread-safe and fork-aware).

A client that blocks on purpose (`blpop`) asks for `blocking=True`: no read
timeout, since the server holds the call, with keepalive and the health
check noticing a dead server.

Workers keep RQ's own settings for their connection (RQ raises the read
timeout above its blocking time itself); this module is for the clients.
"""

import threading

import redis
from redis.backoff import ExponentialBackoff
from redis.retry import Retry

SOCKET_CONNECT_SECONDS = 2.0
"""float: How long a connect may take."""

SOCKET_SECONDS = 5.0
"""float: How long one command may take to answer."""

RETRIES = 2
"""int: Extra tries after a dropped connection or a timeout, each after a short, growing pause."""

HEALTH_CHECK_SECONDS = 30
"""int: A connection idle this long is pinged before its next command."""

_clients = {}
_lock = threading.Lock()


def default_url():
    """The Redis server in the settings: `PLANETGEN_REDIS_URL`, else `config.json`'s `redis.url`."""
    from planetgen.util.settings import get_settings

    return get_settings().redis.url


def _build(url, timeout, retries, decode_responses, blocking):
    connect_timeout = SOCKET_CONNECT_SECONDS if timeout is None else min(SOCKET_CONNECT_SECONDS, timeout)
    read_timeout = None if blocking else (SOCKET_SECONDS if timeout is None else timeout)
    options = {
        "socket_connect_timeout": connect_timeout, "socket_timeout": read_timeout, "socket_keepalive": True,
        "health_check_interval": HEALTH_CHECK_SECONDS, "decode_responses": decode_responses,
    }
    if retries:
        options["retry"] = Retry(ExponentialBackoff(cap=0.5, base=0.05), retries)
        options["retry_on_error"] = [redis.ConnectionError, redis.TimeoutError]
    return redis.Redis.from_url(url, **options)


def client(url=None, *, timeout=None, retries=RETRIES, decode_responses=False, blocking=False):
    """
    The shared client for the server at `url` with these settings.

    Args:
        url (str, optional): The server (default `default_url()`).
        timeout (float, optional): Connect and read timeout in seconds
            (default `SOCKET_CONNECT_SECONDS` and `SOCKET_SECONDS`).
        retries (int): Extra tries after a dropped connection or a timeout; 0 for none.
        decode_responses (bool): Strings back instead of bytes.
        blocking (bool): For a caller that waits on the server (`blpop`): no read timeout.

    Returns:
        redis.Redis: The same object for the same arguments.
    """
    url = url or default_url()
    key = (url, timeout, retries, decode_responses, blocking)
    found = _clients.get(key)
    if found is None:
        with _lock:
            found = _clients.get(key)
            if found is None:
                found = _clients[key] = _build(url, timeout, retries, decode_responses, blocking)
    return found


def healthy(url=None, **options):
    """Whether the server answers a ping (with the retries `client` gives)."""
    try:
        return bool(client(url, **options).ping())
    except (redis.RedisError, OSError):
        return False


def forget_all():
    """Drops every cached client and closes its pool (for tests that point the settings elsewhere)."""
    with _lock:
        clients = list(_clients.values())
        _clients.clear()
    for found in clients:
        try:
            found.close()
        except Exception:  # noqa: BLE001 -- a pool that will not close is dropped anyway
            pass
