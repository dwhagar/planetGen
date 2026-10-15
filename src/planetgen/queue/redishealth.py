# planetgen/queue/redishealth.py

"""
What the admin Stats page shows about Redis (ADM.54): whether it answers, and the settings and figures that
decide whether the work queue survives a restart (docs/design/redis-and-the-work-queue.md, section 2.3): the
append-only file, the memory policy, the memory used, the last save, and how many workers are alive.

Read with `INFO` only, so it works against a server that has `CONFIG` switched off, and with a short timeout and
no retry: a Redis that is slow or gone is reported, never waited for.
"""

import time
from urllib.parse import urlsplit

import redis
import rq

from planetgen.queue import api_jobs, redisclient, redisqueue

TIMEOUT_SECONDS = 1.0
"""float: How long the health check waits for Redis."""


def display_url(url):
    """`url` without a password (`redis://:secret@host:6379/0` becomes `redis://host:6379/0`)."""
    parts = urlsplit(url)
    host = parts.hostname or ""
    if ":" in host:
        host = f"[{host}]"
    netloc = host + (f":{parts.port}" if parts.port else "")
    return f"{parts.scheme}://{netloc}{parts.path}"


def _live_workers(connection):
    """How many RQ workers are registered and alive (their process exists on this host; one on another host
    counts, as all there is to go by is its key)."""
    here = None
    alive = 0
    for worker in rq.Worker.all(connection=connection):
        if here is None:
            import socket
            here = socket.gethostname()
        pid = getattr(worker, "pid", None)
        if worker.hostname and worker.hostname != here or not isinstance(pid, int) or api_jobs._process_exists(pid):
            alive += 1
    return alive


def redis_info(url=None, now=None):
    """
    The Redis server's health.

    Returns:
        dict: `url` (no password), `answering`, `error` (the reason, when it does not answer), and, when it does:
            `version`, `uptime_seconds`, `aof_enabled`, `maxmemory_policy`, `used_memory_bytes`,
            `maxmemory_bytes` (0 for no limit), `last_save_ok`, `last_save_seconds_ago` (`None` before the
            first save) and `workers`.
    """
    url = url or redisqueue.redis_url()
    now = time.time() if now is None else now
    info = {"url": display_url(url), "answering": False, "error": None}
    try:
        client = redisclient.client(url, timeout=TIMEOUT_SECONDS, retries=0)
        server = client.info()
        info.update(
            answering=True,
            version=server.get("redis_version"),
            uptime_seconds=server.get("uptime_in_seconds"),
            aof_enabled=bool(int(server.get("aof_enabled", 0))),
            maxmemory_policy=server.get("maxmemory_policy"),
            used_memory_bytes=server.get("used_memory"),
            maxmemory_bytes=int(server.get("maxmemory", 0)),
            last_save_ok=server.get("rdb_last_bgsave_status") in (None, "ok"),
            last_save_seconds_ago=(max(0, int(now - int(server["rdb_last_save_time"])))
                                   if server.get("rdb_last_save_time") else None),
        )
        try:
            info["workers"] = _live_workers(client)
        except (redis.RedisError, OSError):
            info["workers"] = None
    except (redis.RedisError, OSError) as exc:
        info["error"] = str(exc) or type(exc).__name__
    return info
