# tests/private_redis.py

"""
A Redis database of the test worker's own, for every test it runs (TEST.266, TEST.280).

A test that starts `planetgen` as a process, or real queue workers, runs it on the machine's Redis with the default
key prefix, which names the one live run of the machine (`planetgen:run:active`) and the queues the update guard looks
at (`planetgen-bulk`); the in-process tests get a prefix of their own (`inline_queue`), but a process cannot be given
one. Two such tests on different xdist workers then refused each other as "Busy", and `update_guard` and the update
step refused to begin because another worker's job was running. Each worker takes the Redis database numbered after
it (`tests/conftest.py` points `PLANETGEN_REDIS_URL` at it before any test runs and empties it before and after the
session), so their runs never meet.
"""

import os
import re
from urllib.parse import urlsplit, urlunsplit

DATABASES = 16
"""int: How many databases a stock Redis server has (`databases 16`); worker numbers wrap around it."""


def url_for_worker(url, worker=None):
    """
    `url` with its database changed to one private to the xdist worker, never the database `url` itself names.

    Args:
        url (str): The test server's Redis URL.
        worker (str, optional): The xdist worker id (`gw3`); by default this process's own (`master` without xdist).
    """
    worker = worker or os.environ.get("PYTEST_XDIST_WORKER", "master")
    found = re.search(r"\d+", worker)
    number = int(found.group()) if found else 0
    parts = urlsplit(url)
    own = int(parts.path.strip("/") or 0)
    database = 1 + number % (DATABASES - 1)
    if database == own:
        database = 1 + database % (DATABASES - 1)
    return urlunsplit(parts._replace(path=f"/{database}"))


def flush(url):
    """Empties the database `url` names; nothing happens without a server that answers."""
    try:
        import redis

        redis.Redis.from_url(url, socket_timeout=5, socket_connect_timeout=5).flushdb()
    except Exception:  # noqa: BLE001 -- no redis package or no server: the tests that need one skip
        pass
