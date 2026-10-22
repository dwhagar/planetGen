# tests/private_redis.py

"""
A Redis database of the test worker's own for a test that starts a process of the machine's runs (TEST.266).

A test that starts `planetgen` as a process runs it on the machine's Redis with the default key prefix, which names the
one live run of the machine (`planetgen:run:active`); the in-process tests get a prefix of their own (`inline_queue`),
but a process cannot be given one. Two such tests on different xdist workers then refused each other as "Busy". Each
worker takes the Redis database numbered after it, so their runs never meet.
"""

import os
import re
from urllib.parse import urlsplit, urlunsplit

import pytest

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


@pytest.fixture
def private_redis(monkeypatch):
    """
    The test, and any process it starts, uses a Redis database of the worker's own, emptied before the test and after
    it. Does nothing without a Redis server. Import it into the test module to use it.
    """
    url = os.environ.get("PLANETGEN_REDIS_URL")
    if not url:
        yield None
        return
    import redis

    private = url_for_worker(url)
    client = redis.Redis.from_url(private)
    client.flushdb()
    monkeypatch.setenv("PLANETGEN_REDIS_URL", private)
    try:
        yield private
    finally:
        client.flushdb()
