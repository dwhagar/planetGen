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

Two runs on one machine would still take the same databases (worker 0 of each is database 1). TEST.724: a run started with a
slot (`PLANETGEN_TEST_REDIS_SLOT`, `scripts/test.py --slot N`, `scripts/train_gate.py run --slot N`) takes the `SLOT_WIDTH` databases of its slot
only, so up to `SLOTS` runs at once share a stock Redis and meet nowhere; a slot's run has at most `SLOT_WIDTH` workers. A run with no slot
keeps all fifteen, as before, and is the machine's only run.
"""

import os
import re
from urllib.parse import urlsplit, urlunsplit

DATABASES = 16
"""int: How many databases a stock Redis server has (`databases 16`); worker numbers wrap around it."""

SLOT_ENVIRONMENT = "PLANETGEN_TEST_REDIS_SLOT"
"""str: The environment variable that gives a run its slot (a number from 0)."""

SLOT_WIDTH = 3
"""int: How many databases, and so how many test workers, one slot has (the gate's cap of three workers, TEST.724)."""

SLOTS = (DATABASES - 1) // SLOT_WIDTH
"""int: How many slots a stock Redis server holds: database 0 is the machine's own."""


def slot_databases(slot):
    """The `SLOT_WIDTH` databases of slot number `slot`: slot 0 is 1 to 3, slot 1 is 4 to 6, and so on."""
    if not 0 <= slot < SLOTS:
        raise ValueError(f"a slot is a number from 0 to {SLOTS - 1}, not {slot}")
    return list(range(1 + slot * SLOT_WIDTH, 1 + (slot + 1) * SLOT_WIDTH))


def slot_of(environ=None):
    """The slot the environment names, or `None` for a run with none (it keeps every database)."""
    value = (os.environ if environ is None else environ).get(SLOT_ENVIRONMENT, "").strip()
    return int(value) if value.isdigit() else None


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
    slot = slot_of()
    if slot is not None:
        mine = slot_databases(slot)
        database = mine[number % SLOT_WIDTH]
        if database == own:      # the server's own database inside the slot: the next one of the slot (the workers wrap, so two may share)
            database = mine[(mine.index(database) + 1) % SLOT_WIDTH]
        return urlunsplit(parts._replace(path=f"/{database}"))
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
