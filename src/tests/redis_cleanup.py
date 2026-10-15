# tests/redis_cleanup.py

"""
Clearing what earlier test runs left on the test Redis server (TEST.152).

A test that queues a job leaves the job, its result and its place in the queue's finished, failed or canceled
registry on the server when the test ends. Those keys have no expiry, so a test server that is never emptied
grows by a few hundred keys a run (4,000 had piled up on the one this was written on) and every `scripts/test.py`
run starts with them.

`scripts/test.py` calls `clean` before it starts pytest, on every command, a resumed run included. It only ever
touches the server named by `PLANETGEN_TEST_REDIS_URL` (never `PLANETGEN_REDIS_URL`), and only keys that a test
makes: keys under the `planetgen-test:` prefix, the API test queues (`planetgen-api-test-*`), the run queues
(`planetgen-run-*`, named at random, so a test's cannot be told from any other) and queues called `test-*`, with
the job records those queues list. It does not flush the database, and it leaves the rate-limit keys alone (the
limiter counts of every test expire on their own).

Because a run queue of the application has the same name as a test's, `clean` does nothing while a worker has
sent a heartbeat in the last `WORKER_GRACE_SECONDS`: a real worker may be running jobs on this server.
"""

import datetime
import os

WORKER_GRACE_SECONDS = 90
TEST_KEY_PATTERNS = ("planetgen-test:*", "*planetgen-api-test-*", "rq:*:planetgen-run-*", "rq:*:test-*")
REGISTRIES = ("finished", "failed", "canceled", "started", "deferred", "scheduled")
QUEUE_PREFIX = "rq:queue:"
_BATCH = 500


def _live_worker(client, now):
    """The name of a worker that has sent a heartbeat lately (or whose heartbeat cannot be read), else None."""
    for key in client.smembers("rq:workers"):
        key = key.decode() if isinstance(key, bytes) else key
        beat = client.hget(key, "last_heartbeat")
        if beat is None:
            continue  # the worker's record is gone or empty: a dead one
        beat = beat.decode() if isinstance(beat, bytes) else beat
        try:
            when = datetime.datetime.fromisoformat(beat.replace("Z", "+00:00"))
        except ValueError:
            return key
        if when.tzinfo is None:
            when = when.replace(tzinfo=datetime.timezone.utc)
        if (now - when).total_seconds() < WORKER_GRACE_SECONDS:
            return key
    return None


def _job_ids(client, key):
    kind = client.type(key)
    kind = kind.decode() if isinstance(kind, bytes) else kind
    if kind == "zset":
        found = client.zrange(key, 0, -1)
    elif kind == "list":
        found = client.lrange(key, 0, -1)
    else:
        found = []
    return [i.decode() if isinstance(i, bytes) else i for i in found]


def _delete(client, keys):
    keys = list(keys)
    return sum(client.delete(*keys[start:start + _BATCH]) for start in range(0, len(keys), _BATCH))


def clean(url=None, client=None, now=None):
    """(keys removed, a sentence for the run's output or ""). `url` defaults to `PLANETGEN_TEST_REDIS_URL`; with
    neither, or when the server cannot be reached, nothing happens."""
    url = url if url is not None else os.environ.get("PLANETGEN_TEST_REDIS_URL")
    if client is None:
        if not url:
            return 0, ""
        try:
            import redis
            client = redis.Redis.from_url(url)
            client.ping()
        except Exception:  # noqa: BLE001 -- no redis package or no server: the tests that need one will skip
            return 0, ""
    now = now or datetime.datetime.now(datetime.timezone.utc)
    try:
        worker = _live_worker(client, now)
        if worker:
            return 0, (f"the test Redis has a live worker ({worker}) that sent a heartbeat in the last "
                       f"{WORKER_GRACE_SECONDS} s, so its leftover keys were not cleared")
        keys = set()
        for pattern in TEST_KEY_PATTERNS:
            keys.update(k.decode() if isinstance(k, bytes) else k for k in client.scan_iter(match=pattern, count=1000))
        jobs = set()
        queues = set()
        for key in keys:
            if key.startswith(QUEUE_PREFIX):
                queues.add(key)
            if key.startswith(tuple(f"rq:{name}:" for name in REGISTRIES)) or key.startswith(QUEUE_PREFIX):
                jobs.update(_job_ids(client, key))
        for job in jobs:
            keys.update((f"rq:job:{job}", f"rq:results:{job}"))
        removed = _delete(client, keys)
        if queues:
            client.srem("rq:queues", *queues)
    except Exception as exc:  # noqa: BLE001 -- a server that fails mid-way must not stop the run
        return 0, f"the test Redis could not be cleared ({type(exc).__name__}: {exc})"
    return removed, (f"removed {removed} keys earlier test runs left on the test Redis" if removed else "")
