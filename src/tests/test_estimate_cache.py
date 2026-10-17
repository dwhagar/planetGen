# tests/test_estimate_cache.py

"""
The Generate page's cached estimate (PERF.109, `planetgen/web/lib/estimatecache.py`): a repeat answers from
Redis, a change of the versions it was made under misses, concurrent identical requests run it once, and a
request that cannot read the versions simply runs it.
"""

import threading
import time
import uuid

import pytest
import redis

from planetgen.web.lib import estimatecache, singleflight

ARGV = ["planetgen", "galaxy", "--shell", "40"]


@pytest.fixture
def database(redis_server):
    name = f"estimate_cache_{uuid.uuid4().hex[:10]}"
    yield name
    client = redis.Redis.from_url(redis_server)
    for key in client.scan_iter(match=f"planetgen:estimate:{name}:*"):
        client.delete(key)


class _Counter:
    def __init__(self, seconds=0.0):
        self.calls = 0
        self.seconds = seconds
        self.lock = threading.Lock()

    def __call__(self):
        with self.lock:
            self.calls += 1
            number = self.calls
        time.sleep(self.seconds)
        return {"sectors": 9, "run": number}


VERSIONS = {"generator_epoch": 10, "data": ["a", {"_n": 1, "epoch": 0}], "speeds": "2@t/1@t"}


def test_a_second_identical_estimate_answers_from_the_cache(database):
    compute = _Counter()
    first = estimatecache.estimate(database, "galaxy", ARGV, lambda: VERSIONS, compute)
    started = time.monotonic()
    second = estimatecache.estimate(database, "galaxy", ARGV, lambda: VERSIONS, compute)
    assert compute.calls == 1 and second == first == {"sectors": 9, "run": 1}
    assert time.monotonic() - started < 0.2


def test_a_changed_plan_epoch_data_or_statistics_misses(database):
    compute = _Counter()
    estimate = lambda argv, versions, kind="galaxy": estimatecache.estimate(  # noqa: E731
        database, kind, argv, lambda: versions, compute)
    estimate(ARGV, VERSIONS)
    estimate(ARGV + ["--shell-band", "1"], VERSIONS)
    estimate(ARGV, VERSIONS, kind="check_db")
    estimate(ARGV, dict(VERSIONS, generator_epoch=11))
    estimate(ARGV, dict(VERSIONS, data=["a", {"_n": 2, "epoch": 0}]))
    estimate(ARGV, dict(VERSIONS, data=["b", {"_n": 1, "epoch": 0}]))
    estimate(ARGV, dict(VERSIONS, speeds="3@t/1@t"))
    assert compute.calls == 7
    estimate(ARGV, dict(VERSIONS, speeds="3@t/1@t"))
    assert compute.calls == 7


def test_two_concurrent_identical_requests_run_the_estimate_once(database):
    compute = _Counter(seconds=0.6)
    results = []

    def ask():
        results.append(estimatecache.estimate(database, "galaxy", ARGV, lambda: VERSIONS, compute))

    threads = [threading.Thread(target=ask) for _ in range(4)]
    for thread in threads:
        thread.start()
        time.sleep(0.05)
    for thread in threads:
        thread.join(30)
    assert compute.calls == 1
    assert results == [{"sectors": 9, "run": 1}] * 4


def test_an_estimate_that_fails_is_not_stored_and_frees_the_claim(database):
    def failing():
        raise ValueError("Shell 900 holds 5000 sectors")

    with pytest.raises(ValueError):
        estimatecache.estimate(database, "galaxy", ARGV, lambda: VERSIONS, failing)
    compute = _Counter()
    assert estimatecache.estimate(database, "galaxy", ARGV, lambda: VERSIONS, compute)["run"] == 1


def test_a_claimant_that_never_delivers_is_built_after_the_wait(database, monkeypatch):
    monkeypatch.setattr(estimatecache, "WAIT_SECONDS", 0.3)
    key = estimatecache.cache_key(database, "galaxy", ARGV)
    owned, _held = singleflight.claim([key])  # claimed by somebody who then dies
    try:
        compute = _Counter()
        assert estimatecache.estimate(database, "galaxy", ARGV, lambda: VERSIONS, compute)["run"] == 1
    finally:
        singleflight.release(owned)


def test_versions_that_cannot_be_read_run_the_estimate_uncached(database):
    def down():
        raise OSError("no redis")

    compute = _Counter()
    estimatecache.estimate(database, "galaxy", ARGV, down, compute)
    estimatecache.estimate(database, "galaxy", ARGV, down, compute)
    assert compute.calls == 2


def test_the_stored_estimate_is_dropped_after_its_lifetime(database, redis_server):
    compute = _Counter()
    estimatecache.estimate(database, "galaxy", ARGV, lambda: VERSIONS, compute)
    key = estimatecache.cache_key(database, "galaxy", ARGV)
    ttl = redis.Redis.from_url(redis_server).ttl(key)
    assert 0 < ttl <= estimatecache.KEEP_SECONDS
