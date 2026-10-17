# tests/test_xfetch.py

"""
Costly entries that expire by age are fetched again shortly before they expire (PERF.110, `xfetch_early`,
`AgedCache` and `ResponseCache` in `planetgen/web/lib/pagecache.py`): one reader refreshes them, the burst of
readers at the expiry second does not all miss.
"""

import math
import secrets
import threading
import time

import pytest

from planetgen.queue import redisclient
from planetgen.web.lib import xfetch
from planetgen.web.lib import pagecache


class _Clock:
    def __init__(self, now=1000.0):
        self.now = now
        self.lock = threading.Lock()

    def __call__(self):
        return self.now

    def advance(self, seconds):
        with self.lock:
            self.now += seconds


def test_the_refresh_probability_rises_toward_the_expiry():
    # now - cost * beta * ln(rand) >= expires: with rand 0.5 a 2 s entry refreshes from 1.386 s before it expires
    half = lambda: 0.5  # noqa: E731
    assert not xfetch.xfetch_early(1000.0, 1002.0, 2.0, rand=half)
    assert xfetch.xfetch_early(1000.7, 1002.0, 2.0, rand=half)
    # a costlier entry starts earlier, and a larger beta earlier still
    assert xfetch.xfetch_early(1000.7, 1002.0, 4.0, rand=half) and not xfetch.xfetch_early(
        1000.0, 1002.0, 0.5, rand=half)
    assert xfetch.xfetch_early(1000.0, 1002.0, 1.0, beta=3.0, rand=half)
    # the draw decides: rand near 1 is the latest start, rand 0 the earliest
    assert not xfetch.xfetch_early(1000.0, 1002.0, 2.0, rand=lambda: 0.999)
    assert xfetch.xfetch_early(0.0, 1e9, 2.0, rand=lambda: 0.0)
    assert xfetch.xfetch_early(1002.0, 1002.0, 2.0, rand=lambda: 0.999999) and math.isfinite(1.0)


def _many(count, call):
    threads = [threading.Thread(target=call) for _ in range(count)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join(30)


def test_fifty_readers_at_the_expiry_rebuild_a_costly_entry_once_before_it_and_not_after():
    clock = _Clock()
    builds = []

    def build():
        builds.append(clock.now)
        time.sleep(0.05)  # the other readers arrive while it is rebuilt
        clock.advance(2.0)  # it costs 2 s on the entry's clock
        return f"build-{len(builds)}"

    cache = xfetch.AgedCache(60, min_cost=0.5, clock=clock, rand=lambda: 0.5)
    assert cache.get_or_build("stats", build) == "build-1"  # stored at 1002, expires at 1062
    answers = []
    _many(50, lambda: answers.append(cache.get_or_build("stats", build)))
    assert len(builds) == 1 and answers == ["build-1"] * 50  # 1002: far from the expiry

    clock.now = 1061.0  # inside the refresh window (from 1060.6), before the expiry
    answers.clear()
    _many(50, lambda: answers.append(cache.get_or_build("stats", build)))
    assert len(builds) == 2
    assert set(answers) <= {"build-1", "build-2"}  # nobody waited for the rebuild or built again
    assert cache.get_or_build("stats", build) == "build-2"  # stored at 1063, expires at 1123

    clock.now = 1062.5  # the second after the first copy would have expired
    answers.clear()
    _many(50, lambda: answers.append(cache.get_or_build("stats", build)))
    assert len(builds) == 2 and answers == ["build-2"] * 50


def test_an_entry_that_expired_is_rebuilt_by_one_reader_while_the_others_wait():
    clock = _Clock()
    builds = []

    def build():
        builds.append(1)
        time.sleep(0.1)
        clock.advance(2.0)
        return f"build-{len(builds)}"

    cache = xfetch.AgedCache(60, min_cost=0.5, clock=clock, rand=lambda: 0.999)
    cache.get_or_build("k", build)
    clock.now = 1070.0
    answers = []
    _many(20, lambda: answers.append(cache.get_or_build("k", build)))
    assert len(builds) == 2 and answers == ["build-2"] * 20


def test_a_build_faster_than_the_stated_time_is_never_kept():
    clock = _Clock()
    builds = []

    def build():
        builds.append(1)
        clock.advance(0.01)
        return len(builds)

    cache = xfetch.AgedCache(60, min_cost=0.5, clock=clock)
    assert [cache.get_or_build("k", build) for _ in range(3)] == [1, 2, 3]


def test_a_failed_rebuild_keeps_the_old_copy():
    clock = _Clock()
    calls = []

    def build():
        calls.append(1)
        clock.advance(2.0)
        if len(calls) > 1:
            raise OSError("the API is down")
        return "kept"

    cache = xfetch.AgedCache(60, min_cost=0.5, clock=clock, rand=lambda: 0.5)
    cache.get_or_build("k", build)
    clock.now = 1061.0
    with pytest.raises(OSError):
        cache.get_or_build("k", build)
    clock.now = 1061.2
    # the entry is still there and not expired: the next reader tries again, the failure is not cached
    with pytest.raises(OSError):
        cache.get_or_build("k", build)
    clock.now = 1001.0
    assert cache.get_or_build("k", build) == "kept"


@pytest.fixture
def db(redis_server):
    name = f"xfetch_{secrets.token_hex(6)}"
    yield name
    client = redisclient.client(decode_responses=True)
    for key in client.scan_iter(f"planetgen:*:{name}*"):
        client.delete(key)


def _response_cache(clock, **settings):
    """A cache that does not take a copy back from Redis (that copy would stand in for the one that aged out)."""
    cache = pagecache.ResponseCache(lambda db, since: {"stamp": "a", "state": "state-a"},
                                   dict(pagecache.DEFAULTS, stamp_seconds=1e9, max_age_seconds=300, **settings),
                                   clock=clock, rand=lambda: 0.5)
    cache._recall = lambda db, target: None
    return cache


def test_a_costly_page_answer_is_fetched_again_by_one_reader_before_it_ages_out(db):
    clock = _Clock()
    cache = _response_cache(clock)
    assert cache.get(db, "/search?x=1") is None
    token = cache.begin(db)
    clock.advance(2.0)  # the fetch took 2 s
    cache.put(db, "/search?x=1", "body", token)  # stored at 1002, ages out at 1302
    clock.now = 1300.0
    assert cache.get(db, "/search?x=1") == "body"  # 1300.6 is the start of the window
    clock.now = 1301.0
    assert cache.get(db, "/search?x=1") is None  # the first reader inside the window fetches it again
    assert cache.get(db, "/search?x=1") == "body"  # the others keep getting the old copy
    assert cache.get(db, "/search?x=1") == "body"
    token = cache.begin(db)
    clock.advance(2.0)
    cache.put(db, "/search?x=1", "body-2", token)  # stored at 1303, ages out at 1603
    assert cache.get(db, "/search?x=1") == "body-2"


def test_a_cheap_page_answer_just_ages_out(db):
    clock = _Clock()
    cache = _response_cache(clock)
    token = cache.begin(db)
    clock.advance(0.05)
    cache.get(db, "/sectors/1")
    cache.put(db, "/sectors/1", "body", token)
    clock.now = 1299.9
    assert cache.get(db, "/sectors/1") == "body"
    assert cache.get(db, "/sectors/1") == "body"
    clock.now = 1301.0
    assert cache.get(db, "/sectors/1") is None


def test_an_early_refresh_that_never_arrives_is_claimed_again_later(db):
    clock = _Clock()
    cache = _response_cache(clock)
    cache.get(db, "/search")
    token = cache.begin(db)
    clock.advance(2.0)
    cache.put(db, "/search", "body", token)
    clock.now = 1301.0
    assert cache.get(db, "/search") is None  # claimed by a reader whose fetch then fails
    assert cache.get(db, "/search") == "body"
    clock.now = 1301.0 + 2 * 2.0 + 1.5
    assert cache.get(db, "/search") is None  # the claim lapsed: another reader tries
