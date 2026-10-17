# planetgen/util/xfetch.py

"""
Fetching a costly, age-expiring value again shortly before it expires (PERF.110).

When many readers want a value that expires by age and took a while to build (the admin Stats page's database
figures, a search's facets), they all find it gone in the same second and all rebuild it. XFetch (Vattani,
Chierichetti and Lowenstein, "Optimal Probabilistic Cache Stampede Prevention", PVLDB 8(8), 2015) has each
reader decide on its own, with a probability that rises toward the expiry and starts earlier the longer the
value took to build, whether to rebuild it now; with one reader rebuilding at a time the rest keep getting
the old copy until the new one is in. `beta` is 1 as in the paper.

`xfetch_early` is the decision, with an injectable random source. `AgedCache` is a small in-process cache of
such values; the page cache (`planetgen/web/lib/pagecache.py`) uses `xfetch_early` for its own entries.
"""

import math
import random
import threading
import time

XFETCH_MIN_SECONDS = 0.5
"""float: An answer that took this long to fetch (or more) is refreshed before it expires (`xfetch_early`); a
cheaper one is just fetched again when it has aged out."""

XFETCH_BETA = 1.0
"""float: XFetch's beta (Vattani, Chierichetti and Lowenstein, PVLDB 8(8), 2015): above 1 refreshes earlier."""



def xfetch_early(now, expires_at, cost, beta=XFETCH_BETA, rand=random.random):
    """
    Whether an entry that expires at `expires_at` should be rebuilt at `now`, ahead of its expiry (XFetch,
    PERF.110): the probability rises toward the expiry, and the longer the entry took to build (`cost`
    seconds) the earlier it starts, so one reader rebuilds it while the rest still get the old copy.

    Args:
        now (float): Seconds on the same clock as `expires_at`.
        expires_at (float): When the entry expires.
        cost (float): How long the entry took to build, in seconds.
        beta (float): XFetch's beta; 1 is the paper's default.
        rand (callable): `rand()` -> a float in [0, 1), replaceable in tests.

    Returns:
        bool: True when `now - cost * beta * ln(rand()) >= expires_at`; always True once expired.
    """
    draw = rand()
    if draw <= 0.0:
        return True  # ln(0) is minus infinity: the earliest possible refresh
    return now - cost * beta * math.log(draw) >= expires_at


class AgedCache:
    """
    A small in-process cache for one costly, age-expiring value per key that no write invalidates (the admin
    Stats page's database figures, PERF.110), with XFetch: an entry that took `min_cost` seconds or more to
    build is kept `ttl` seconds and rebuilt by one reader shortly before then while the others keep getting
    the old copy; one that built faster is not kept at all, so a small galaxy always shows live figures.

    Args:
        ttl (float): Seconds a kept value is served.
        min_cost (float): The build time from which a value is kept.
        beta (float): XFetch's beta.
        clock (callable): Seconds, monotonic; replaceable in tests.
        rand (callable): `rand()` -> a float in [0, 1); replaceable in tests.
    """

    def __init__(self, ttl, min_cost=XFETCH_MIN_SECONDS, beta=XFETCH_BETA, clock=time.monotonic, rand=random.random):
        self._ttl = float(ttl)
        self._min_cost = float(min_cost)
        self._beta = float(beta)
        self._clock = clock
        self._rand = rand
        self._lock = threading.Lock()
        self._entries = {}  # key -> (value, stored_at, cost)
        self._builds = {}   # key -> the lock held while that key is built

    def _fresh(self, key, now):
        with self._lock:
            entry = self._entries.get(key)
        return entry if entry is not None and now < entry[1] + self._ttl else None

    def _early(self, entry, now):
        return xfetch_early(now, entry[1] + self._ttl, entry[2], self._beta, self._rand)

    def clear(self):
        """Forgets every value."""
        with self._lock:
            self._entries.clear()

    def get_or_build(self, key, build):
        """
        The value for `key`: the kept one, or `build()`'s result (an exception from `build` propagates and keeps
        the old copy). Concurrent callers share one build: when an entry is close to expiry only the reader that
        wins its lock rebuilds it, and once it has expired the others wait for that build.
        """
        entry = self._fresh(key, self._clock())
        if entry is not None and not self._early(entry, self._clock()):
            return entry[0]
        with self._lock:
            build_lock = self._builds.setdefault(key, threading.Lock())
        if entry is not None:
            if not build_lock.acquire(blocking=False):
                return entry[0]  # somebody is already rebuilding it: keep serving the old copy
        else:
            build_lock.acquire()
        try:
            again = self._fresh(key, self._clock())
            if again is not None and (again is not entry or not self._early(again, self._clock())):
                return again[0]  # rebuilt while this caller waited
            started = self._clock()
            value = build()
            cost = self._clock() - started
            with self._lock:
                if cost >= self._min_cost:
                    self._entries[key] = (value, started + cost, cost)
                else:
                    self._entries.pop(key, None)
            return value
        finally:
            build_lock.release()


def is_cacheable(path):
    """Whether a GET of `path` (under `/api`, no query) may be cached."""
    return not path.startswith(UNCACHED_PREFIXES)


def _path_of(target):
    return target.partition("?")[0]


def shared_key(db, target):
    """The Redis key of the answer to `target` in `db`."""
    return f"planetgen:page:{db}:{hashlib.sha256(target.encode('utf-8')).hexdigest()[:32]}"
