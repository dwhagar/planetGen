# html/lib/pagecache.py

"""
An in-memory cache of the API's public read responses, so the pages
don't query the database on every request (TODO 8).

Every page reads through `apiclient`'s plain GET wrappers (`get_sector`,
`get_system`, `get_phenomena`, ...), which send no cookie, so their
answers are the same for every visitor and can be shared. `apiclient`
asks the cache installed with `apiclient.set_response_cache` before each
such call, and stores the raw JSON body after a miss (each hit is parsed
afresh, so a page that edits the dict it gets never changes the cached
copy).

A cached answer is dropped when any of these says it may be stale:

- **A write through the API.** Any successful `POST`/`PATCH`/`PUT`/
  `DELETE` under `/api` in this process (a page's admin form, or an API
  client) clears the cache (`web.init_app` hooks this up), so an edit
  shows on the very next page.
- **The galaxy's content stamp.** At most once per `stamp_seconds` per
  database, one cheap call (`GET /api/galaxy/changes`, the stamp the
  Galaxy Map's tile cache uses, from the v27 `modified_at` columns) says
  whether sectors or systems were added, edited or deleted, by another
  process too (a `generate.py` job, another web worker). A new stamp
  drops that database's entries.
- **Age.** Nothing is served older than `max_age_seconds`, the backstop
  for edits the stamp can't see (a rename of a single system from the
  command line, say).

Every failure fails open: if the stamp can't be read the call simply
isn't cached.

Settings live under `page_cache` in `config.json` (`enabled`,
`max_entries`, `max_mb`, `stamp_seconds`, `max_age_seconds`), and
`PLANETGEN_PAGE_CACHE=off` turns it off.
"""

import collections
import os
import threading
import time

DEFAULTS = {
    "enabled": True,
    "max_entries": 2000,
    "max_mb": 64,
    "stamp_seconds": 15,
    "max_age_seconds": 300,
}
"""dict: The settings used when `config.json` doesn't give them."""

UNCACHED_PREFIXES = ("/galaxy/changes", "/galaxy/tiles", "/galaxy/stage", "/galaxy/stamp", "/auth", "/admin")
"""tuple: Paths never cached here: the stamp itself, the Galaxy Map's tile
and stage data (`tilecache.py` keeps those on disk), and anything that
depends on who is asking."""

MAX_BODY_FRACTION = 0.25
"""float: One body may use at most this share of `max_mb`, so a single
huge answer can't push everything else out."""


def settings_from(config, environ=None):
    """The cache's settings: `DEFAULTS`, then `config["page_cache"]`,
    then `PLANETGEN_PAGE_CACHE=off`."""
    environ = os.environ if environ is None else environ
    merged = dict(DEFAULTS)
    merged.update((config or {}).get("page_cache") or {})
    if str(environ.get("PLANETGEN_PAGE_CACHE", "")).strip().lower() in ("0", "off", "false", "no"):
        merged["enabled"] = False
    return merged


def is_cacheable(path):
    """Whether a GET of `path` (under `/api`, no query) may be cached."""
    return not path.startswith(UNCACHED_PREFIXES)


class ResponseCache:
    """
    The cache itself: raw response bodies keyed by request target, least
    recently used first out, grouped by database for the stamp check.

    Args:
        stamp_for (callable): `stamp_for(db)` -> the database's current
            content stamp (a short string). Raising anything means "can't
            tell", and that call goes uncached.
        settings (dict): `settings_from(...)`'s keys.
        clock (callable): Seconds, monotonic; replaceable in tests.
    """

    def __init__(self, stamp_for, settings=None, clock=time.monotonic):
        settings = dict(DEFAULTS, **(settings or {}))
        self._stamp_for = stamp_for
        self._max_entries = int(settings["max_entries"])
        self._max_bytes = int(float(settings["max_mb"]) * 1024 * 1024)
        self._stamp_seconds = float(settings["stamp_seconds"])
        self._max_age = float(settings["max_age_seconds"])
        self._clock = clock
        self._lock = threading.Lock()
        self._entries = collections.OrderedDict()  # (db, target) -> (body, stored_at)
        self._bytes = 0
        self._stamps = {}  # db -> (stamp, checked_at)
        self._generation = 0  # bumped by clear()

    def __len__(self):
        return len(self._entries)

    @property
    def generation(self):
        """Read before fetching a miss and hand to `put`, so an answer
        fetched while an edit cleared the cache is never stored."""
        return self._generation

    def clear(self):
        """Drops everything (an edit went through the API)."""
        with self._lock:
            self._generation += 1
            self._entries.clear()
            self._bytes = 0
            self._stamps.clear()

    def _fresh_stamp(self, db):
        """True when `db`'s entries are still good, re-checking its stamp
        when the last check is older than `stamp_seconds`. False when the
        stamp can't be read (the call then goes uncached)."""
        now = self._clock()
        with self._lock:
            known = self._stamps.get(db)
        if known is not None and now - known[1] < self._stamp_seconds:
            return True
        if db is None:
            return True  # Nothing database-specific (the list of databases): age alone.
        try:
            stamp = self._stamp_for(db)
        except Exception:  # noqa: BLE001 -- fail open: just don't cache
            return False
        with self._lock:
            if known is not None and known[0] != stamp:
                self._drop_db(db)
            self._stamps[db] = (stamp, now)
        return True

    def _drop_db(self, db):
        for key in [key for key in self._entries if key[0] == db]:
            body, _ = self._entries.pop(key)
            self._bytes -= len(body)

    def get(self, db, target):
        """The cached body for `target` in `db`, or `None`."""
        if not self._fresh_stamp(db):
            return None
        key = (db, target)
        with self._lock:
            entry = self._entries.get(key)
            if entry is None:
                return None
            body, stored_at = entry
            if self._clock() - stored_at >= self._max_age:
                del self._entries[key]
                self._bytes -= len(body)
                return None
            self._entries.move_to_end(key)
            return body

    def put(self, db, target, body, generation):
        """Stores `body` for `target` in `db` -- skipped when it's too
        big, when `db`'s stamp couldn't be checked, or when the cache was
        cleared since `generation` was read."""
        if len(body) > self._max_bytes * MAX_BODY_FRACTION:
            return
        with self._lock:
            if generation != self._generation or (db is not None and db not in self._stamps):
                return
            key = (db, target)
            old = self._entries.pop(key, None)
            if old is not None:
                self._bytes -= len(old[0])
            self._entries[key] = (body, self._clock())
            self._bytes += len(body)
            while self._entries and (len(self._entries) > self._max_entries or self._bytes > self._max_bytes):
                _, (evicted, _) = self._entries.popitem(last=False)
                self._bytes -= len(evicted)
