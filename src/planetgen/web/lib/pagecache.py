# planetgen/web/lib/pagecache.py

"""
A cache of the API's public read responses, so the pages don't query the database on every
request (PERF.2), shared by every web process and the warm-up through Redis (PERF.104;
docs/design/cache-warm-up.md).

Every page reads through `apiclient`'s plain GET wrappers (`get_sector`, `get_system`,
`get_phenomena`, ...), which send no cookie, so their answers are the same for every visitor and
can be shared. `apiclient` asks the cache installed with `apiclient.set_response_cache` before
each such call, and stores the raw JSON body after a miss (each hit is parsed afresh, so a page
that edits the dict it gets never changes the cached copy).

**An answer is served only while the data it read has not changed.** It is stored with the tags
of the data it read (`pagetags.read_tags`) and the counters those tags had
(`planetgen.db.dataversion`, in Redis); a copy whose counters have moved on, or whose Redis
incarnation is not the current one, is never served. So a stale copy is never served as current,
by this process or another, and a change rebuilds only the answers that read it. Whatever changes
data bumps the tags it touches:

- **A write through the API** (`web.init_app` hooks this up): `pagetags.change_for_write` names
  the change and `dataversion.record` bumps its tags before the response goes out, so an edit
  shows on the very next page.
- **The command line** (the population pass, the orbital update, generation): they record their
  change when they finish.
- **The galaxy's content stamp**, the backstop for everything else (another process, a direct
  database edit). At most once per `stamp_seconds` per database, one cheap call (`GET
  /api/galaxy/changes`, the stamp the Galaxy Map's tile cache uses, from the v27 `modified_at`
  columns) says whether sectors or systems were added, edited or deleted. A new stamp bumps the
  tags of the sectors named (`sector:<id>`, `sector-data`, and the set tags when something was
  added) or, when the change cannot be placed, `epoch`. One process wins the claim on each state
  (`dataversion.claim_state`) and bumps; the others skip. While a fill changes the galaxy faster
  than any page could be rebuilt (`busy`, PERF.38) nothing is bumped and the entries are kept until
  they age out, the one case where an answer can be older than the data (bounded by
  `max_age_seconds`).
- **Age.** Nothing is served older than `max_age_seconds` from this process's copy, the backstop for
  edits the stamp can't see (a rename of a single system from the command line, say).

Every failure fails open: if Redis or the stamp can't be read the call simply isn't cached.

A process keeps the answers it has served in a `cachetools.TTLCache` sized in characters of body
(least recently used out first), and Redis (key `planetgen:page:<database>:<digest of the target>`,
kept `shared_seconds`) holds the copy every process and the warm-up share. cachetools isn't
thread-safe, so every use is under the cache's lock.

Settings live under `page_cache` in `config.json` (`enabled`, `max_entries`, `max_mb`,
`stamp_seconds`, `max_age_seconds`, `shared_seconds`), and `PLANETGEN_PAGE_CACHE=off` turns it off
(all in `planetgen.util.settings`).
"""

import hashlib
import json
import threading
import time

import cachetools
import redis

from planetgen.db import dataversion
from planetgen.queue import redisclient
from planetgen.util.settings import PageCache
from planetgen.web.lib import pagetags

DEFAULTS = PageCache().model_dump()
"""dict: The model's defaults (`settings.PageCache`), the settings a cache
built without any uses."""

UNCACHED_PREFIXES = ("/galaxy/changes", "/galaxy/tiles", "/galaxy/stage", "/galaxy/window", "/galaxy/stamp", "/auth", "/admin")
"""tuple: Paths never cached here: the stamp itself, the Galaxy Map's tile
and stage data (`tilecache.py` keeps those on disk), and anything that
depends on who is asking."""

MAX_BODY_FRACTION = 0.25
"""float: One body may use at most this share of `max_mb`, so a single
huge answer can't push everything else out."""

SOCKET_SECONDS = 0.5
"""float: How long one Redis call may take before the answer goes uncached."""


def is_cacheable(path):
    """Whether a GET of `path` (under `/api`, no query) may be cached."""
    return not path.startswith(UNCACHED_PREFIXES)


def _path_of(target):
    return target.partition("?")[0]


def shared_key(db, target):
    """The Redis key of the answer to `target` in `db`."""
    return f"planetgen:page:{db}:{hashlib.sha256(target.encode('utf-8')).hexdigest()[:32]}"


class ResponseCache:
    """
    The cache itself: raw response bodies keyed by request target, least
    recently used first out, grouped by database for the stamp check, each
    valid only while the data it read is unchanged.

    Args:
        changes_for (callable): `changes_for(db, since)` -> what
            `GET /api/galaxy/changes` answers: `stamp` (the database's
            current content stamp, a short string), `state` (to pass as
            `since` next time) and, when a fill is changing too much at
            once, `busy`. Raising anything means "can't tell", and that
            call goes uncached.
        settings (dict): `settings.PageCache`'s fields as a dict.
        clock (callable): Seconds, monotonic; replaceable in tests.
    """

    def __init__(self, changes_for, settings=None, clock=time.monotonic):
        settings = dict(DEFAULTS, **(settings or {}))
        self._changes_for = changes_for
        self._max_entries = int(settings["max_entries"])
        self._max_bytes = int(float(settings["max_mb"]) * 1024 * 1024)
        self._stamp_seconds = float(settings["stamp_seconds"])
        self._max_age = float(settings["max_age_seconds"])
        self._shared_seconds = int(float(settings["shared_seconds"]))
        self._clock = clock
        self._lock = threading.Lock()
        # (db, target) -> (incarnation, {tag: counter}, body); the timer reads self._clock each time, so
        # a test can swap the clock after construction.
        self._entries = cachetools.TTLCache(maxsize=max(self._max_bytes, 1), ttl=self._max_age,
                                            timer=lambda: self._clock(), getsizeof=lambda entry: len(entry[2]))
        self._stamps = {}  # db -> (stamp, checked_at, state, accepted_at, busy)
        self._generation = 0  # bumped by clear()

    def __len__(self):
        with self._lock:
            self._entries.expire()
            return len(self._entries)

    def begin(self, db):
        """
        Read before fetching a miss and hand to `put`, so an answer fetched while data changed (or the
        cache was cleared) is never stored.

        Returns:
            tuple or None: `None` when the versions can't be read (the answer then goes uncached).
        """
        generation = self._generation
        if db is None:
            return generation, None, None
        try:
            incarnation, found = dataversion.versions(db, [dataversion.COUNT_FIELD])
        except dataversion.Unavailable:
            return None
        return generation, incarnation, found[dataversion.COUNT_FIELD]

    def clear(self):
        """Drops this process's copies (a test, or a write nothing places)."""
        with self._lock:
            self._generation += 1
            self._entries.clear()
            self._stamps.clear()

    def _redis(self):
        return redisclient.client(timeout=SOCKET_SECONDS, retries=1)

    def _bump_from(self, db, changes, previous_state):
        """Turns what `changes` says into tag bumps, once across the web processes."""
        state = changes.get("state") or changes.get("stamp")
        try:
            if not dataversion.claim_state(db, state, previous_state):
                return  # another process has bumped for this state already
        except dataversion.Unavailable:
            return
        if changes.get("full") or "sectors" not in changes:
            dataversion.bump(db, dataversion.TAG_EPOCH)
            return
        tags = [dataversion.sector_tag(printed) for printed in changes["sectors"]]
        if changes["sectors"]:
            tags.append(dataversion.TAG_SECTOR_DATA)
        if changes.get("grew"):
            tags += [*dataversion.SET_TAGS, *dataversion.DATA_TAGS]
        dataversion.bump(db, *tags)

    def _fresh_stamp(self, db):
        """True when `db`'s entries can be checked, re-reading the galaxy's state when the last check is
        older than `stamp_seconds` (a change found bumps the tags it touched). False when it can't be
        read (the call then goes uncached)."""
        now = self._clock()
        with self._lock:
            known = self._stamps.get(db)
        if known is not None and now - known[1] < self._stamp_seconds:
            return True
        if db is None:
            return True  # Nothing database-specific (the list of databases): age alone.
        try:
            recorded = dataversion.recorded_state(db)
            changes = self._changes_for(db, recorded)
            stamp = changes["stamp"]
            state = changes.get("state") or stamp
            busy = bool(changes.get("busy"))
            if recorded is None:
                # No process has looked at this database since Redis started: nothing cached is older than
                # a state anyone saw, so this state is simply the first.
                dataversion.claim_state(db, state, None)
            elif state != recorded:
                if busy and known is not None and known[4] and now - known[3] < self._max_age:
                    # A fill is adding sectors faster than pages could be rebuilt (PERF.38) -- this is the second
                    # check in a row that says so: keep the entries (they age out on their own) and the state
                    # they were accepted under, so the check that follows the fill still sees everything since;
                    # just look again later.
                    with self._lock:
                        self._stamps[db] = (known[0], now, known[2], known[3], True)
                    return True
                self._bump_from(db, changes, recorded)
        except Exception:  # noqa: BLE001 -- fail open: just don't cache
            return False
        with self._lock:
            self._stamps[db] = (stamp, now, state, now, busy)
        return True

    def _valid(self, db, entry):
        """Whether `entry`'s tag counters still stand."""
        incarnation, tagged, _body = entry
        try:
            current_id, current = dataversion.versions(db, tagged)
        except dataversion.Unavailable:
            return False
        return current_id == incarnation and current == tagged

    def get(self, db, target):
        """The cached body for `target` in `db`, or `None`."""
        if not self._fresh_stamp(db):
            return None
        with self._lock:
            entry = self._entries.get((db, target))
        if db is None:
            return None if entry is None else entry[2]
        if entry is None:
            entry = self._recall(db, target)
            if entry is None:
                return None
        if not self._valid(db, entry):
            with self._lock:
                self._entries.pop((db, target), None)
            return None
        return entry[2]

    def _recall(self, db, target):
        """The copy Redis holds (made by another process or the warm-up), taken into this process."""
        try:
            raw = self._redis().hmget(shared_key(db, target), "i", "v", "b")
        except (redis.RedisError, OSError):
            return None
        if raw[0] is None or raw[1] is None or raw[2] is None:
            return None
        try:
            entry = (raw[0].decode("utf-8"), json.loads(raw[1]), raw[2].decode("utf-8"))
        except (ValueError, UnicodeDecodeError):
            return None
        with self._lock:
            self._entries[(db, target)] = entry
            while len(self._entries) > self._max_entries:
                self._entries.popitem()
        return entry

    def put(self, db, target, body, token):
        """Stores `body` for `target` in `db` -- skipped when it's too
        big, when `db`'s stamp couldn't be checked, when data changed (or the cache was cleared) since
        `begin` was called, or when the versions can't be read."""
        if len(body) > self._max_bytes * MAX_BODY_FRACTION or token is None:
            return
        generation, incarnation, count = token
        if db is None:
            with self._lock:
                if generation == self._generation:
                    self._entries[(db, target)] = ("", {}, body)
            return
        tags = pagetags.read_tags(_path_of(target), body)
        try:
            now_id, found = dataversion.versions(db, [dataversion.COUNT_FIELD, *tags])
        except dataversion.Unavailable:
            return
        if now_id != incarnation or found.pop(dataversion.COUNT_FIELD) != count:
            return  # data changed while the answer was fetched: it may be older than the data
        entry = (now_id, found, body)
        with self._lock:
            if generation != self._generation or db not in self._stamps:
                return
            self._entries[(db, target)] = entry
            while len(self._entries) > self._max_entries:
                self._entries.popitem()
        try:
            key = shared_key(db, target)
            pipe = self._redis().pipeline(transaction=False)
            pipe.hset(key, mapping={"i": now_id, "v": json.dumps(found, sort_keys=True), "b": body})
            pipe.expire(key, self._shared_seconds)
            pipe.execute()
        except (redis.RedisError, OSError):
            pass
