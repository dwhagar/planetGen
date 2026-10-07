# tests/test_caches_threads.py

"""
TEST.54: the web layer's two caches used from real threads at once, as
a threaded WSGI server (or several Apache threads) does.

- The page cache (`planetgen/web/lib/pagecache.py`'s `ResponseCache`): many threads
  filling, reading and clearing it together never raise, never hand back
  one target's body for another (a torn or mixed entry), keep its size
  bookkeeping right, and a `clear()` really empties it -- including of
  answers fetched before the clear and stored after it.
- The tile cache (`planetgen/web/lib/tilecache.py`): two writers of the same tile file
  at once both finish, and the file always ends as one complete, valid
  tile -- never a mix of the two or a truncated file -- whether they
  meet in `_write_json` itself or in two `fetch_tiles` requests for the
  same tile. A reader running alongside only ever sees a whole tile.

No database or API: the page cache gets a stamp function, and the tile
cache's API calls are fakes, as in `test_tilecache.py`.
"""

import json
import os
import sys
import threading

import pytest

_SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _SRC_DIR)

from planetgen.web.lib import pagecache  # noqa: E402
from planetgen.web.lib import tilecache  # noqa: E402

THREADS = 8
ROUNDS = 400
JOIN_TIMEOUT = 30
WRITES = 10


@pytest.fixture
def _switch_often():
    """Hands the GIL over far more often than the default 5 ms, so the
    threads really interleave inside the caches' methods."""
    old = sys.getswitchinterval()
    sys.setswitchinterval(1e-5)
    try:
        yield
    finally:
        sys.setswitchinterval(old)


def _run_threads(targets):
    """Starts one thread per callable, all released together, and
    returns every exception any of them raised."""
    errors = []
    start = threading.Barrier(len(targets))

    def wrap(target):
        def run():
            try:
                start.wait()
                target()
            except BaseException as exc:  # noqa: BLE001 -- reported by the test
                errors.append(exc)
        return run

    threads = [threading.Thread(target=wrap(target), daemon=True) for target in targets]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join(JOIN_TIMEOUT)
    assert not any(thread.is_alive() for thread in threads), "a thread hung"
    return errors


# --- The page cache ---------------------------------------------------------------------

def _body(db, target):
    """Each (db, target) has exactly one possible body, long enough that a
    mix of two would show."""
    return json.dumps({"db": db, "target": target, "pad": target * 40})


def _stored_bytes(cache):
    return sum(len(body) for body, _stored_at in cache._entries.values())


def test_page_cache_fills_reads_and_clears_from_many_threads(_switch_often):
    stamp_calls = []
    stamp_lock = threading.Lock()

    def stamp_for(db):
        # A new stamp every few calls, so entries are dropped by stamp
        # changes too while other threads fill and read.
        with stamp_lock:
            stamp_calls.append(db)
            return f"{db}-{len(stamp_calls) // 25}"

    cache = pagecache.ResponseCache(stamp_for, dict(pagecache.DEFAULTS, max_entries=50, stamp_seconds=0))
    dbs = ("galaxy_a", "galaxy_b", None)
    targets = [f"/sectors/{n}" for n in range(80)]
    served = []
    served_lock = threading.Lock()

    def filler(seed):
        def run():
            hits = []
            for n in range(ROUNDS):
                db = dbs[(seed + n) % len(dbs)]
                target = targets[(seed * 7 + n * 3) % len(targets)]
                generation = cache.generation
                body = cache.get(db, target)
                if body is None:
                    cache.put(db, target, _body(db, target), generation)
                else:
                    hits.append((db, target, body))
            with served_lock:
                served.extend(hits)
        return run

    def clearer():
        for n in range(ROUNDS // 4):
            cache.clear()
            len(cache)
            cache.generation

    errors = _run_threads([filler(seed) for seed in range(THREADS)] + [clearer, clearer])
    assert errors == []
    assert served, "nothing was ever served from the cache"
    for db, target, body in served:
        assert body == _body(db, target), (db, target)
    assert len(cache) <= 50
    assert cache._bytes == _stored_bytes(cache)
    for (db, target), (body, _stored_at) in cache._entries.items():
        assert body == _body(db, target)

    cache.clear()
    assert len(cache) == 0 and cache._bytes == 0
    for db in dbs:
        for target in targets:
            assert cache.get(db, target) is None


def test_page_cache_clear_beats_fills_in_flight_from_many_threads(_switch_often):
    """Answers fetched before a `clear()` (their generation read before
    it) are never stored after it, however the threads interleave."""
    cache = pagecache.ResponseCache(lambda db: "stamp", dict(pagecache.DEFAULTS, stamp_seconds=1e9))
    fetched = threading.Barrier(THREADS + 1)
    cleared = threading.Event()

    def filler(seed):
        def run():
            pending = []
            for n in range(ROUNDS // 4):
                target = f"/systems/{seed}/{n}"
                cache.get("g", target)  # a miss: the answer is "fetched"
                pending.append((target, cache.generation))
            fetched.wait(JOIN_TIMEOUT)  # every "fetch" is in flight...
            cleared.wait(JOIN_TIMEOUT)  # ...when an edit clears the cache
            for target, generation in pending:
                cache.put("g", target, _body("g", target), generation)
        return run

    def editor():
        fetched.wait(JOIN_TIMEOUT)
        cache.clear()
        cleared.set()

    errors = _run_threads([filler(seed) for seed in range(THREADS)] + [editor])
    assert errors == []
    assert len(cache) == 0 and cache._bytes == 0


# --- The tile cache ---------------------------------------------------------------------

def _tile(label, count=1500):
    """A big tile (about 100 KB) whose every entry says who wrote it,
    so half of one and half of another can't pass for either."""
    return {
        "placed": [{"id": n, "name": f"{label}-{n}", "x": n * 0.5, "y": -n * 0.25, "z": 1.0} for n in range(count)],
        "planned": [], "filled": {"g": 1, "cells": []}, "clouds": [],
    }


def _whole_tile(value, tiles):
    return any(value == tile for tile in tiles)


def test_two_writers_of_one_tile_file_leave_one_whole_tile(tmp_path, _switch_often):
    path = os.path.join(str(tmp_path), "db", "gen", "t12_1_2_3.json")
    tiles = [_tile("first"), _tile("second")]
    results = {0: [], 1: []}
    torn = []
    stop = threading.Event()

    def writer(index):
        def run():
            for _ in range(WRITES):
                results[index].append(tilecache._write_json(path, tiles[index]))
        return run

    def reader():
        while not stop.is_set():
            seen = tilecache._read_json(path)
            if seen is not None and not _whole_tile(seen, tiles):
                torn.append(seen)

    reader_thread = threading.Thread(target=reader, daemon=True)
    reader_thread.start()
    try:
        errors = _run_threads([writer(0), writer(1)])
    finally:
        stop.set()
        reader_thread.join(JOIN_TIMEOUT)
    assert errors == []
    assert torn == []
    assert len(results[0]) == len(results[1]) == WRITES
    if os.name != "nt":
        # (Windows may refuse a rename onto a file another thread has
        # open; `_write_json` then reports False and the tile is simply
        # fetched again next time.)
        assert all(results[0]) and all(results[1])
    assert _whole_tile(tilecache._read_json(path), tiles)  # None for a torn file
    leftovers = [name for name in os.listdir(os.path.dirname(path)) if name != os.path.basename(path)]
    assert leftovers == []  # no temp file left behind


class _RacingApi:
    """Two `fetch_tiles` requests for the same tile: each one's
    `get_galaxy_tiles` waits until both are fetching, then answers with
    its own version of the tile."""

    def __init__(self):
        self.stamp = "00000000000000aa"
        self.both_fetching = threading.Barrier(2)
        self.lock = threading.Lock()
        self.answered = []

    def get_galaxy_changes(self, db, since=None):
        return {"stamp": self.stamp, "state": "state", "full": since is None, "tiles": [], "stages": []}

    def get_galaxy_tiles(self, db, tile_keys):
        self.both_fetching.wait(JOIN_TIMEOUT)
        with self.lock:
            label = f"request{len(self.answered)}"
            self.answered.append(label)
        return {"tiles": {key: _tile(label) for key in tile_keys}, "edge_pc": 3.526, "has_shape": True}


def test_two_requests_writing_one_tile_leave_one_whole_tile(monkeypatch, tmp_path):
    api = _RacingApi()
    monkeypatch.setattr(tilecache, "get_galaxy_changes", api.get_galaxy_changes)
    monkeypatch.setattr(tilecache, "get_galaxy_tiles", api.get_galaxy_tiles)
    monkeypatch.setattr(tilecache, "PRUNE_PROBABILITY", 0.0)
    monkeypatch.setenv("PLANETGEN_TILE_CACHE_DIR", str(tmp_path / "tiles"))
    monkeypatch.delenv("PLANETGEN_TILE_CACHE_MAX_MB", raising=False)
    key = "12/1/2/3"
    tilecache.current_stamp("mydb", tilecache.cache_dir())  # stamp.json first, as a warm cache has
    results = []

    def request():
        results.append(tilecache.fetch_tiles("mydb", [key]))

    errors = _run_threads([request, request])
    assert errors == []
    assert sorted(api.answered) == ["request0", "request1"]
    versions = [_tilecache_rounded(_tile(label)) for label in api.answered]
    assert len(results) == 2
    for result in results:
        assert result["cached"] == 0 and _whole_tile(result["tiles"][key], versions)

    # The file both wrote is one of the two, whole; the next request is
    # served from it.
    root = tilecache.cache_dir()
    path = os.path.join(tilecache._db_dir(root, "mydb"), api.stamp, tilecache._tile_filename("t", key))
    with open(path, "r", encoding="utf-8") as f:
        assert _whole_tile(json.load(f), versions)
    again = tilecache.fetch_tiles("mydb", [key])
    assert again["cached"] == 1 and _whole_tile(again["tiles"][key], versions)
    assert len(api.answered) == 2


def _tilecache_rounded(tile):
    return tilecache._round_floats(tile)
