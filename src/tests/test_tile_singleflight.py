# tests/test_tile_singleflight.py

"""
One build of each Galaxy Map tile or stage at a time (PERF.38,
`planetgen/web/lib/singleflight.py`): the claims on a real Redis server, the
tile cache's behaviour around them with the API replaced by a fake.
"""

import threading
import time

import pytest

from planetgen.web.lib import singleflight, tilecache

STAMP = "00000000000000aa"
KEY = "12/1/2/3"


class _Api:
    def __init__(self, slow=0.0):
        self.slow = slow
        self.tile_calls = []
        self.stage_calls = []
        self.lock = threading.Lock()

    def get_galaxy_changes(self, db, since=None):
        return {"stamp": STAMP, "state": "state", "full": since is None, "tiles": [], "stages": []}

    def get_galaxy_tiles(self, db, tile_keys):
        with self.lock:
            self.tile_calls.append(list(tile_keys))
        time.sleep(self.slow)
        tiles = {key: {"stars": [], "clouds": [], "generated": [], "points": []} for key in tile_keys}
        return {"tiles": tiles, "edge_pc": 3.526, "has_shape": True}

    def get_galaxy_stage(self, db, at=None):
        with self.lock:
            self.stage_calls.append(at)
        time.sleep(self.slow)
        return {"at": at, "child_m": 27, "children": [], "sectors": None}


@pytest.fixture
def api(monkeypatch, tmp_path, redis_server):
    fake = _Api(slow=0.4)
    monkeypatch.setattr(tilecache, "get_galaxy_changes", fake.get_galaxy_changes)
    monkeypatch.setattr(tilecache, "get_galaxy_tiles", fake.get_galaxy_tiles)
    monkeypatch.setattr(tilecache, "get_galaxy_stage", fake.get_galaxy_stage)
    monkeypatch.setattr(tilecache, "PRUNE_PROBABILITY", 0.0)
    monkeypatch.setenv("PLANETGEN_TILE_CACHE_DIR", str(tmp_path / "tiles"))
    tilecache.current_stamp("sfdb", tilecache.cache_dir())
    return fake


def _together(*calls):
    results = [None] * len(calls)
    errors = []

    def run(index, call):
        try:
            results[index] = call()
        except Exception as exc:  # noqa: BLE001
            errors.append(exc)

    threads = [threading.Thread(target=run, args=(i, c)) for i, c in enumerate(calls)]
    for thread in threads:
        thread.start()
        time.sleep(0.05)
    for thread in threads:
        thread.join(30)
    assert errors == []
    return results


def test_concurrent_requests_for_one_tile_build_it_once(api):
    results = _together(*[lambda: tilecache.fetch_tiles("sfdb", [KEY])] * 4)
    assert api.tile_calls == [[KEY]]
    assert all(KEY in result["tiles"] for result in results)
    assert sorted(result["cached"] for result in results) == [0, 1, 1, 1]
    assert tilecache.fetch_tiles("sfdb", [KEY])["cached"] == 1


def test_overlapping_requests_build_only_the_tiles_nobody_has_claimed(api):
    other = "12/1/2/4"
    results = _together(lambda: tilecache.fetch_tiles("sfdb", [KEY]),
                        lambda: tilecache.fetch_tiles("sfdb", [KEY, other]))
    assert sorted(sum(api.tile_calls, [])) == [KEY, other]
    assert set(results[1]["tiles"]) == {KEY, other}


def test_concurrent_requests_for_one_stage_build_it_once(api):
    results = _together(*[lambda: tilecache.fetch_stage("sfdb", None)] * 3)
    assert api.stage_calls == [None]
    assert all(result["child_m"] == 27 for result in results)


def test_a_claim_is_released_after_the_build(api):
    tilecache.fetch_tiles("sfdb", [KEY])
    name = tilecache._flight_name(tilecache.cache_dir(), "sfdb", STAMP, "t", KEY)
    owned, held = singleflight.claim([name])
    assert owned == [name] and held == []
    singleflight.release(owned)


def test_a_tile_whose_claimant_never_delivers_is_built_after_the_wait(api, monkeypatch):
    monkeypatch.setattr(singleflight, "WAIT_SECONDS", 0.3)
    name = tilecache._flight_name(tilecache.cache_dir(), "sfdb", STAMP, "t", KEY)
    owned, _held = singleflight.claim([name])  # claimed by somebody who then dies
    try:
        result = tilecache.fetch_tiles("sfdb", [KEY])
    finally:
        singleflight.release(owned)
    assert api.tile_calls[-1] == [KEY] and KEY in result["tiles"]


def test_without_a_lock_service_every_request_builds(api, monkeypatch):
    def down():
        raise OSError("no redis")

    monkeypatch.setattr(singleflight, "_client", down)
    assert singleflight.claim(["a", "b"]) == (["a", "b"], [])
    singleflight.release(["a"])  # no error
    assert KEY in tilecache.fetch_tiles("sfdb", [KEY])["tiles"]
    assert api.tile_calls == [[KEY]]
