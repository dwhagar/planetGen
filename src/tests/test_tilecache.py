# tests/test_tilecache.py

"""
Tests for `html/lib/tilecache.py`, the web layer's on-disk cache of 3D
Galaxy Map tiles. The API calls (`get_galaxy_stamp`/`get_galaxy_tiles`)
are replaced with fakes that count calls, so no API or database is
needed.
"""

import json
import os
import sys

_SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_SRC_DIR, "html", "lib"))
sys.path.insert(0, _SRC_DIR)

import pytest  # noqa: E402

import tilecache  # noqa: E402


class FakeApi:
    def __init__(self, stamp="00000000000000aa"):
        self.stamp = stamp
        self.stamp_calls = 0
        self.tile_calls = []

    def get_galaxy_stamp(self, db):
        self.stamp_calls += 1
        return self.stamp

    def get_galaxy_tiles(self, db, tile_keys, density_key=None):
        self.tile_calls.append((list(tile_keys), density_key))
        return {
            "tiles": {key: {"placed": [{"id": 1, "x": 1.23456789, "y": 0.0, "z": 0.0}], "planned": []} for key in tile_keys},
            "density": {"key": density_key, "points": [{"x": 1.0, "y": 2.0, "z": 3.0, "relative_density": 0.5}]} if density_key else None,
            "edge_pc": 3.526, "has_shape": True,
        }


@pytest.fixture
def api(monkeypatch, tmp_path):
    fake = FakeApi()
    monkeypatch.setattr(tilecache, "get_galaxy_stamp", fake.get_galaxy_stamp)
    monkeypatch.setattr(tilecache, "get_galaxy_tiles", fake.get_galaxy_tiles)
    monkeypatch.setattr(tilecache, "PRUNE_PROBABILITY", 0.0)
    monkeypatch.setenv("PLANETGEN_TILE_CACHE_DIR", str(tmp_path / "tiles"))
    monkeypatch.delenv("PLANETGEN_TILE_CACHE_MAX_MB", raising=False)
    return fake


def test_second_request_is_served_from_disk(api):
    first = tilecache.fetch_tiles("mydb", ["12/1/2/3", "12/1/2/4"], "5/1/1/1")
    assert first["cached"] == 0
    assert first["stamp"] == api.stamp
    assert first["tiles"]["12/1/2/3"]["placed"][0]["x"] == 1.2346  # rounded for size
    assert first["density"]["points"]
    assert len(api.tile_calls) == 1

    second = tilecache.fetch_tiles("mydb", ["12/1/2/4", "12/1/2/3"], "5/1/1/1")
    assert second["cached"] == 3
    assert len(api.tile_calls) == 1  # no API call for tiles at all
    assert second["tiles"] == first["tiles"]
    assert second["density"] == first["density"]
    assert second["edge_pc"] == 3.526 and second["has_shape"] is True


def test_only_missing_tiles_reach_the_api(api):
    tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    tilecache.fetch_tiles("mydb", ["12/1/2/3", "12/9/9/9"])
    assert api.tile_calls[-1] == (["12/9/9/9"], None)


def test_stamp_is_remembered_then_rechecked(api, monkeypatch):
    tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    assert api.stamp_calls == 1
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    assert api.stamp_calls == 2


def test_new_stamp_drops_old_tiles(api, monkeypatch):
    tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    root = tilecache.cache_dir()
    db_dir = tilecache._db_dir(root, "mydb")
    assert os.path.isdir(os.path.join(db_dir, api.stamp))

    old_stamp = api.stamp
    api.stamp = "00000000000000bb"
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    result = tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    assert result["stamp"] == "00000000000000bb"
    assert result["cached"] == 0
    assert not os.path.exists(os.path.join(db_dir, old_stamp))


def test_databases_are_cached_separately(api):
    tilecache.fetch_tiles("db_one", ["12/1/2/3"])
    tilecache.fetch_tiles("db_two", ["12/1/2/3"])
    assert len(api.tile_calls) == 2


@pytest.mark.parametrize("keys,density", [
    (["nope"], None),
    (["13/0/0/0"], None),
    (["1/2/0/0"], None),
    ([], "garbage"),
    ([f"12/{i}/0/0" for i in range(129)], None),
])
def test_bad_requests_are_rejected_before_any_api_call(api, keys, density):
    with pytest.raises(tilecache.TileRequestError):
        tilecache.fetch_tiles("mydb", keys, density)
    assert api.tile_calls == [] and api.stamp_calls == 0


def test_disabled_cache_still_works(api, monkeypatch):
    monkeypatch.setenv("PLANETGEN_TILE_CACHE_MAX_MB", "0")
    assert tilecache.cache_dir() is None
    tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    second = tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    assert second["cached"] == 0
    assert len(api.tile_calls) == 2


def test_unwritable_cache_dir_falls_back_to_the_api(api, monkeypatch, tmp_path):
    blocker = tmp_path / "a-file"
    blocker.write_text("not a directory")
    monkeypatch.setenv("PLANETGEN_TILE_CACHE_DIR", str(blocker / "tiles"))
    assert tilecache.cache_dir() is None
    result = tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    assert result["tiles"]["12/1/2/3"]["placed"]


def test_corrupt_cache_file_is_refetched(api):
    tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    path = os.path.join(tilecache._db_dir(tilecache.cache_dir(), "mydb"), api.stamp, "t12_1_2_3.json")
    with open(path, "w", encoding="utf-8") as f:
        f.write("{not json")
    result = tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    assert result["cached"] == 0
    with open(path, encoding="utf-8") as f:
        assert json.load(f)["placed"]


def test_prune_deletes_oldest_files_down_to_budget(tmp_path):
    root = tmp_path / "cache"
    stamp_dir = root / "db" / "00000000000000aa"
    stamp_dir.mkdir(parents=True)
    (root / "db" / "stamp.json").write_text('{"stamp": "00000000000000aa"}')
    for i in range(10):
        path = stamp_dir / f"t12_{i}_0_0.json"
        path.write_text("x" * 1000)
        os.utime(path, (1000 + i, 1000 + i))
    tilecache.prune(str(root), max_bytes=5000)
    remaining = sorted(p.name for p in stamp_dir.iterdir())
    assert remaining == [f"t12_{i}_0_0.json" for i in range(6, 10)]
    assert (root / "db" / "stamp.json").exists()
