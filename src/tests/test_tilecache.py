# tests/test_tilecache.py

"""
Tests for `planetgen/web/lib/tilecache.py`, the web layer's on-disk cache of 3D
Galaxy Map tiles. The API calls (`get_galaxy_changes`/`get_galaxy_tiles`)
are replaced with fakes that count calls, so no API or database is
needed.
"""

from planetgen.util import settings as settings_model
import json
import os


import pytest  # noqa: E402

from planetgen.web.lib import tilecache  # noqa: E402


FAKE_STAR = {
    "id": 7, "name": "Named", "x": 1.23456789, "y": -2.0, "z": 0.00049, "luminosity_sol": 1234.56789,
    "temperature_k": 5778.4, "radius_sol": 1.234567, "star_type": "G2V Yellow Main Sequence Star",
    "population": "young", "yerkes_class": "V", "ring_index": 3, "layer_index": 1, "ring_slot_index": 4,
    "system_id": None,
}


class FakeApi:
    """`stamp` is the database's current stamp; `changed` the tiles
    `get_galaxy_changes` reports since the last stamp, or `None` for a
    `full` answer."""

    def __init__(self, stamp="00000000000000aa"):
        self.stamp = stamp
        self.changed = []
        self.stamp_calls = 0
        self.since = []
        self.tile_calls = []
        self.changed_stages = []
        self.stage_calls = []

    def get_galaxy_changes(self, db, since=None):
        self.stamp_calls += 1
        self.since.append(since)
        full = since is None or self.changed is None
        return {
            "stamp": self.stamp, "state": "state-" + self.stamp,
            "full": full, "tiles": [] if full else list(self.changed),
            "stages": [] if full else list(self.changed_stages),
        }

    def get_galaxy_stage(self, db, at=None):
        self.stage_calls.append(at)
        return {"at": at, "child_m": 27, "children": [{"ring": 1, "wedge": 2, "slab": 0, "generated": 5}], "sectors": None}

    def get_galaxy_tiles(self, db, tile_keys):
        self.tile_calls.append(list(tile_keys))
        return {
            "tiles": {key: {"placed": [{"id": 1, "x": 1.23456789, "y": 0.0, "z": 0.0}], "planned": [], "filled": {"g": 1},
                            "stars": [dict(FAKE_STAR)], "generated": [], "points": [], "clouds": [
                                {"id": 3, "x": 1.23456789, "y": 0.0, "z": 0.0, "radius_pc": 2.34567891}]}
                      for key in tile_keys},
            "edge_pc": 3.526, "has_shape": True,
        }


@pytest.fixture
def api(monkeypatch, tmp_path):
    fake = FakeApi()
    monkeypatch.setattr(tilecache, "get_galaxy_changes", fake.get_galaxy_changes)
    monkeypatch.setattr(tilecache, "get_galaxy_tiles", fake.get_galaxy_tiles)
    monkeypatch.setattr(tilecache, "get_galaxy_stage", fake.get_galaxy_stage)
    monkeypatch.setattr(tilecache, "PRUNE_PROBABILITY", 0.0)
    monkeypatch.setenv("PLANETGEN_TILE_CACHE_DIR", str(tmp_path / "tiles"))
    monkeypatch.delenv("PLANETGEN_TILE_CACHE_MAX_MB", raising=False)
    return fake


def test_second_request_is_served_from_disk(api):
    first = tilecache.fetch_tiles("mydb", ["12/1/2/3", "12/1/2/4"])
    assert first["cached"] == 0
    assert first["stamp"] == api.stamp
    clouds = first["tiles"]["12/1/2/3"]["clouds"]
    assert clouds[0]["x"] == 1.2346 and clouds[0]["radius_pc"] == 2.3457  # rounded for size
    assert "density" not in first
    assert len(api.tile_calls) == 1

    second = tilecache.fetch_tiles("mydb", ["12/1/2/4", "12/1/2/3"])
    assert second["cached"] == 2
    assert len(api.tile_calls) == 1  # no API call for tiles at all
    assert second["tiles"] == first["tiles"]
    assert second["edge_pc"] == 3.526 and second["has_shape"] is True


def test_only_missing_tiles_reach_the_api(api):
    tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    tilecache.fetch_tiles("mydb", ["12/1/2/3", "12/9/9/9"])
    assert api.tile_calls[-1] == ["12/9/9/9"]


def test_stamp_is_remembered_then_rechecked(api, monkeypatch):
    tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    assert api.stamp_calls == 1
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    assert api.stamp_calls == 2


def test_full_change_drops_every_old_tile(api, monkeypatch):
    tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    root = tilecache.cache_dir()
    db_dir = tilecache._db_dir(root, "mydb")
    assert os.path.isdir(os.path.join(db_dir, api.stamp))

    old_stamp = api.stamp
    api.stamp = "00000000000000bb"
    api.changed = None
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    result = tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    assert result["stamp"] == result["generation"] == "00000000000000bb"
    assert result["cached"] == 0
    assert not os.path.exists(os.path.join(db_dir, old_stamp))


def test_only_changed_tiles_are_dropped(api, monkeypatch):
    first = tilecache.fetch_tiles("mydb", ["12/1/2/3", "12/1/2/4"])
    assert api.since == [None]

    api.stamp = "00000000000000bb"
    api.changed = ["12/1/2/3", "11/0/1/1"]
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    result = tilecache.fetch_tiles("mydb", ["12/1/2/3", "12/1/2/4"])
    assert api.since[-1] == "state-00000000000000aa"
    assert result["stamp"] == "00000000000000bb"
    assert result["generation"] == first["generation"] == "00000000000000aa"
    assert result["cached"] == 1
    assert api.tile_calls[-1] == ["12/1/2/3"]


def test_browser_is_told_which_tiles_changed_since_its_stamp(api, monkeypatch):
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    api.stamp, api.changed = "00000000000000bb", ["12/1/2/3"]
    tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    api.stamp, api.changed = "00000000000000cc", ["12/9/9/9"]

    current = tilecache.fetch_tiles("mydb", ["12/1/2/3"], known_stamp="00000000000000cc")
    assert "history" not in current

    behind = tilecache.fetch_tiles("mydb", ["12/1/2/3"], known_stamp="00000000000000aa")
    assert behind["history"] == [
        {"from": "00000000000000aa", "to": "00000000000000bb", "tiles": ["12/1/2/3"]},
        {"from": "00000000000000bb", "to": "00000000000000cc", "tiles": ["12/9/9/9"]},
    ]
    assert tilecache.fetch_tiles("mydb", ["12/1/2/3"], known_stamp="00000000000000bb")["history"] == behind["history"][1:]
    # The first page render doesn't know the browser's stamp yet.
    assert tilecache.fetch_tiles("mydb", ["12/1/2/3"])["history"] == behind["history"]
    # Too old for the history: the browser drops everything.
    assert tilecache.fetch_tiles("mydb", ["12/1/2/3"], known_stamp="00000000000000ff")["history"] == []


def test_history_is_trimmed_to_its_key_budget(monkeypatch):
    monkeypatch.setattr(tilecache, "HISTORY_MAX_KEYS", 5)
    history = [{"from": str(i), "to": str(i + 1), "tiles": ["k"] * 2} for i in range(4)]
    assert tilecache._trim_history(history) == history[2:]
    assert tilecache._trim_history([{"from": "a", "to": "b", "tiles": ["k"] * 6}]) == []


def test_old_stamp_file_format_starts_afresh(api):
    root = tilecache.cache_dir()
    db_dir = tilecache._db_dir(root, "mydb")
    os.makedirs(os.path.join(db_dir, "00000000000000aa"))
    with open(os.path.join(db_dir, "stamp.json"), "w", encoding="utf-8") as f:
        json.dump({"stamp": "00000000000000aa"}, f)
    result = tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    assert api.since == [None]
    assert result["cached"] == 0 and result["generation"] == "00000000000000aa"


def test_databases_are_cached_separately(api):
    tilecache.fetch_tiles("db_one", ["12/1/2/3"])
    tilecache.fetch_tiles("db_two", ["12/1/2/3"])
    assert len(api.tile_calls) == 2


@pytest.mark.parametrize("keys", [
    ["nope"],
    ["13/0/0/0"],
    ["1/2/0/0"],
    [f"12/{i}/0/0" for i in range(129)],
])
def test_bad_requests_are_rejected_before_any_api_call(api, keys):
    with pytest.raises(tilecache.TileRequestError):
        tilecache.fetch_tiles("mydb", keys)
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
    assert result["tiles"]["12/1/2/3"]["stars"]


def test_corrupt_cache_file_is_refetched(api):
    tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    path = os.path.join(tilecache._db_dir(tilecache.cache_dir(), "mydb"), api.stamp, "t12_1_2_3.json")
    with open(path, "w", encoding="utf-8") as f:
        f.write("{not json")
    result = tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    assert result["cached"] == 0
    with open(path, encoding="utf-8") as f:
        assert json.load(f)["stars"]


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


def test_configured_cache_dir_reports_where_the_cache_belongs(monkeypatch, tmp_path):
    # What examples/apache/create-cache-dir.sh creates -- never created or
    # checked by the helper itself.
    target = tmp_path / "not-yet" / "tiles"
    monkeypatch.setenv("PLANETGEN_TILE_CACHE_DIR", str(target))
    monkeypatch.delenv("PLANETGEN_TILE_CACHE_MAX_MB", raising=False)
    assert tilecache.configured_cache_dir() == str(target)
    assert not target.exists()

    monkeypatch.delenv("PLANETGEN_TILE_CACHE_DIR")
    monkeypatch.setattr(tilecache, "_config", lambda: settings_model.TileCache())
    assert tilecache.configured_cache_dir() == tilecache.DEFAULT_CACHE_DIR

    monkeypatch.setattr(tilecache, "_config", lambda: settings_model.TileCache(max_mb=0))
    assert tilecache.configured_cache_dir() is None


def test_stages_are_cached_until_their_chain_changes(api, monkeypatch):
    first = tilecache.fetch_stage("mydb", "243.7.14.0")
    assert first["children"][0]["generated"] == 5 and first["stamp"] == api.stamp
    tilecache.fetch_stage("mydb", "243.7.14.0")
    tilecache.fetch_stage("mydb", "27.63.115.-1")
    tilecache.fetch_stage("mydb", None)
    tilecache.fetch_stage("mydb", None)
    assert api.stage_calls == ["243.7.14.0", "27.63.115.-1", None]

    api.stamp = "00000000000000bb"
    api.changed = []
    api.changed_stages = ["galaxy", "243.7.14.0"]
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    tilecache.fetch_stage("mydb", "27.63.115.-1")
    tilecache.fetch_stage("mydb", "243.7.14.0")
    tilecache.fetch_stage("mydb", None)
    assert api.stage_calls[3:] == ["243.7.14.0", None]


def test_stage_keys_are_checked(api):
    for bad in ("81.0.0.0", "243.0.99.0", "../x", "1.0.0.0"):
        with pytest.raises(tilecache.TileRequestError):
            tilecache.fetch_stage("mydb", bad)
    assert api.stage_calls == []


def test_a_busy_galaxy_keeps_its_cached_tiles_for_a_while_then_refreshes(api, monkeypatch):
    # PERF.34: a fill changes so much that every check says "full"; the cache is kept (and the stored state
    # with it) until the generation is BUSY_KEEP_SECONDS old, then thrown away once.
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    api.stamp = "00000000000000bb"
    real = api.get_galaxy_changes
    api.get_galaxy_changes = lambda db, since=None: {**real(db, since), "full": True, "busy": True, "tiles": []}
    monkeypatch.setattr(tilecache, "get_galaxy_changes", api.get_galaxy_changes)

    again = tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    assert again["cached"] == 1 and again["stamp"] == "00000000000000aa"
    assert len(api.tile_calls) == 1

    monkeypatch.setattr(tilecache, "BUSY_KEEP_SECONDS", -1)
    refreshed = tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    assert refreshed["cached"] == 0 and refreshed["stamp"] == "00000000000000bb"


def test_a_check_that_fails_under_load_serves_the_cache(api, monkeypatch):
    # PERF.34: a database too busy to answer the freshness check in time must not fail a request the cache can serve.
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    tilecache._failed_checks.clear()
    tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    calls = []

    def failing(db, since=None):
        calls.append(since)
        raise tilecache.apiclient.ApiError("planetGen API error (503): database unavailable", status_code=503)

    monkeypatch.setattr(tilecache, "get_galaxy_changes", failing)
    again = tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    assert again["cached"] == 1
    tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    assert len(calls) == 1  # not asked again straight away
    tilecache._failed_checks.clear()


# --- MAP.157: trimmed tiles, served as stored bytes --------------------------------------

def test_a_tile_keeps_only_what_the_page_reads():
    tile = tilecache.trim_tile({"placed": [{"id": 1}], "planned": [], "filled": {"g": 1}, "clouds": [],
                                "stars": [dict(FAKE_STAR)], "generated": [dict(FAKE_STAR, system_id=12)], "points": []})
    assert set(tile) == {"clouds", "stars", "generated", "points"}
    assert tile["stars"][0] == {
        "id": 7, "x": 1.235, "y": -2.0, "z": 0.0, "luminosity_sol": 1235.0, "temperature_k": 5780,
        "radius_sol": 1.23, "star_type": "G"}
    assert tile["generated"][0]["system_id"] == 12 and "name" not in tile["generated"][0]


def test_a_star_without_a_radius_or_temperature_stays_without(api):
    star = {key: value for key, value in FAKE_STAR.items() if key not in ("radius_sol", "temperature_k")}
    trimmed = tilecache.trim_tile({"stars": [star]})["stars"][0]
    assert "radius_sol" not in trimmed and "temperature_k" not in trimmed


def test_the_wire_body_is_the_stored_bytes_joined(api):
    keys = ["12/1/2/3", "12/1/2/4"]
    cold = tilecache.fetch_tiles_wire("mydb", keys)
    warm = tilecache.fetch_tiles_wire("mydb", keys)
    assert json.loads(warm)["cached"] == 2 and json.loads(cold)["cached"] == 0
    body = json.loads(cold)
    assert list(body["tiles"]) == keys and body["stamp"] == api.stamp and body["edge_pc"] == 3.526
    assert body["tiles"]["12/1/2/3"]["stars"][0]["star_type"] == "G"
    assert len(api.tile_calls) == 1, "the second request was served from disk"


def test_an_empty_request_still_makes_valid_json(api):
    body = tilecache.fetch_tiles_wire("mydb", [])
    assert json.loads(body)["tiles"] == {}


def test_a_torn_tile_file_is_refetched_not_sent(api):
    tilecache.fetch_tiles_wire("mydb", ["12/1/2/3"])
    path = os.path.join(tilecache._db_dir(tilecache.cache_dir(), "mydb"), api.stamp, "t12_1_2_3.json")
    with open(path, "wb") as f:
        f.write(b'{"stars":[')
    body = tilecache.fetch_tiles_wire("mydb", ["12/1/2/3"])
    assert json.loads(body)["tiles"]["12/1/2/3"]["stars"]
    assert len(api.tile_calls) == 2
