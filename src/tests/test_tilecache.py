# tests/test_tilecache.py

"""
Tests for `planetgen/web/lib/tilecache.py`, the web layer's on-disk cache of the 3D Galaxy Map's sector-block nodes
(MAP.260). The API calls (`get_galaxy_changes`/`get_galaxy_nodes`) are replaced with fakes that count calls, so no API
or database is needed.
"""

from planetgen.util import settings as settings_model
import json
import os


import pytest  # noqa: E402

from planetgen.web.lib import tilecache, tilepack  # noqa: E402


FAKE_STAR = {
    "id": 7, "name": "Named", "x": 1.23456789, "y": -2.0, "z": 0.00049, "luminosity_sol": 1234.56789,
    "temperature_k": 5778.4, "radius_sol": 1.234567, "star_type": "G2V Yellow Main Sequence Star",
    "population": "young", "yerkes_class": "V", "shell_index": 3, "band_index": 1, "wedge_index": 4,
    "system_id": None,
}


class FakeApi:
    """`stamp` is the database's current stamp; `changed` the nodes
    `get_galaxy_changes` reports since the last stamp, or `None` for a
    `full` answer."""

    def __init__(self, stamp="00000000000000aa"):
        self.stamp = stamp
        self.changed = []
        self.stamp_calls = 0
        self.since = []
        self.node_calls = []
        self.changed_stages = []
        self.stage_calls = []

    def get_galaxy_changes(self, db, since=None):
        self.stamp_calls += 1
        self.since.append(since)
        full = since is None or self.changed is None
        return {
            "stamp": self.stamp, "state": "state-" + self.stamp,
            "full": full, "nodes": [] if full else list(self.changed),
            "stages": [] if full else list(self.changed_stages),
        }

    def get_galaxy_stage(self, db, at=None):
        self.stage_calls.append(at)
        return {"at": at, "child_m": 27, "children": [{"shell": 1, "band": 0, "wedge": 2, "generated": 5}], "sectors": None}

    def get_galaxy_nodes(self, db, node_keys):
        self.node_calls.append(list(node_keys))
        return {
            "nodes": {key: {"aggregate": {"generated": 5, "stats": {"systems": 3.0},
                                          "lights": {"stars": 4, "star_luminosity_sol": 1234.56789}},
                            "stars": [dict(FAKE_STAR)], "generated": [], "points": [{"x": 1.23456789, "y": 0, "z": 0}],
                            "clouds": [{"id": 3, "x": 1.23456789, "y": 0.0, "z": 0.0, "radius_pc": 2.34567891}],
                            "placed": [{"id": 1, "x": 1.0, "y": 0.0, "z": 0.0}], "planned": [{"a": 1}]}
                      for key in node_keys},
            "edge_pc": 3.526, "has_shape": True,
        }


@pytest.fixture
def api(monkeypatch, tmp_path):
    fake = FakeApi()
    monkeypatch.setattr(tilecache, "get_galaxy_changes", fake.get_galaxy_changes)
    monkeypatch.setattr(tilecache, "get_galaxy_nodes", fake.get_galaxy_nodes)
    monkeypatch.setattr(tilecache, "get_galaxy_stage", fake.get_galaxy_stage)
    monkeypatch.setattr(tilecache, "PRUNE_PROBABILITY", 0.0)
    monkeypatch.setenv("PLANETGEN_TILE_CACHE_DIR", str(tmp_path / "tiles"))
    monkeypatch.delenv("PLANETGEN_CACHE_MAX_MB", raising=False)
    return fake


def test_second_request_is_served_from_disk(api):
    first = tilecache.fetch_nodes("mydb", ["243.4.1.2", "243.4.1.3"])
    assert first["cached"] == 0
    assert first["stamp"] == api.stamp
    clouds = first["nodes"]["243.4.1.2"]["clouds"]
    assert clouds[0]["x"] == 1.2346 and clouds[0]["radius_pc"] == 2.3457  # rounded for size
    assert "density" not in first
    assert len(api.node_calls) == 1

    second = tilecache.fetch_nodes("mydb", ["243.4.1.3", "243.4.1.2"])
    assert second["cached"] == 2
    assert len(api.node_calls) == 1  # no API call for nodes at all
    assert second["nodes"] == first["nodes"]
    assert second["edge_pc"] == 3.526 and second["has_shape"] is True


def test_only_missing_nodes_reach_the_api(api):
    tilecache.fetch_nodes("mydb", ["243.4.1.2"])
    tilecache.fetch_nodes("mydb", ["243.4.1.2", "81.12.3.0"])
    assert api.node_calls[-1] == ["81.12.3.0"]


def test_stamp_is_remembered_then_rechecked(api, monkeypatch):
    tilecache.fetch_nodes("mydb", ["243.4.1.2"])
    tilecache.fetch_nodes("mydb", ["243.4.1.2"])
    assert api.stamp_calls == 1
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    tilecache.fetch_nodes("mydb", ["243.4.1.2"])
    assert api.stamp_calls == 2


def test_full_change_drops_every_old_node(api, monkeypatch):
    tilecache.fetch_nodes("mydb", ["243.4.1.2"])
    root = tilecache.cache_dir()
    db_dir = tilecache._db_dir(root, "mydb")
    assert os.path.isdir(os.path.join(db_dir, api.stamp))

    old_stamp = api.stamp
    api.stamp = "00000000000000bb"
    api.changed = None
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    result = tilecache.fetch_nodes("mydb", ["243.4.1.2"])
    assert result["stamp"] == result["generation"] == "00000000000000bb"
    assert result["cached"] == 0
    assert not os.path.exists(os.path.join(db_dir, old_stamp))


def test_only_changed_nodes_are_dropped(api, monkeypatch):
    first = tilecache.fetch_nodes("mydb", ["243.4.1.2", "243.4.1.3"])
    assert api.since == [None]

    api.stamp = "00000000000000bb"
    api.changed = ["243.4.1.2", "27.1.1.1"]
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    result = tilecache.fetch_nodes("mydb", ["243.4.1.2", "243.4.1.3"])
    assert api.since[-1] == "state-00000000000000aa"
    assert result["stamp"] == "00000000000000bb"
    assert result["generation"] == first["generation"] == "00000000000000aa"
    assert result["cached"] == 1
    assert api.node_calls[-1] == ["243.4.1.2"]


def test_browser_is_told_which_nodes_changed_since_its_stamp(api, monkeypatch):
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    tilecache.fetch_nodes("mydb", ["243.4.1.2"])
    api.stamp, api.changed = "00000000000000bb", ["243.4.1.2"]
    tilecache.fetch_nodes("mydb", ["243.4.1.2"])
    api.stamp, api.changed = "00000000000000cc", ["81.12.3.0"]

    current = tilecache.fetch_nodes("mydb", ["243.4.1.2"], known_stamp="00000000000000cc")
    assert "history" not in current

    behind = tilecache.fetch_nodes("mydb", ["243.4.1.2"], known_stamp="00000000000000aa")
    assert behind["history"] == [
        {"from": "00000000000000aa", "to": "00000000000000bb", "nodes": ["243.4.1.2"]},
        {"from": "00000000000000bb", "to": "00000000000000cc", "nodes": ["81.12.3.0"]},
    ]
    assert tilecache.fetch_nodes("mydb", ["243.4.1.2"], known_stamp="00000000000000bb")["history"] == behind["history"][1:]
    # The first page render doesn't know the browser's stamp yet.
    assert tilecache.fetch_nodes("mydb", ["243.4.1.2"])["history"] == behind["history"]
    # Too old for the history: the browser drops everything.
    assert tilecache.fetch_nodes("mydb", ["243.4.1.2"], known_stamp="00000000000000ff")["history"] == []


def test_history_is_trimmed_to_its_key_budget(monkeypatch):
    monkeypatch.setattr(tilecache, "HISTORY_MAX_KEYS", 5)
    history = [{"from": str(i), "to": str(i + 1), "nodes": ["k"] * 2} for i in range(4)]
    assert tilecache._trim_history(history) == history[2:]
    assert tilecache._trim_history([{"from": "a", "to": "b", "nodes": ["k"] * 6}]) == []


def test_old_stamp_file_format_starts_afresh(api):
    root = tilecache.cache_dir()
    db_dir = tilecache._db_dir(root, "mydb")
    os.makedirs(os.path.join(db_dir, "00000000000000aa"))
    with open(os.path.join(db_dir, "stamp.json"), "w", encoding="utf-8") as f:
        json.dump({"stamp": "00000000000000aa"}, f)
    result = tilecache.fetch_nodes("mydb", ["243.4.1.2"])
    assert api.since == [None]
    assert result["cached"] == 0 and result["generation"] == "00000000000000aa"


def test_databases_are_cached_separately(api):
    tilecache.fetch_nodes("db_one", ["243.4.1.2"])
    tilecache.fetch_nodes("db_two", ["243.4.1.2"])
    assert len(api.node_calls) == 2


@pytest.mark.parametrize("keys", [
    ["nope"],
    [f"243.4.1.{i}" for i in range(129)],
    ["243.4.1.2", "243.4.1.x"],
    ["5.1.2.3"],
    ["243.1.2"],
])
def test_bad_requests_are_rejected_before_any_api_call(api, keys):
    with pytest.raises(tilecache.TileRequestError):
        tilecache.fetch_nodes("mydb", keys)
    assert api.node_calls == [] and api.stamp_calls == 0


def test_disabled_cache_still_works(api, monkeypatch):
    monkeypatch.setenv("PLANETGEN_CACHE_MAX_MB", "0")
    assert tilecache.cache_dir() is None
    tilecache.fetch_nodes("mydb", ["243.4.1.2"])
    second = tilecache.fetch_nodes("mydb", ["243.4.1.2"])
    assert second["cached"] == 0
    assert len(api.node_calls) == 2


def test_unwritable_cache_dir_falls_back_to_the_api(api, monkeypatch, tmp_path):
    blocker = tmp_path / "a-file"
    blocker.write_text("not a directory")
    monkeypatch.setenv("PLANETGEN_TILE_CACHE_DIR", str(blocker / "tiles"))
    assert tilecache.cache_dir() is None
    result = tilecache.fetch_nodes("mydb", ["243.4.1.2"])
    assert result["nodes"]["243.4.1.2"]["stars"]


def test_corrupt_cache_file_is_refetched(api):
    tilecache.fetch_nodes("mydb", ["243.4.1.2"])
    path = os.path.join(tilecache._db_dir(tilecache.cache_dir(), "mydb"), api.stamp, "n243_4_1_2.bin")
    with open(path, "wb") as f:
        f.write(b"not a node")
    result = tilecache.fetch_nodes("mydb", ["243.4.1.2"])
    assert result["cached"] == 0
    with open(path, "rb") as f:
        assert tilepack.unpack_node(f.read())["stars"]


def _stored_nodes(tmp_path, count=10, ahead=()):
    """A cache of `count` 1,000-byte node files, the i-th last visited at time 1000 + i; those in `ahead` are marked
    as built ahead."""
    root = tmp_path / "cache"
    stamp_dir = root / "db" / "00000000000000aa"
    stamp_dir.mkdir(parents=True)
    (root / "db" / "stamp.json").write_text('{"stamp": "00000000000000aa"}')
    for i in range(count):
        path = stamp_dir / f"n1_{i}_1_2.bin"
        path.write_text("x" * 1000)
        path.chmod(tilecache.AHEAD_MODE if i in ahead else 0o600)
        os.utime(path, (1000 + i, 1000 + i))
    return root, stamp_dir


def test_prune_deletes_the_least_recently_visited_files_down_to_budget(tmp_path):
    root, stamp_dir = _stored_nodes(tmp_path)
    tilecache.prune(str(root), max_bytes=5000)
    remaining = sorted(p.name for p in stamp_dir.iterdir())
    assert remaining == [f"n1_{i}_1_2.bin" for i in range(6, 10)]
    assert (root / "db" / "stamp.json").exists()


def test_prune_never_deletes_a_node_built_ahead(tmp_path):
    """PERF.133: the oldest nodes here were built ahead; only the ones a visit asked for go, and when only
    built-ahead tiles are left over the budget nothing more goes."""
    root, stamp_dir = _stored_nodes(tmp_path, ahead={0, 1, 2})
    tilecache.prune(str(root), max_bytes=5000)
    remaining = sorted(p.name for p in stamp_dir.iterdir())
    assert remaining == [f"n1_{i}_1_2.bin" for i in (0, 1, 2, 9)]
    tilecache.prune(str(root), max_bytes=1000)
    assert sorted(p.name for p in stamp_dir.iterdir()) == [f"n1_{i}_1_2.bin" for i in (0, 1, 2)]
    assert tilecache.cache_is_full(str(root), max_bytes=1000), "over the budget with only built-ahead nodes left"
    assert not tilecache.cache_is_full(str(root), max_bytes=10_000)


def test_a_visit_refreshes_a_node_and_a_warm_up_does_not(api, monkeypatch):
    first = tilecache.fetch_nodes("mydb", ["243.4.1.2", "243.4.1.3"])
    directory = os.path.join(tilecache._db_dir(tilecache.cache_dir(), "mydb"), first["generation"])
    visited, ahead = (os.path.join(directory, name) for name in ("n243_4_1_2.bin", "n243_4_1_3.bin"))
    for path in (visited, ahead):
        os.utime(path, (1000, 1000))
    tilecache.fetch_nodes("mydb", ["243.4.1.2"])
    tilecache.fetch_nodes("mydb", ["243.4.1.3"], ahead=True)
    assert os.stat(visited).st_mtime > 1000, "a visit counts as the last visit"
    assert os.stat(ahead).st_mtime == 1000, "building ahead does not"


def test_what_is_built_ahead_is_marked_and_what_a_visit_asked_for_is_not(api):
    ahead = tilecache.fetch_nodes("mydb", ["243.4.1.2"], ahead=True)
    asked = tilecache.fetch_nodes("mydb", ["243.4.1.3"])
    directory = os.path.join(tilecache._db_dir(tilecache.cache_dir(), "mydb"), ahead["generation"])
    assert os.stat(os.path.join(directory, "n243_4_1_2.bin")).st_mode & 0o777 == tilecache.AHEAD_MODE
    assert os.stat(os.path.join(directory, "n243_4_1_3.bin")).st_mode & 0o777 == 0o600
    assert asked["generation"] == ahead["generation"]


def test_a_full_cache_drops_visited_nodes_and_keeps_those_built_ahead(api, monkeypatch):
    """The decided behaviour at the budget (Boss, 2026-10-18, 5:56 PM EDT): the warm-up stops, the nodes built ahead
    stay, the least recently visited node a visit asked for goes."""
    keys = [f"1.{i + 4}.1.2" for i in range(8)]
    tilecache.fetch_nodes("mydb", keys[:3], ahead=True)
    for key in keys[3:]:
        tilecache.fetch_nodes("mydb", [key])
    root = tilecache.cache_dir()
    sizes = {path: size for path, size, _mtime, _ahead in tilecache.cache_files(root)}
    one = max(sizes.values())
    generation = tilecache.current_stamp("mydb", root)["generation"]
    directory = os.path.join(tilecache._db_dir(root, "mydb"), generation)
    for i in range(3, 8):
        os.utime(os.path.join(directory, f"n1_{i + 4}_1_2.bin"), (1000 + i, 1000 + i))
    tilecache.prune(root, max_bytes=int(one * 6.5))                # room for six whole nodes at most
    left = sorted((name for name in os.listdir(directory) if name.startswith("n1_")), key=lambda name: int(name.split("_")[1]))
    assert left[:3] == ["n1_4_1_2.bin", "n1_5_1_2.bin", "n1_6_1_2.bin"], "built ahead, kept"
    assert "n1_7_1_2.bin" not in left and "n1_11_1_2.bin" in left, "the least recently visited went, the latest stayed"


def test_configured_cache_dir_reports_where_the_cache_belongs(monkeypatch, tmp_path):
    # What examples/apache/create-cache-dir.sh creates -- never created or
    # checked by the helper itself.
    target = tmp_path / "not-yet" / "tiles"
    monkeypatch.setenv("PLANETGEN_TILE_CACHE_DIR", str(target))
    monkeypatch.delenv("PLANETGEN_CACHE_MAX_MB", raising=False)
    assert tilecache.configured_cache_dir() == str(target)
    assert not target.exists()

    monkeypatch.delenv("PLANETGEN_TILE_CACHE_DIR")
    monkeypatch.setattr(tilecache, "_config", lambda: settings_model.TileCache())
    assert tilecache.configured_cache_dir() == tilecache.DEFAULT_CACHE_DIR

    monkeypatch.setenv("PLANETGEN_CACHE_MAX_MB", "0")
    assert tilecache.configured_cache_dir() is None


def test_stages_are_cached_until_their_chain_changes(api, monkeypatch):
    first = tilecache.fetch_stage("mydb", "243.7.3.0")
    assert first["children"][0]["generated"] == 5 and first["stamp"] == api.stamp
    tilecache.fetch_stage("mydb", "243.7.3.0")
    tilecache.fetch_stage("mydb", "27.63.-2.5")
    tilecache.fetch_stage("mydb", None)
    tilecache.fetch_stage("mydb", None)
    assert api.stage_calls == ["243.7.3.0", "27.63.-2.5", None]

    api.stamp = "00000000000000bb"
    api.changed = []
    api.changed_stages = ["galaxy", "243.7.3.0"]
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    tilecache.fetch_stage("mydb", "27.63.-2.5")
    tilecache.fetch_stage("mydb", "243.7.3.0")
    tilecache.fetch_stage("mydb", None)
    assert api.stage_calls[3:] == ["243.7.3.0", None]


def test_stage_keys_are_checked(api):
    for bad in ("729.0.0.0", "243.0.99.0", "../x", "1.0.0.0"):
        with pytest.raises(tilecache.TileRequestError):
            tilecache.fetch_stage("mydb", bad)
    assert api.stage_calls == []


def test_a_busy_galaxy_keeps_its_cached_nodes_for_a_while_then_refreshes(api, monkeypatch):
    # PERF.34: a fill changes so much that every check says "full"; the cache is kept (and the stored state
    # with it) until the generation is BUSY_KEEP_SECONDS old, then thrown away once.
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    tilecache.fetch_nodes("mydb", ["243.4.1.2"])
    api.stamp = "00000000000000bb"
    real = api.get_galaxy_changes
    api.get_galaxy_changes = lambda db, since=None: {**real(db, since), "full": True, "busy": True, "nodes": []}
    monkeypatch.setattr(tilecache, "get_galaxy_changes", api.get_galaxy_changes)

    again = tilecache.fetch_nodes("mydb", ["243.4.1.2"])
    assert again["cached"] == 1 and again["stamp"] == "00000000000000aa"
    assert len(api.node_calls) == 1

    monkeypatch.setattr(tilecache, "BUSY_KEEP_SECONDS", -1)
    refreshed = tilecache.fetch_nodes("mydb", ["243.4.1.2"])
    assert refreshed["cached"] == 0 and refreshed["stamp"] == "00000000000000bb"


def test_a_check_that_fails_under_load_serves_the_cache(api, monkeypatch):
    # PERF.34: a database too busy to answer the freshness check in time must not fail a request the cache can serve.
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    tilecache._failed_checks.clear()
    tilecache.fetch_nodes("mydb", ["243.4.1.2"])
    calls = []

    def failing(db, since=None):
        calls.append(since)
        raise tilecache.apiclient.ApiError("planetGen API error (503): database unavailable", status_code=503)

    monkeypatch.setattr(tilecache, "get_galaxy_changes", failing)
    again = tilecache.fetch_nodes("mydb", ["243.4.1.2"])
    assert again["cached"] == 1
    tilecache.fetch_nodes("mydb", ["243.4.1.2"])
    assert len(calls) == 1  # not asked again straight away
    tilecache._failed_checks.clear()


# --- MAP.157: trimmed nodes, served as stored bytes --------------------------------------

def test_a_node_keeps_what_the_page_reads():
    node = tilecache.trim_node({"placed": [{"id": 1}], "planned": [{"a": 1}], "filled": {"g": 1}, "clouds": [],
                                "stars": [dict(FAKE_STAR)], "generated": [dict(FAKE_STAR, system_id=12)], "points": [],
                                "aggregate": {"lights": {"star_luminosity_sol": 1234.56789}}})
    assert set(node) == {"aggregate", "clouds", "stars", "generated", "points", "placed"}
    assert node["aggregate"]["lights"]["star_luminosity_sol"] == 1234.5679
    assert node["stars"][0] == {
        "id": 7, "x": 1.235, "y": -2.0, "z": 0.0, "luminosity_sol": 1235.0, "temperature_k": 5780,
        "radius_sol": 1.23, "star_type": "G"}
    assert node["generated"][0]["system_id"] == 12 and "name" not in node["generated"][0]


def test_a_star_without_a_radius_or_temperature_stays_without(api):
    star = {key: value for key, value in FAKE_STAR.items() if key not in ("radius_sol", "temperature_k")}
    trimmed = tilecache.trim_node({"stars": [star]})["stars"][0]
    assert "radius_sol" not in trimmed and "temperature_k" not in trimmed


def test_the_wire_body_is_the_stored_bytes_joined(api):
    keys = ["243.4.1.2", "243.4.1.3"]
    cold_body = tilecache.fetch_nodes_wire("mydb", keys)
    warm_body = tilecache.fetch_nodes_wire("mydb", keys)
    warm, _ = tilepack.unpack_response(warm_body)
    cold, nodes = tilepack.unpack_response(cold_body)
    assert warm["cached"] == 2 and cold["cached"] == 0
    assert list(nodes) == keys and cold["stamp"] == api.stamp and cold["edge_pc"] == 3.526
    assert nodes["243.4.1.2"]["stars"][0]["star_type"] == "G"
    assert len(api.node_calls) == 1, "the second request was served from disk"
    path = os.path.join(tilecache._db_dir(tilecache.cache_dir(), "mydb"), api.stamp, "n243_4_1_2.bin")
    with open(path, "rb") as f:
        assert f.read() in warm_body, "a node goes out as the bytes the cache stores"


def test_an_empty_request_still_makes_a_valid_response(api):
    header, nodes = tilepack.unpack_response(tilecache.fetch_nodes_wire("mydb", []))
    assert nodes == {} and header["stamp"] == api.stamp


def test_a_torn_node_file_is_refetched_not_sent(api):
    tilecache.fetch_nodes_wire("mydb", ["243.4.1.2"])
    path = os.path.join(tilecache._db_dir(tilecache.cache_dir(), "mydb"), api.stamp, "n243_4_1_2.bin")
    with open(path, "rb") as f:
        whole = f.read()
    with open(path, "wb") as f:
        f.write(whole[:-5])
    _, nodes = tilepack.unpack_response(tilecache.fetch_nodes_wire("mydb", ["243.4.1.2"]))
    assert nodes["243.4.1.2"]["stars"]
    assert len(api.node_calls) == 2


def test_a_sector_key_is_a_node(api):
    assert tilecache.fetch_nodes("mydb", ["1.4.1.2"])["nodes"]["1.4.1.2"]["stars"]


def test_the_wire_body_carries_the_aggregate(api):
    header, nodes = tilepack.unpack_response(tilecache.fetch_nodes_wire("mydb", ["243.4.1.2"]))
    assert header["stamp"] == api.stamp and nodes["243.4.1.2"]["aggregate"]["generated"] == 5


def test_a_dropped_node_is_owed_a_rebuild(api, monkeypatch):
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    tilecache.fetch_nodes("mydb", ["243.4.1.2", "81.12.3.0"])
    api.stamp, api.changed = "00000000000000bb", ["243.4.1.2"]
    tilecache.fetch_nodes("mydb", ["81.12.3.0"])
    generation, keys = tilecache.take_rebuilds("mydb")
    assert keys == ["243.4.1.2"] and generation == "00000000000000aa"
    assert tilecache.take_rebuilds("mydb") == (None, [])
