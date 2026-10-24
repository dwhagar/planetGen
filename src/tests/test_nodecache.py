# tests/test_nodecache.py

"""
MAP.260: the node cache in `planetgen/web/lib/tilecache.py` (`fetch_nodes`, `fetch_nodes_wire`). Nodes are stored on
disk beside the tiles, once per key and stamp; a change to a sector drops the nodes of its chain and says so in
`history`; the API's calls are replaced with fakes that count calls, so no database is needed.
"""

import os

import pytest

from planetgen.web.lib import tilecache, tilepack

from tests.test_tilecache import FAKE_STAR, FakeApi

NODE = "243.4.1.2"
OTHER = "81.12.3.0"


class NodeApi(FakeApi):
    """`FakeApi` that also answers `get_galaxy_nodes` and reports `changed_nodes` in the change feed."""

    def __init__(self, stamp="00000000000000aa"):
        super().__init__(stamp)
        self.node_calls = []
        self.changed_nodes = []

    def get_galaxy_changes(self, db, since=None):
        answer = super().get_galaxy_changes(db, since)
        answer["nodes"] = [] if answer["full"] else list(self.changed_nodes)
        return answer

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
    fake = NodeApi()
    monkeypatch.setattr(tilecache, "get_galaxy_changes", fake.get_galaxy_changes)
    monkeypatch.setattr(tilecache, "get_galaxy_nodes", fake.get_galaxy_nodes)
    monkeypatch.setattr(tilecache, "get_galaxy_tiles", fake.get_galaxy_tiles)
    monkeypatch.setattr(tilecache, "PRUNE_PROBABILITY", 0.0)
    monkeypatch.setenv("PLANETGEN_TILE_CACHE_DIR", str(tmp_path / "tiles"))
    monkeypatch.delenv("PLANETGEN_CACHE_MAX_MB", raising=False)
    return fake


def test_a_node_keeps_what_the_page_reads_and_is_served_from_disk_next_time(api):
    first = tilecache.fetch_nodes("mydb", [NODE, OTHER])
    assert first["cached"] == 0 and first["stamp"] == api.stamp
    node = first["nodes"][NODE]
    assert node["aggregate"]["lights"]["star_luminosity_sol"] == 1234.5679
    assert node["clouds"][0]["radius_pc"] == 2.3457 and node["points"][0]["x"] == 1.2346
    assert node["placed"] == [{"id": 1, "x": 1.0, "y": 0.0, "z": 0.0}]
    assert "planned" not in node and node["stars"][0]["id"] == 7
    second = tilecache.fetch_nodes("mydb", [OTHER, NODE])
    assert second["cached"] == 2 and second["nodes"] == first["nodes"]
    assert len(api.node_calls) == 1


def test_only_the_missing_nodes_reach_the_api(api):
    tilecache.fetch_nodes("mydb", [NODE])
    tilecache.fetch_nodes("mydb", [NODE, OTHER])
    assert api.node_calls[-1] == [OTHER]


def test_the_wire_body_unpacks_to_the_same_nodes(api):
    body = tilecache.fetch_nodes_wire("mydb", [NODE])
    header, nodes = tilepack.unpack_response(body)
    assert header["stamp"] == api.stamp and list(nodes) == [NODE]
    assert nodes[NODE]["aggregate"]["generated"] == 5


def test_node_keys_are_checked_before_any_api_call(api):
    for keys in (["nonsense"], ["5.1.2.3"], ["243.1.2"], [f"243.4.1.{i}" for i in range(129)], ["1.4.1.2", "1.4.1.x"]):
        with pytest.raises(tilecache.TileRequestError):
            tilecache.fetch_nodes("mydb", keys)
    assert api.node_calls == []


def test_a_sector_key_is_a_node(api):
    assert tilecache.fetch_nodes("mydb", ["1.4.1.2"])["nodes"]["1.4.1.2"]["stars"]


def test_only_the_changed_nodes_are_dropped_and_the_browser_is_told(api, monkeypatch):
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    first = tilecache.fetch_nodes("mydb", [NODE, OTHER])
    api.stamp, api.changed_nodes = "00000000000000bb", [NODE, "27.1.1.1"]
    result = tilecache.fetch_nodes("mydb", [NODE, OTHER], known_stamp="00000000000000aa")
    assert result["generation"] == first["generation"] and result["stamp"] == "00000000000000bb"
    assert result["cached"] == 1 and api.node_calls[-1] == [NODE]
    assert result["history"] == [{"from": "00000000000000aa", "to": "00000000000000bb", "tiles": [],
                                  "nodes": [NODE, "27.1.1.1"]}]


def test_a_full_change_drops_every_node(api, monkeypatch):
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    tilecache.fetch_nodes("mydb", [NODE])
    old = os.path.join(tilecache._db_dir(tilecache.cache_dir(), "mydb"), api.stamp)
    assert os.path.isdir(old)
    api.stamp, api.changed = "00000000000000bb", None
    result = tilecache.fetch_nodes("mydb", [NODE])
    assert result["cached"] == 0 and not os.path.exists(old)


def test_a_dropped_node_is_owed_a_rebuild(api, monkeypatch):
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    tilecache.fetch_nodes("mydb", [NODE, OTHER])
    api.stamp, api.changed_nodes = "00000000000000bb", [NODE]
    tilecache.fetch_nodes("mydb", [OTHER])
    generation, keys = tilecache.take_node_rebuilds("mydb")
    assert keys == [NODE] and generation == "00000000000000aa"
    assert tilecache.take_node_rebuilds("mydb") == (None, [])


def test_nodes_and_tiles_are_stored_apart(api):
    tilecache.fetch_tiles("mydb", ["12/1/2/3"])
    tilecache.fetch_nodes("mydb", [NODE])
    generation = os.path.join(tilecache._db_dir(tilecache.cache_dir(), "mydb"), api.stamp)
    names = sorted(os.listdir(generation))
    assert any(name.startswith("n") for name in names) and any(not name.startswith("n") for name in names)
    assert tilecache.fetch_tiles("mydb", ["12/1/2/3"])["cached"] == 1
