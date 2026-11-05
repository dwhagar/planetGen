# tests/test_tilebuild.py

"""
PERF.233: `planetgen/web/tilebuild.py` builds every node of the whole galaxy into a new cache beside the live one and swaps
it in. The API calls are the fakes of `test_tilecache.py` (no API or database); `nodes.all_node_keys` is tried on a real
database below.
"""

import os

import pytest

from planetgen.db import nodes, store
from planetgen.web import tilebuild
from planetgen.web.lib import tilecache
from tests.test_tilecache import FakeApi, api  # noqa: F401 -- the fake API and its fixture
from tests.test_api import seeded_sector  # noqa: F401 -- the fixture

def _keys():
    from planetgen.galaxy import drill, spherical
    keys = []
    for m, shell, wedges in ((243, 3, 2), (81, 12, 2), (27, 40, 1), (9, 200, 1), (3, 700, 1)):
        count = drill.drill_wedge_count(m, shell, 0)
        keys += [drill.format_drill_key(drill.DrillBlock(m, shell, 0, wedge)) for wedge in range(min(wedges, count))]
    count = spherical.wedge_count(20, 0)
    keys += [drill.format_drill_key(drill.DrillBlock(1, 20, 0, wedge)) for wedge in range(min(70, count))]
    return keys


KEYS = _keys()
"""A small galaxy's nodes: five levels of blocks and a sector level of more than one chunk."""


def _live_dir(root):
    return os.path.join(tilecache._db_dir(root, "mydb"), api_generation())


def api_generation():
    return "00000000000000aa"


def _stored(root, keys):
    folder = _live_dir(root)
    return [key for key in keys if os.path.isfile(os.path.join(folder, tilecache._node_filename(key)))]


def _build(**kw):
    kw.setdefault("workers", 1)
    kw.setdefault("keys", KEYS)
    kw.setdefault("free_bytes", 10 ** 12)
    lines = []
    result = tilebuild.build_all("mydb", progress=lines.append, **kw)
    return result, lines


def test_every_node_of_every_level_exists_after_the_build(api, tmp_path):  # noqa: F811
    root = str(tmp_path / "tiles")
    result, lines = _build()
    assert result["swapped"] and result["nodes"] == len(KEYS) and result["built"] == len(KEYS)
    assert _stored(root, KEYS) == KEYS, "every node of every level is in the live cache"
    assert result["levels"] == {243: 2, 81: 2, 27: 1, 9: 1, 3: 1, 1: len(KEYS) - 7}
    for level in result["levels"]:
        assert any(line.startswith(f"level {level}:") for line in lines), (level, lines)
    # What a visit asks for is answered from the disk, not the API.
    api.node_calls.clear()
    visited = tilecache.fetch_nodes("mydb", KEYS[:5])
    assert visited["cached"] == 5 and not api.node_calls
    assert not os.path.isdir(tilebuild.staging_root(root)) or not os.listdir(tilebuild.staging_root(root))


def test_the_old_cache_serves_while_the_new_one_is_built(api, tmp_path, monkeypatch):  # noqa: F811
    root = str(tmp_path / "tiles")
    tilecache.fetch_nodes("mydb", [KEYS[0]])      # the live cache, with one node of its own
    old_file = os.path.join(_live_dir(root), tilecache._node_filename(KEYS[0]))
    assert os.path.isfile(old_file)
    during = []
    real = api.get_galaxy_nodes

    def watching(db, node_keys):
        during.append((os.path.isfile(old_file), os.path.isdir(tilebuild.staging_root(root))))
        return real(db, node_keys)

    monkeypatch.setattr(tilecache, "get_galaxy_nodes", watching)
    _build()
    assert during and all(old_serving and staging for old_serving, staging in during), during
    assert _stored(root, KEYS) == KEYS


def test_a_killed_build_leaves_the_old_cache_and_the_next_run_picks_up(api, tmp_path, monkeypatch):  # noqa: F811
    root = str(tmp_path / "tiles")
    tilecache.fetch_nodes("mydb", [KEYS[0]])
    real = api.get_galaxy_nodes
    calls = []

    def dying(db, node_keys):
        calls.append(list(node_keys))
        if len(calls) > 2:
            raise RuntimeError("killed")
        return real(db, node_keys)

    monkeypatch.setattr(tilecache, "get_galaxy_nodes", dying)
    with pytest.raises(RuntimeError):
        _build()
    assert _stored(root, [KEYS[0]]) == [KEYS[0]], "the old cache still serves"
    assert not _stored(root, KEYS[1:]), "nothing of the half-built cache is live"
    monkeypatch.setattr(tilecache, "get_galaxy_nodes", real)
    result, _lines = _build()
    assert 0 < result["already"] < len(KEYS), "the second run kept what the first built"
    assert _stored(root, KEYS) == KEYS


def test_the_build_ignores_the_cache_budget_and_prunes_nothing(api, tmp_path, monkeypatch):  # noqa: F811
    monkeypatch.setattr(tilecache, "max_cache_bytes", lambda: 1000)     # a budget the first few nodes fill
    monkeypatch.setattr(tilecache, "PRUNE_PROBABILITY", 1.0)
    root = str(tmp_path / "tiles")
    result, _lines = _build()
    assert tilecache.cache_is_full(root), "the budget is far passed"
    assert result["swapped"] and _stored(root, KEYS) == KEYS, "and every node is there all the same"


def test_prune_counts_only_the_nodes_a_visit_asked_for(tmp_path):
    root = str(tmp_path / "cache")
    folder = os.path.join(root, "db", "gen")
    os.makedirs(folder)
    for n in range(5):
        path = os.path.join(folder, f"ahead{n}.bin")
        with open(path, "wb") as handle:
            handle.write(b"x" * 1000)
        os.chmod(path, tilecache.AHEAD_MODE)
    demand = os.path.join(folder, "visited.bin")
    with open(demand, "wb") as handle:
        handle.write(b"x" * 1000)
    tilecache.prune(root, max_bytes=2000)       # 6000 bytes in all, but only 1000 were asked for by a visit
    assert os.path.isfile(demand) and len(os.listdir(folder)) == 6


def test_a_disk_without_twice_the_room_refuses_to_start(api, tmp_path):  # noqa: F811
    with pytest.raises(tilebuild.BuildRefused) as refused:
        tilebuild.build_all("mydb", workers=1, keys=KEYS, free_bytes=1000, progress=lambda line: None)
    assert "GB" in str(refused.value) and "free" in str(refused.value)
    assert not api.node_calls, "nothing was built"
    forced = tilebuild.build_all("mydb", workers=1, keys=KEYS, free_bytes=1000, force=True, progress=lambda line: None)
    assert forced["swapped"]


def test_a_change_in_the_data_meanwhile_is_built_again_before_the_swap(api, tmp_path, monkeypatch):  # noqa: F811
    root = str(tmp_path / "tiles")
    real = api.get_galaxy_nodes
    state = {"changed": False}

    def changing(db, node_keys):
        if not state["changed"] and KEYS[2] in node_keys:
            state["changed"] = True
            api.changed = [KEYS[0]]       # a sector under this block is edited while the build runs
        return real(db, node_keys)

    monkeypatch.setattr(tilecache, "get_galaxy_nodes", changing)
    # The check looks at the database again (not the stamp's 60 s of trust), as a build that takes long does.
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    api.stamp = "00000000000000aa"
    original = api.get_galaxy_changes

    def moving(db, since=None):
        out = original(db, since)
        if state["changed"] and since is not None:
            out["stamp"] = "00000000000000ab"
            out["nodes"] = [KEYS[0]]
        return out

    monkeypatch.setattr(tilecache, "get_galaxy_changes", moving)
    result, lines = _build()
    assert any("taken away" in line for line in lines), lines
    folder = os.path.join(tilecache._db_dir(root, "mydb"))
    assert os.path.isfile(os.path.join(folder, "00000000000000aa", tilecache._node_filename(KEYS[0])))
    tilecache.take_rebuilds("mydb")       # the change owed the dropped node a rebuild: not left for another test


def test_all_node_keys_lists_the_blocks_and_the_sectors_of_a_database(seeded_sector):  # noqa: F811
    from tests.test_api import _place_sector
    config, _sector_id, _systems = seeded_sector
    _place_sector(config, "Placed For Tiles", (40.0, 3.0, 1.0), address=(11, 0, 3))
    conn = store.get_connection(config)
    try:
        keys = nodes.all_node_keys(conn)
    finally:
        conn.close()
    assert len(keys) == len(set(keys)), "each node once"
    assert "1.11.0.3" in keys, "the sector's own node"
    levels = [int(key.split(".", 1)[0]) for key in keys]
    assert levels == sorted(levels, reverse=True), "the biggest blocks first"
