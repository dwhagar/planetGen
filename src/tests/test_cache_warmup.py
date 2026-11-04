# tests/test_cache_warmup.py

"""
PERF.104, PERF.127: the warm-up (`planetgen/web/warmup.py`, `planetgen warm`) fills the shared
caches so the first visit answers from them, one command per family, loudly.

The families visit the pages of a real app on a throwaway database; a second, separate app
(a web worker that never saw the warm-up) is then asked for the same pages with a spy on the API
calls, to show nothing is rebuilt. The wrapper's list of the commands runs on the real Redis queue (the
tests' workers run the units here), with the families themselves faked where only the list is under test.
"""

import time
import types
import uuid

import pytest

from planetgen.admin import maintenance
from planetgen.api.config import Config
from planetgen.cli import generate as generate_cli
from planetgen.cli import warm as warm_cli
from planetgen.cli import warm_progress
from planetgen.db import block_stats, dataversion, store
from planetgen.galaxy import spherical
from planetgen.galaxy.sector import SpaceSector
from planetgen.physics.units import ly_to_pc
from planetgen.tuning import DEFAULT_SECTOR_EDGE_LY
from planetgen.util.settings import get_settings
from planetgen.queue import commands, redisclient, redisqueue, runrecord, units, wrapper
from planetgen.web import warmup
from planetgen.web.app import create_app
from planetgen.web.lib import apiclient, tilecache
from tests.publicids import pid
from tests.test_api import _place_sector, seeded_sector  # noqa: F401 -- the fixture


@pytest.fixture
def galaxy(seeded_sector, redis_server, tmp_path, monkeypatch):  # noqa: F811
    """The seeded sector (two systems) and an empty one, a nebula in the first, and a node cache of its own."""
    config, sector_id, system_ids = seeded_sector
    monkeypatch.setenv("PLANETGEN_TILE_CACHE_DIR", str(tmp_path / "tiles"))
    monkeypatch.setattr(tilecache, "PRUNE_PROBABILITY", 0.0)
    other_id = store.save_sector(SpaceSector(name="Other Sector", edge_ly=10.0), config=config)
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.phenomena.nebula import Nebula

    conn = store.get_connection(config)
    try:
        nebula_id = store.insert_nebula(conn, Nebula(SystemConfig(), name="Warm Cloud"), placement={
            "center_x_pc": 5.0, "center_y_pc": 1.0, "center_z_pc": 2.0, "galactic_radius_pc": 5.5})
        conn.commit()
        pid("nebula", nebula_id, conn)   # issues the printed ID the pages and the API use
    finally:
        conn.close()

    class Galaxy:
        pass

    found = Galaxy()
    found.config, found.sector, found.other, found.systems = config, sector_id, other_id, system_ids
    found.sector_pid, found.other_pid = pid("sector", sector_id, config), pid("sector", other_id, config)
    found.system_pids = [pid("system", system_id, config) for system_id in system_ids]
    yield found
    for key in redisclient.client(decode_responses=True).scan_iter(f"planetgen:*:{config.database}*"):
        redisclient.client().delete(key)


class _Visitor(Config):
    SESSION_COOKIE_SECURE = False
    SECRET_KEY = "test-secret"


def _web_worker(config):
    """A web process that has not seen the warm-up: its own empty in-process copies, the shared Redis."""

    class Site(_Visitor):
        MYSQL_CONFIG = WRITE_MYSQL_CONFIG = CONTROL_MYSQL_CONFIG = config
        WEB_DATABASE = ""

    app = create_app(Site)
    app.testing = True
    return app.test_client()


@pytest.fixture
def spy(monkeypatch):
    seen = []
    real = apiclient._send

    def counting(method, target, *args, **kwargs):
        if method == "GET":
            seen.append(target.partition("?")[0])
        return real(method, target, *args, **kwargs)
    monkeypatch.setattr(apiclient, "_send", counting)
    return seen


def _warm(galaxy, family, **options):
    return warmup.warm_family(galaxy.config.database, family, app=warmup.build_app(galaxy.config), **options)


# --- The families ------------------------------------------------------------------------------


def test_the_first_visit_to_a_sector_and_a_system_is_answered_from_the_cache(galaxy, spy):
    sectors = _warm(galaxy, "sectors")
    systems = _warm(galaxy, "systems")
    assert sectors["warmed"] == 2 and sectors["left_out"] == 0 and sectors["built"] > 0
    assert systems["warmed"] == 2 and "systems: 2 systems warmed" in systems["line"]
    web = _web_worker(galaxy.config)
    del spy[:]
    for path in (f"/sector/{galaxy.sector_pid}", f"/sector/{galaxy.other_pid}", f"/system/{galaxy.system_pids[0]}",
                 f"/system/{galaxy.system_pids[1]}", f"/system/{galaxy.system_pids[0]}/scene",
                 f"/sector/{galaxy.sector_pid}/scene"):
        assert web.get(path).status_code == 200, path
    rebuilt = [call for call in spy if call.startswith(("/sectors/", "/systems/"))]
    assert rebuilt == [], f"the first visit rebuilt {rebuilt}"


def test_the_first_visit_to_a_list_and_the_search_is_answered_from_the_cache(galaxy, spy):
    lists = _warm(galaxy, "lists")
    search = _warm(galaxy, "search")
    from planetgen.web.object_lists import KIND_LABELS
    assert lists["warmed"] == 7 + len(KIND_LABELS) and search["warmed"] == 1
    web = _web_worker(galaxy.config)
    del spy[:]
    for path in ("/sectors", "/systems", "/phenomena", "/", "/search", "/lists", "/lists/planet", "/lists/black_hole"):
        assert web.get(path).status_code == 200, path
    rebuilt = [call for call in spy if call != "/galaxy/changes"]   # the stamp check is never cached
    assert rebuilt == [], f"the first visit rebuilt {rebuilt}"


def test_the_first_visit_to_a_phenomenon_is_answered_from_the_cache(galaxy, spy):
    result = _warm(galaxy, "phenomena")
    assert result["warmed"] == 1
    web = _web_worker(galaxy.config)
    envelope = web.get(f"/api/phenomena?db={galaxy.config.database}").get_json()
    [nebula] = [row for row in envelope["items"] if row["name"] == "Warm Cloud"]
    del spy[:]
    assert web.get(f"/phenomenon/{nebula['type']}/{nebula['id']}").status_code == 200
    assert [call for call in spy if call.startswith("/phenomena/")] == []


def test_the_galaxy_map_family_builds_the_opening_view_and_the_page(galaxy, spy):
    result = _warm(galaxy, "galaxy-map")
    assert result["noun"] == "nodes" and result["warmed"] > 0
    web = _web_worker(galaxy.config)
    del spy[:]
    assert web.get("/galaxy").status_code == 200
    assert [call for call in spy if call.startswith(("/galaxy/sectors", "/galaxy/shape", "/galaxy/phenomena"))] == []


def test_the_galaxy_map_family_warms_while_the_site_is_held_in_maintenance(galaxy, tmp_path, monkeypatch):
    """The update warms the caches with the site closed; the page's maintenance redirect is not the warm-up's failure."""
    monkeypatch.setenv(maintenance.FILE_ENV, str(tmp_path / "maintenance.json"))
    maintenance.turn_on(by="update")
    assert _web_worker(galaxy.config).get("/galaxy").status_code == 302, "a visitor is still sent to the maintenance page"
    result = _warm(galaxy, "galaxy-map")
    assert result["warmed"] > 0


def test_the_galaxy_map_family_warms_an_empty_galaxy_in_maintenance(mysql_config, redis_server, tmp_path, monkeypatch):
    """After a reset the galaxy holds nothing and the update still has to get through stage 1."""
    monkeypatch.setenv(maintenance.FILE_ENV, str(tmp_path / "maintenance.json"))
    monkeypatch.setenv("PLANETGEN_TILE_CACHE_DIR", str(tmp_path / "tiles"))
    maintenance.turn_on(by="update")
    result = warmup.warm_family(mysql_config.database, "galaxy-map", app=warmup.build_app(mysql_config))
    assert result["family"] == "galaxy-map"


def test_a_second_warm_up_finds_it_all_built(galaxy):
    first = _warm(galaxy, "sectors")
    second = _warm(galaxy, "sectors")
    assert first["built"] > 0 and second["built"] == 0 and second["cached"] >= first["built"]


def test_a_cap_leaves_the_rest_and_says_so(galaxy):
    result = _warm(galaxy, "sectors", cap=1)
    assert result["warmed"] == 1 and result["left_out"] == 1
    assert "1 left out by the cap" in result["line"]


def test_a_page_that_does_not_answer_fails_the_family(galaxy, monkeypatch):
    real = apiclient.get_system

    def broken(db, system_id):
        if system_id == galaxy.system_pids[0]:
            raise apiclient.ApiError("boom", status_code=500)
        return real(db, system_id)
    monkeypatch.setattr(apiclient, "get_system", broken)
    with pytest.raises(warmup.WarmupFailed) as failure:
        _warm(galaxy, "systems")
    assert "did not answer" in str(failure.value) and galaxy.system_pids[0] in str(failure.value)


def test_an_unknown_family_is_refused(galaxy):
    with pytest.raises(warmup.WarmupFailed):
        _warm(galaxy, "nonsense")


def test_changing_one_sector_after_the_warm_up_rebuilds_only_it(galaxy, spy):
    _warm(galaxy, "sectors")
    dataversion.record(galaxy.config.database, "sector-edited", sector_id=galaxy.other_pid)
    web = _web_worker(galaxy.config)
    del spy[:]
    for pid_ in (galaxy.sector_pid, galaxy.other_pid):
        assert web.get(f"/sector/{pid_}").status_code == 200
    assert f"/sectors/{galaxy.other_pid}" in spy and f"/sectors/{galaxy.sector_pid}" not in spy


# --- The commands ----------------------------------------------------------------------------------


def test_a_family_command_runs_alone_and_warms_only_its_family(galaxy, monkeypatch):
    real = warmup.build_app
    monkeypatch.setattr(warmup, "build_app", lambda config=None: real(galaxy.config))
    result = commands.run_alone("warm-search", {"database": galaxy.config.database, "cap": None})
    assert result["family"] == "search" and result["line"].startswith("search: 1 pages warmed")
    # Every registered command is cut into units now, so the refusal is shown with an inline one made for the test.
    monkeypatch.setitem(commands.REGISTRY, "inline-demo", commands.CommandSpec("inline-demo", "Not cut into units"))
    with pytest.raises(ValueError, match="cannot be cut"):
        commands.expand_function("inline-demo")


def test_every_cache_family_is_its_own_registered_command_with_one_unit():
    assert list(commands.KIND_OF_WARM_FAMILY) == list(warmup.FAMILIES)
    for family, kind in commands.KIND_OF_WARM_FAMILY.items():
        assert not commands.spec(kind).inline
        if family in ("cluster-field", "galaxy-nodes", "galaxy-deep-nodes", "sectors", "systems", "phenomena"):
            continue   # many units, one a part of the cluster field, a chunk of nodes or a slice of a listing: their own tests
        expansion = commands.expand_function(kind)({"database": "x", "cap": 4})
        assert expansion.units == [({"database": "x", "cap": 4}, 1.0)], "one unit that warms that family and no other"
    assert commands.spec(commands.END_MAP_REBUILD).summary


@pytest.fixture
def object_site(galaxy, monkeypatch):
    real = warmup.build_app
    monkeypatch.setattr(warmup, "build_app", lambda config=None: real(galaxy.config))
    return galaxy


@pytest.mark.parametrize("family, kind", [("sectors", "warm-sectors"), ("systems", "warm-systems"), ("phenomena", "warm-phenomena")])
def test_an_object_family_is_cut_into_slices_that_warm_every_object_once_and_none_is_longer_than_the_ceiling(
        object_site, monkeypatch, family, kind):
    """PERF.180: the sector, system and phenomenon pages were one unit each, minutes on a large galaxy; each is now slices of its
    listing, `OBJECT_CHUNK` objects (a whole multiple past `units.MAX_UNITS` slices), and the longest the command declares fits the ceiling."""
    database = object_site.config.database
    monkeypatch.setattr(warmup, "OBJECT_CHUNK", 1)
    expansion = commands.expand_function(kind)({"database": database, "cap": None})
    slices = [args for args, _cost in expansion.units]
    assert [args["offset"] for args in slices] == list(range(len(slices))) and all(args["count"] in (0, 1) for args in slices)
    results = [commands.run_alone(kind, args) for args in slices]
    assert sum(result["warmed"] for result in results) == sum(args["count"] for args in slices)
    assert sum(result["left_out"] for result in results) == 0
    spec = commands.spec(kind)
    assert not spec.whole and 0 < spec.unit_seconds <= 300
    assert spec.unit_seconds >= warmup.OBJECT_SECONDS * 250      # the declared longest is for the real slice size


def test_a_cap_below_the_listing_is_left_out_by_the_first_slice_alone(object_site, monkeypatch):
    database = object_site.config.database
    monkeypatch.setattr(warmup, "OBJECT_CHUNK", 1)
    expansion = commands.expand_function("warm-sectors")({"database": database, "cap": 1})
    [(first, _cost)] = expansion.units
    assert (first["offset"], first["count"]) == (0, 1) and first["left_out"] >= 1
    result = commands.run_alone("warm-sectors", first)
    assert result["warmed"] == 1 and result["left_out"] == first["left_out"] and "left out by the cap" in result["line"]


def test_a_slice_warms_the_objects_at_its_place_in_the_listing(object_site):
    database = object_site.config.database
    app = warmup.build_app(object_site.config)
    everything = warmup.warm_family(database, "systems", app=app)
    second = warmup.warm_family(database, "systems", app=app, window={"offset": 1, "count": 1, "left_out": 0})
    assert everything["warmed"] >= 2 and second["warmed"] == 1 and second["left_out"] == 0


# --- The Galaxy Map's nodes (PERF.130, MAP.260) --------------------------------------------------

DEEP_SLICE = 40
"""int: How many of the zoom run's nodes a test builds (the real run is some thousands)."""


@pytest.fixture
def node_levels(galaxy, monkeypatch):
    """`warm-galaxy-nodes` and `warm-galaxy-deep-nodes` on the seeded galaxy: a zoom run to blocks of 81 shells, cut to
    its first `DEEP_SLICE` nodes."""
    real = warmup.build_app
    monkeypatch.setattr(warmup, "build_app", lambda config=None: real(galaxy.config))
    levels = types.SimpleNamespace(warm_zoom_floor=81)
    monkeypatch.setattr(warmup, "get_settings", lambda: types.SimpleNamespace(tile_cache=levels))
    whole = warmup.zoom_run_nodes
    monkeypatch.setattr(warmup, "zoom_run_nodes", lambda *args: whole(*args)[:DEEP_SLICE])
    return levels


def _run_deep_units(galaxy):
    expansion = commands.expand_function("warm-galaxy-deep-nodes")({"database": galaxy.config.database, "cap": None})
    return expansion, [commands.run_alone("warm-galaxy-deep-nodes", args) for args, _cost in expansion.units]


def test_the_zoom_run_is_the_nodes_between_the_floor_and_the_built_levels_coarsest_first():
    keys = warmup.deep_nodes_to_warm(None, 81)
    levels = [int(key.split(".")[0]) for key in keys]
    assert keys and len(set(keys)) == len(keys)
    assert levels == sorted(levels, reverse=True) and max(levels) < warmup.NODE_WARM_FLOOR and min(levels) >= 81
    finer = warmup.deep_nodes_to_warm(None, 27)
    assert {int(key.split(".")[0]) for key in finer} >= {81, 27} and set(keys) <= set(finer)
    assert warmup.deep_nodes_to_warm(None, warmup.NODE_WARM_FLOOR) == [], "no zoom run below the built levels"


def test_a_floor_that_would_build_too_many_nodes_is_refused_with_the_count(monkeypatch):
    monkeypatch.setattr(warmup, "MAX_WARM_NODES", 10)
    with pytest.raises(warmup.WarmupFailed, match=r"blocks of 27 shells would build \d+ nodes.*raise tile_cache.warm_zoom_floor"):
        warmup.deep_nodes_to_warm(None, 27)


def _place_addressed(config, name, center_pc):
    """A sector with the grid address its position falls in, as the generation commands save it (the nodes over it are
    named by the address)."""
    edge_pc = ly_to_pc(DEFAULT_SECTOR_EDGE_LY)
    sector_id = _place_sector(config, name, center_pc, address=spherical.sector_address_at(center_pc, edge_pc))
    conn = store.get_connection(config)
    try:
        block_stats.rebuild(conn)           # the blocks over the sector, which `warm-galaxy-nodes` lists
        conn.commit()
    finally:
        conn.close()
    return sector_id


def _rename_sector(galaxy, sector_id):
    conn = store.get_connection(galaxy.config)
    try:
        conn.execute("UPDATE sectors SET name = ? WHERE id = ?", ("Renamed By A Test", sector_id))
        conn.commit()
    finally:
        conn.close()


def test_the_warm_up_stops_when_the_cache_is_full_and_what_it_built_stays(galaxy, node_levels, monkeypatch):
    """PERF.139 / PERF.133 (Boss, 2026-10-18, 5:56 PM EDT): at the `cache.max_mb` budget the warm-up builds nothing
    more, says how many nodes it left out, and `prune` keeps the nodes it built."""
    monkeypatch.setattr(warmup, "NODE_CHUNK", 4)
    monkeypatch.setenv("PLANETGEN_CACHE_MAX_MB", "0.0004")               # about 400 bytes: the first chunk fills it
    expansion, results = _run_deep_units(galaxy)
    assert len(expansion.units) > 3
    built = [result for result in results if result["built"]]
    left_out = [result for result in results if result["left_out"]]
    assert built and left_out and results.index(left_out[0]) > results.index(built[0]), "full after some, not before"
    assert all(not result["warmed"] for result in left_out) and "cache is full" in left_out[0]["line"]
    assert sum(result["built"] for result in results) < sum(len(args["keys"]) for args, _cost in expansion.units)
    root = tilecache.cache_dir()
    before = sorted(path for path, _size, _mtime, _ahead in tilecache.cache_files(root))
    assert tilecache.cache_is_full(root)
    tilecache.prune(root)
    assert sorted(path for path, _size, _mtime, _ahead in tilecache.cache_files(root)) == before, "built ahead, kept"


def test_a_stage_with_more_chunks_than_the_unit_cap_is_cut_into_bigger_units_and_still_warms_every_node(
        galaxy, node_levels, monkeypatch):
    """Boss, 2026-10-24: "unit 0 of stage 3 failed: a stage queues at most 4000 units, not 4051"."""
    monkeypatch.setattr(warmup, "NODE_CHUNK", 4)
    plain = commands.expand_function("warm-galaxy-deep-nodes")({"database": galaxy.config.database, "cap": None})
    keys = [key for args, _cost in plain.units for key in args["keys"]]
    assert len(plain.units) > 6
    monkeypatch.setattr(units, "MAX_UNITS", 6)
    capped = commands.expand_function("warm-galaxy-deep-nodes")({"database": galaxy.config.database, "cap": None})
    assert len(capped.units) <= 6 and [key for args, _cost in capped.units for key in args["keys"]] == keys
    assert all(len(args["keys"]) % 4 == 0 for args, _cost in capped.units[:-1])
    asked = []
    real = warmup.fetch_nodes
    monkeypatch.setattr(warmup, "fetch_nodes",
                        lambda db, wanted, **kwargs: asked.append(len(wanted)) or real(db, wanted, **kwargs))
    results = [commands.run_alone("warm-galaxy-deep-nodes", args) for args, _cost in capped.units]
    assert sum(result["warmed"] for result in results) == len(keys) and sum(result["left_out"] for result in results) == 0
    assert max(asked) <= 4, "a request never holds more than NODE_CHUNK nodes, however big the unit"
    assert sum(asked) == len(keys)


def test_a_big_unit_stops_part_way_when_the_cache_fills_and_counts_the_rest_as_left_out(galaxy, node_levels, monkeypatch):
    monkeypatch.setattr(warmup, "NODE_CHUNK", 2)
    monkeypatch.setenv("PLANETGEN_CACHE_MAX_MB", "0.0004")
    keys = warmup.deep_nodes_to_warm(None, 81)[:8]
    result = commands.run_alone("warm-galaxy-deep-nodes", {"database": galaxy.config.database, "keys": keys})
    assert 0 < result["warmed"] < len(keys) and result["warmed"] + result["left_out"] == len(keys)
    assert result["warmed"] % 2 == 0 and "left out as the cache filled" in result["line"]


def test_a_cache_with_room_warms_everything_and_says_nothing_was_left_out(galaxy, node_levels):
    expansion, results = _run_deep_units(galaxy)
    assert sum(result["left_out"] for result in results) == 0
    assert sum(result["warmed"] for result in results) == DEEP_SLICE == expansion.plan["nodes"]


def test_a_second_warm_up_with_nothing_changed_builds_no_node(galaxy, node_levels, monkeypatch):
    """PERF.133 / Boss 2026-10-18 6:05 PM EDT: a node is not compiled again when nothing it holds has changed."""
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)   # every call asks the API what changed
    _expansion, first = _run_deep_units(galaxy)
    assert sum(result["built"] for result in first) > 0
    _expansion, again = _run_deep_units(galaxy)
    assert sum(result["built"] for result in again) == 0
    assert sum(result["cached"] for result in again) == sum(result["warmed"] for result in again) > 0
    assert "0 built" in again[0]["line"] and "already stored" in again[0]["line"]


def test_a_change_inside_one_sector_rebuilds_only_the_nodes_holding_it(galaxy, node_levels, monkeypatch):
    placed = _place_addressed(galaxy.config, "Placed", (5.0, 1.0, 2.0))
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    _expansion, first = _run_deep_units(galaxy)
    time.sleep(0.01)
    _rename_sector(galaxy, placed)
    _expansion, again = _run_deep_units(galaxy)
    rebuilt = sum(result["built"] for result in again)
    assert 0 < rebuilt < sum(result["warmed"] for result in first), "the nodes over the sector, not the rest"
    _expansion, settled = _run_deep_units(galaxy)
    assert sum(result["built"] for result in settled) == 0


def test_a_release_alone_does_not_rebuild_a_node(galaxy, node_levels, monkeypatch):
    """PERF.133 / PERF.41: the stored nodes survive a release that changes neither the node format nor the epoch."""
    from planetgen.db import query
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    _expansion, first = _run_deep_units(galaxy)
    assert sum(result["built"] for result in first) > 0
    monkeypatch.setattr(query, "__version__", "99.0.0")
    _expansion, after_release = _run_deep_units(galaxy)
    assert sum(result["built"] for result in after_release) == 0
    # A new node format, though, drops them all once.
    monkeypatch.setattr(query, "TILE_FORMAT", query.TILE_FORMAT + 1)
    _expansion, after_format = _run_deep_units(galaxy)
    assert sum(result["built"] for result in after_format) == sum(result["warmed"] for result in first)


def _all_settings(monkeypatch, levels):
    """The warm-up's settings as they are, except the zoom floor of `node_levels` (which fakes the whole lot)."""
    real = get_settings()

    class Settings:
        tile_cache = levels

        def __getattr__(self, name):
            return getattr(real, name)

    monkeypatch.setattr(warmup, "get_settings", lambda: Settings())


def _phase_lines(caplog):
    return [record.getMessage() for record in caplog.records if record.getMessage().startswith("Phase ")]


def test_a_small_galaxy_runs_through_every_phase_and_a_second_run_changes_nothing(galaxy, node_levels, monkeypatch,
                                                                                 caplog):
    """PERF.141 / Boss's no-recompile rule: a warm-up over a galaxy that has not changed builds nothing, and every
    phase says so."""
    import logging
    _place_addressed(galaxy.config, "Warm", (5.0, 1.0, 2.0))      # the built levels have a block to warm
    caplog.set_level(logging.INFO)
    _all_settings(monkeypatch, node_levels)
    monkeypatch.setattr(warm_cli, "_site_database", lambda: galaxy.config.database)
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    warm_cli.run_warm(_args())
    first = _phase_lines(caplog)
    assert [line.split(":")[0] for line in first] == [
        "Phase 1 of 5, opening", "Phase 2 of 5, built-levels", "Phase 3 of 5, deep-nodes", "Phase 4 of 5, objects",
        "Phase 5 of 5, search"]
    assert "built" in first[1] and "unchanged" not in first[1] and "unchanged" not in first[2]
    caplog.clear()
    warm_cli.run_warm(_args())
    again = _phase_lines(caplog)
    assert len(again) == 5
    # The sector page of a sector with a grid address asks for its Profile, which this galaxy (no plan) refuses with a 404 every
    # time: no answer is built, so the objects phase is unchanged like the others (MAP.289).
    assert all(": unchanged, nothing built" in line for line in again), again


def _plan_galaxy(config):
    """The galaxy's plan (its stored shape), which a sector's Profile is read from (`GET /api/galaxy/profile`)."""
    import math
    from planetgen import tuning
    from planetgen.galaxy.density import build_galaxy_shape
    shape = build_galaxy_shape(disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0,
                               bulge_amplitude=3.11, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4)
    store.save_galaxy_shape(shape, edge_pc=float(tuning.DEFAULT_SECTOR_EDGE_PC), outer_shell_index=1500,
                            expected_system_count_at_density_1=600.0, config=config, galaxy_seed=bytes(range(16)))


def test_a_planned_galaxy_with_addressed_sectors_warms_their_profiles_once_and_a_second_run_builds_nothing(
        galaxy, node_levels, monkeypatch, caplog):
    """MAP.289: the sector page reads its Profile, an answer the warm-up caches like any other, so a second run over the
    same galaxy finds the objects phase unchanged (the profile counted as an answer built every run before)."""
    import logging
    caplog.set_level(logging.INFO)
    _place_addressed(galaxy.config, "Warm", (5.0, 1.0, 2.0))
    _place_addressed(galaxy.config, "Warmer", (9.0, -3.0, 4.0))
    _plan_galaxy(galaxy.config)
    _all_settings(monkeypatch, node_levels)
    monkeypatch.setattr(warm_cli, "_site_database", lambda: galaxy.config.database)
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    warm_cli.run_warm(_args("objects"))
    first = _phase_lines(caplog)
    assert first and "unchanged" not in first[0], first
    caplog.clear()
    warm_cli.run_warm(_args("objects"))
    assert any(": unchanged, nothing built" in line for line in _phase_lines(caplog)), _phase_lines(caplog)


def test_a_change_in_one_sector_rebuilds_only_the_nodes_holding_it_in_the_phases_that_hold_them(
        galaxy, node_levels, monkeypatch, caplog):
    import logging
    caplog.set_level(logging.INFO)
    placed = _place_addressed(galaxy.config, "Placed", (5.0, 1.0, 2.0))
    monkeypatch.setattr(warm_cli, "_site_database", lambda: galaxy.config.database)
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    warm_cli.run_warm(_args("built-levels", "deep-nodes"))
    caplog.clear()
    time.sleep(0.01)
    _rename_sector(galaxy, placed)
    warm_cli.run_warm(_args("built-levels", "deep-nodes"))
    lines = {line.split(":")[0]: line for line in _phase_lines(caplog)}
    deep = lines["Phase 3 of 5, deep-nodes"]
    assert not deep.split(":")[1].strip().startswith("0 built"), lines


def test_a_run_stops_cleanly_at_a_small_budget_and_the_node_phases_say_so(galaxy, node_levels, monkeypatch, caplog):
    import logging
    caplog.set_level(logging.INFO)
    _all_settings(monkeypatch, node_levels)
    monkeypatch.setattr(warm_cli, "_site_database", lambda: galaxy.config.database)
    monkeypatch.setattr(warmup, "NODE_CHUNK", 4)
    monkeypatch.setenv("PLANETGEN_CACHE_MAX_MB", "0.0004")      # about 400 bytes: the first chunk fills it
    warm_cli.run_warm(_args("deep-nodes"))
    lines = {line.split(":")[0]: line for line in _phase_lines(caplog)}
    assert "left out" in lines["Phase 3 of 5, deep-nodes"] and lines["Phase 3 of 5, deep-nodes"].endswith(
        "stopped at the cache budget.")
    assert "Cache warm-up done" in caplog.text, "stopping at the budget is not a failure"


def test_the_phases_hold_every_family_once_in_the_order_they_run():
    assert [name for name, _text, _families in warmup.PHASES] == [
        "opening", "built-levels", "deep-nodes", "objects", "search"]
    assert [family for _name, _text, families in warmup.PHASES for family in families] == list(warmup.FAMILIES)
    assert list(commands.KIND_OF_WARM_FAMILY) == list(warmup.FAMILIES)
    assert warmup.phase_of("systems") == "objects" and warmup.phase_of("galaxy-deep-nodes") == "deep-nodes"


def test_a_phase_name_stands_for_all_its_families_in_the_order_of_the_run():
    assert warmup.expand_phases(["objects"]) == ["sectors", "systems", "phenomena"]
    assert warmup.expand_phases(["search", "lists", "deep-nodes", "lists"]) == [
        "lists", "galaxy-deep-nodes", "search"]
    steps = wrapper.steps_for("warm", _args("deep-nodes", "opening"))
    assert [step.kind for step in steps] == ["warm-cluster-field", "warm-object-index", "warm-galaxy-map", "warm-lists",
                                                  "warm-galaxy-deep-nodes"]
    with pytest.raises(ValueError, match="no such cache family or phase"):
        warmup.expand_phases(["phase-seven"])


def test_a_phase_that_built_nothing_says_unchanged_and_one_cut_short_says_how_many_were_left_out():
    tally = warmup.PhaseTally(["galaxy-map", "lists", "galaxy-nodes", "galaxy-deep-nodes", "search"])
    tally.add("galaxy-map", {"built": 0, "cached": 4, "warmed": 4, "seconds": 0.5})
    tally.add("lists", {"built": 0, "cached": 20, "warmed": 20, "seconds": 1.0})
    assert tally.start("lists") == []
    assert tally.start("galaxy-nodes") == ["Phase 1 of 5, opening: unchanged, nothing built (24 already stored) in 1.5 s."]
    tally.add("galaxy-nodes", {"built": 3, "cached": 125, "warmed": 128, "seconds": 2.0})
    tally.add("galaxy-deep-nodes", {"built": 0, "cached": 0, "warmed": 0, "left_out": 128, "seconds": 0.0})
    assert tally.start("galaxy-deep-nodes") == ["Phase 2 of 5, built-levels: 3 built, 125 already stored, 0 left out in 2.0 s."]
    tally.add("search", {"built": 0, "cached": 1, "warmed": 1, "seconds": 0.1})
    assert tally.start("search") == [
        "Phase 3 of 5, deep-nodes: 0 built, 0 already stored, 128 left out in 0.0 s; stopped at the cache budget."]
    assert tally.finish() == ["Phase 5 of 5, search: unchanged, nothing built (1 already stored) in 0.1 s."]
    assert tally.finish() == [], "each phase is reported once"


def test_the_deep_node_command_is_cut_into_chunks_and_the_map_then_finds_every_node_stored(galaxy, node_levels,
                                                                                         monkeypatch):
    expansion, results = _run_deep_units(galaxy)
    keys = [key for args, _cost in expansion.units for key in args["keys"]]
    assert len(keys) == expansion.plan["nodes"] == len(set(keys)) > 0
    assert all(0 < len(args["keys"]) <= warmup.NODE_CHUNK and cost > 0 for args, cost in expansion.units)
    assert sum(result["warmed"] for result in results) == len(keys)
    assert sum(result["built"] for result in results) == len(keys) and results[0]["noun"] == "nodes"
    # A second run finds every node stored and current: it builds nothing.
    _expansion, again = _run_deep_units(galaxy)
    assert sum(result["built"] for result in again) == 0 and sum(result["warmed"] for result in again) == len(keys)
    # A web worker that never saw the warm-up asks the API for none of them.
    asked = []
    real = tilecache.get_galaxy_nodes
    monkeypatch.setattr(tilecache, "get_galaxy_nodes", lambda db, wanted: asked.append(list(wanted)) or real(db, wanted))
    for start in range(0, len(keys), tilecache.MAX_NODES_PER_REQUEST):
        got = tilecache.fetch_nodes(galaxy.config.database, keys[start:start + tilecache.MAX_NODES_PER_REQUEST])
        assert got["cached"] == len(got["nodes"])
    assert [wanted for wanted in asked if wanted] == []


def test_the_views_the_map_asks_for_on_the_way_in_are_stored(galaxy, node_levels):
    with warmup.build_app().test_request_context("/"):
        shape = apiclient.get_galaxy_shape(galaxy.config.database)
    opening = warmup.opening_request(shape)
    assert opening and len(opening) <= tilecache.MAX_NODES_PER_REQUEST
    with warmup.build_app().test_request_context("/"):
        warmed = warmup.warm_opening_view(galaxy.config.database)
    assert warmed["nodes"] == len(opening) and warmed["cached"] == 0


# --- Building around a visit (PERF.132, MAP.260) ---------------------------------------------------

ASKED = "81.2.1.1"
"""str: A node the prefetch tests visit; 60 nodes lie around it (its children, its parent and the parent's children)."""


@pytest.fixture
def prefetch(galaxy, monkeypatch):
    """The background build on the seeded galaxy, with the nodes it queues built in this process on a queue of its
    own, and the limits as told."""
    from planetgen.web.lib import tileprefetch
    real = warmup.build_app
    monkeypatch.setattr(warmup, "build_app", lambda config=None: real(galaxy.config))
    real(galaxy.config)                 # a web process: the app installs the in-process transport the node cache reads by

    class Rules:
        prefetch_per_visit, prefetch_per_minute, prefetch_max_queued = 20, 100, 8
    monkeypatch.setattr(tileprefetch, "get_settings", lambda: types.SimpleNamespace(tile_cache=Rules))
    ran = []
    monkeypatch.setattr(tileprefetch.workers, "ensure_workers", lambda *args, **kwargs: ran.append(kwargs) or 0)
    token = uuid.uuid4().hex[:12]           # this test's own claims and minute counters on the shared Redis
    monkeypatch.setattr(tileprefetch, "CLAIM_KEY", f"planetgen-test:{token}:claim:")
    monkeypatch.setattr(tileprefetch, "MINUTE_KEY", f"planetgen-test:{token}:minute:")
    Rules.module, Rules.ran, Rules.db = tileprefetch, ran, galaxy.config.database

    def schedule(keys):
        with real(galaxy.config).test_request_context("/"):
            return tileprefetch.schedule(galaxy.config.database, keys)
    Rules.schedule = staticmethod(schedule)
    yield Rules
    client = redisclient.client(decode_responses=True)
    for key in client.scan_iter(f"planetgen-test:{token}:*"):
        client.delete(key)


def _run_queued(prefetch):
    import rq
    connection = redisqueue.connect()
    done = rq.SimpleWorker([redisqueue.queue(redisqueue.MAINTENANCE_QUEUE, connection)], connection=connection,
                           serializer=redisqueue.JsonSerializer)
    done.work(burst=True, logging_level="WARNING")


def test_a_visit_queues_the_nodes_around_it_and_they_are_then_stored(galaxy, prefetch):
    with warmup.build_app().test_request_context("/"):
        tilecache.fetch_nodes(prefetch.db, [ASKED])          # the visit: the map asked for this node
    queued = prefetch.schedule([ASKED])
    assert queued == 20 and prefetch.ran and "count" not in prefetch.ran[0], "the maintenance lane sets its own cap"
    _run_queued(prefetch)
    around = list(prefetch.module.candidates([ASKED]))[:20]
    with warmup.build_app().test_request_context("/"):
        stored = tilecache.fetch_nodes(prefetch.db, around)
    assert stored["cached"] == 20, "a visitor asking for them next is served from the disk"


def test_a_visit_to_a_built_area_queues_nothing_and_a_second_visitor_causes_no_build(galaxy, prefetch):
    assert prefetch.schedule([ASKED]) == 20
    assert prefetch.schedule([ASKED]) == 0, "the same area is claimed already"
    _run_queued(prefetch)
    assert prefetch.schedule([ASKED]) == 20, "the next 20 of the area's 60 nodes, not the 20 built"
    _run_queued(prefetch)
    assert prefetch.schedule([ASKED]) == 20
    _run_queued(prefetch)
    every = list(prefetch.module.candidates([ASKED]))
    with warmup.build_app().test_request_context("/"):
        assert tilecache.fetch_nodes(prefetch.db, every)["cached"] == len(every) == 60
    assert prefetch.schedule([ASKED]) == 0, "built: nothing left to queue"


def test_the_background_build_keeps_to_its_limits(galaxy, prefetch):
    prefetch.prefetch_per_visit, prefetch.prefetch_per_minute = 50, 30
    assert prefetch.schedule([ASKED]) == 30, "the minute's cap trims the visit"
    assert prefetch.schedule(["243.1.1.1"]) == 0, "and the minute is spent"
    prefetch.prefetch_per_visit = 0
    assert prefetch.schedule(["81.3.4.15"]) == 0, "off"


def test_the_background_build_waits_while_the_cache_is_full_or_the_queue_busy(galaxy, prefetch):
    root = tilecache.cache_dir()
    prefetch.module.mark_full(root)
    assert prefetch.module.is_full(root) and prefetch.schedule([ASKED]) == 0
    import os
    os.unlink(os.path.join(root, prefetch.module.FULL_FILE))
    prefetch.prefetch_max_queued = 1
    assert prefetch.schedule([ASKED]) == 20      # one job per 16 nodes: two wait now
    assert prefetch.schedule(["243.1.1.1"]) == 0, "bulk work is waiting"


def test_an_answered_request_for_nodes_queues_the_nodes_around_it(galaxy, prefetch, monkeypatch):
    asked = []
    monkeypatch.setattr(prefetch.module, "schedule", lambda db, keys: asked.append((db, list(keys))) or 0)
    response = _web_worker(galaxy.config).get(f"/galaxy/nodes?keys={ASKED},243.1.1.1")
    assert response.status_code == 200
    response.close()
    assert asked == [(galaxy.config.database, [ASKED, "243.1.1.1"])]


def test_a_change_builds_the_stored_nodes_it_touched_again_at_once(galaxy, prefetch, monkeypatch):
    """PERF.133: the nodes a change deleted are queued for a rebuild, so the next visit finds them built."""
    placed = _place_addressed(galaxy.config, "Placed", (5.0, 1.0, 2.0))
    keys = warmup.deep_nodes_to_warm(None, 81)[:60]
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    app = warmup.build_app()
    with app.test_request_context("/"):
        tilecache.fetch_nodes(prefetch.db, keys)
    assert tilecache.take_rebuilds(prefetch.db) == (None, []), "nothing changed, nothing is owed"
    time.sleep(0.01)
    _rename_sector(galaxy, placed)
    with app.test_request_context("/"):
        tilecache.fetch_nodes(prefetch.db, keys[:1])             # the request that notices the change
    generation, owed = tilecache.take_rebuilds(prefetch.db)
    assert generation and owed and set(owed) <= set(keys), "every stored node that held the sector"
    assert tilecache.take_rebuilds(prefetch.db) == (None, []), "owed once"
    with app.test_request_context("/"):
        assert prefetch.module.rebuild(prefetch.db, generation, owed) == len(owed)
        assert prefetch.module.rebuild(prefetch.db, generation, owed) == 0, "a node already queued is not queued twice"
    _run_queued(prefetch)
    with app.test_request_context("/"):
        assert tilecache.fetch_nodes(prefetch.db, owed)["cached"] == len(owed)


def test_the_rebuild_of_a_change_is_off_with_the_prefetch_and_during_an_update(galaxy, prefetch, monkeypatch):
    keys = [ASKED, "243.1.1.1"]
    with warmup.build_app().test_request_context("/"):
        generation = tilecache.current_stamp(prefetch.db, tilecache.cache_dir())["generation"]
        monkeypatch.setattr(prefetch, "prefetch_per_minute", 0)
        assert prefetch.module.rebuild(prefetch.db, generation, keys) == 0
        monkeypatch.setattr(prefetch, "prefetch_per_minute", 100)
        monkeypatch.setattr(prefetch.module.updateguard, "is_updating", lambda: True)
        assert prefetch.module.rebuild(prefetch.db, generation, keys) == 0


def test_an_answered_request_builds_what_a_change_touched_again(galaxy, prefetch, monkeypatch):
    rebuilt = []
    monkeypatch.setattr(prefetch.module, "schedule", lambda db, keys: 0)
    monkeypatch.setattr(prefetch.module, "rebuild", lambda db, generation, keys: rebuilt.append((db, generation, keys)))
    tilecache._owe_rebuild(galaxy.config.database, "00000000000000aa", [ASKED, "243.1.1.1"])
    response = _web_worker(galaxy.config).get(f"/galaxy/nodes?keys={ASKED}")
    response.close()
    assert rebuilt == [(galaxy.config.database, "00000000000000aa", [ASKED, "243.1.1.1"])]


def test_prune_marks_the_cache_full_for_the_background_build(tmp_path, monkeypatch):
    from planetgen.web.lib import tileprefetch
    (tmp_path / "t1.bin").write_bytes(b"x" * 100)
    (tmp_path / "t2.bin").write_bytes(b"x" * 100)
    tilecache.prune(str(tmp_path), max_bytes=500)
    assert not tileprefetch.is_full(str(tmp_path))
    tilecache.prune(str(tmp_path), max_bytes=100)
    assert tileprefetch.is_full(str(tmp_path))


def _args(*argv):
    parser, _parsers = generate_cli.build_parser()
    args = parser.parse_args(["warm", *argv])
    args.database = "x"
    return args


def test_the_list_is_every_family_in_order_and_then_the_command_that_opens_the_site():
    steps = wrapper.steps_for("warm", _args())
    assert [step.kind for step in steps] == [f"warm-{family}" for family in warmup.FAMILIES] + [commands.END_MAP_REBUILD]
    assert [step.key for step in steps][:-1] == list(warmup.FAMILIES)
    assert all(step.params == {"database": "x", "cap": None} for step in steps[:-1])
    assert steps == wrapper.steps_for("warm", _args()), "the same arguments give the same list"


def test_a_part_of_the_families_is_in_the_order_of_the_families_and_does_not_open_the_site():
    steps = wrapper.steps_for("warm", _args("search", "lists", "--cap", "7"))
    assert [step.kind for step in steps] == ["warm-lists", "warm-search"]
    assert all(step.params == {"database": "x", "cap": 7} for step in steps)
    with pytest.raises(ValueError, match="no such cache family or phase 'nonsense'"):
        wrapper.steps_for("warm", _args("nonsense"))


# --- `planetgen warm` ------------------------------------------------------------------------------


@pytest.fixture
def state_file(tmp_path, monkeypatch):
    path = tmp_path / "maintenance.json"
    monkeypatch.setenv(maintenance.FILE_ENV, str(path))
    return path


@pytest.fixture
def script(monkeypatch, redis_server):
    """`planetgen warm` with the site's database named and the families that succeed or fail as told."""
    class Plan:
        failing = {}
        ran = []
        owners = []

    def warm_family(db, family, cap=None, app=None, config=None):
        conn = redisqueue.connect(redis_server)
        Plan.owners.append((runrecord.active_run(conn), runrecord.active_run(conn, runrecord.MAINTENANCE)))
        Plan.ran.append((family, cap, db))
        if family in Plan.failing:
            raise warmup.WarmupFailed(Plan.failing[family])
        return {"family": family, "warmed": 3, "cached": 3, "built": 0, "left_out": 0, "seconds": 0.1,
                "line": f"{family}: 3 pages warmed (0 answers built, 3 already cached, 0 left out by the cap) in 0.1 s"}

    monkeypatch.setattr(warm_cli, "_site_database", lambda: "x")
    monkeypatch.setattr(warmup, "warm_family", warm_family)
    # The node commands are many units on a real galaxy; here, like the other families, each is one unit that is faked.
    monkeypatch.setattr(warmup, "expand_warm_galaxy_nodes", lambda params: commands.Expansion(
        [({"database": params["database"], "cap": params.get("cap")}, 1.0)], {}))
    monkeypatch.setattr(warmup, "warm_galaxy_nodes_unit", lambda args, report: warmup._warm_unit("galaxy-nodes", args))
    monkeypatch.setattr(warmup, "expand_warm_galaxy_deep_nodes", lambda params: commands.Expansion(
        [({"database": params["database"], "cap": params.get("cap")}, 1.0)], {}))
    monkeypatch.setattr(warmup, "warm_galaxy_deep_nodes_unit",
                        lambda args, report: warmup._warm_unit("galaxy-deep-nodes", args))
    # The object families are slices of a listing the site answers; here, like the node ones, each is one faked unit.
    for family in ("sectors", "systems", "phenomena"):
        monkeypatch.setattr(warmup, f"expand_warm_{family}", lambda params: commands.Expansion(
            [({"database": params["database"], "cap": params.get("cap")}, 1.0)], {}))
    # The object index is built on the database itself, not through the site; here it is one more faked family.
    monkeypatch.setattr(warmup, "warm_object_index_unit", lambda args, report: warmup._warm_unit("object-index", args))
    # The cluster field is made from the database's own rows (PERF.225); here it is one more faked family.
    from planetgen.galaxy import field_job
    monkeypatch.setattr(field_job, "expand_field", lambda params: commands.Expansion(
        [({"database": params["database"], "cap": params.get("cap")}, 1.0)], {}))
    monkeypatch.setattr(field_job, "field_unit", lambda args, report: warmup._warm_unit("cluster-field", args))
    Plan.ran, Plan.failing, Plan.owners = [], {}, []
    return Plan


def test_the_script_prints_one_line_per_family(script, capsys):
    warm_cli.run_warm(_args())
    out = capsys.readouterr().out
    assert [family for family, _cap, _db in script.ran] == list(warmup.FAMILIES)
    assert {db for _family, _cap, db in script.ran} == {"x"}
    for family in warmup.FAMILIES:
        assert f"{family}: 3 pages warmed" in out
    assert "Maintenance mode was not the map rebuild's" in out and f"Cache warm-up done: {len(warmup.FAMILIES)} families." in out


def test_a_failing_command_stops_the_list_exits_non_zero_and_is_named(script, caplog):
    script.failing = {"systems": "2 of 9 pages did not answer"}
    with pytest.raises(SystemExit) as exit_info:
        warm_cli.run_warm(_args())
    assert exit_info.value.code == 1
    err = caplog.text
    assert "the cache warm-up stopped: unit 0 of stage Warm systems: every system page failed: WarmupFailed: 2 of 9 pages did not answer" in err
    assert [family for family, _cap, _db in script.ran] == ["cluster-field", "object-index", "galaxy-map", "lists",
                                                            "galaxy-nodes", "galaxy-deep-nodes", "sectors", "systems"], \
        "the wrapper's list stops at the failing command"


def test_the_families_and_the_cap_are_passed_on(script):
    warm_cli.run_warm(_args("sectors", "lists", "--cap", "7"))
    assert [(family, cap) for family, cap, _db in script.ran] == [("lists", 7), ("sectors", 7)]


def test_an_unknown_family_is_an_error_before_anything_runs(script):
    with pytest.raises(SystemExit, match="no such cache family or phase 'nonsense'"):
        warm_cli.run_warm(_args("nonsense"))
    assert script.ran == []


def test_a_part_of_the_warm_up_does_not_reopen_the_site(script, state_file, capsys):
    maintenance.turn_on(by="migration", cause=maintenance.SECTOR_SWITCH)
    warm_cli.run_warm(_args("lists"))
    assert maintenance.is_on()
    warm_cli.run_warm(_args())
    assert not maintenance.is_on()
    assert "Maintenance mode is off: the Galaxy Map is rebuilt." in capsys.readouterr().out


def test_a_failed_warm_up_leaves_the_site_closed(script, state_file, capsys):
    maintenance.turn_on(by="migration", cause=maintenance.SECTOR_SWITCH)
    script.failing = {"search": "x"}
    with pytest.raises(SystemExit):
        warm_cli.run_warm(_args())
    assert maintenance.is_on()


def test_the_warm_up_is_a_maintenance_run_that_leaves_the_machine_free(script, redis_server):
    conn = redisqueue.connect(redis_server)
    generating = runrecord.start_run(conn, "A generation run", [runrecord.Stage("a", "A")])
    warm_cli.run_warm(_args("lists"))
    assert script.owners == [(generating, script.owners[0][1])] and script.owners[0][1] is not None, \
        "the generation run keeps the machine and the warm-up has a place of its own"
    [warm_id] = [run_id for run_id in runrecord.recent_runs(conn) if run_id != generating]
    run = runrecord.read_run(conn, warm_id)
    assert run["kind"] == runrecord.MAINTENANCE and run["state"] == "done"
    assert runrecord.active_run(conn) == generating and runrecord.active_run(conn, runrecord.MAINTENANCE) is None
    from planetgen.web.lib import galaxystatus
    assert galaxystatus.compute(conn, "x")["active"] is True, "the banner is the generation run's"
    runrecord.finish_run(conn, generating)
    assert galaxystatus.compute(conn, "x")["active"] is False


def test_a_second_warm_up_while_one_runs_is_refused(script, redis_server, caplog):
    conn = redisqueue.connect(redis_server)
    runrecord.start_run(conn, "A warm-up", [runrecord.Stage("a", "A")], kind=runrecord.MAINTENANCE)
    with pytest.raises(SystemExit) as exit_info:
        warm_cli.run_warm(_args("lists"))
    assert exit_info.value.code == 1 and "is running" in caplog.text
    assert script.ran == []


def test_no_redis_is_an_error_not_a_silent_skip(script, monkeypatch, caplog):
    monkeypatch.setattr(redisqueue, "redis_url", lambda: "redis://127.0.0.1:1/0")
    with pytest.raises(SystemExit) as exit_info:
        warm_cli.run_warm(_args("lists"))
    assert exit_info.value.code == 1
    assert "no Redis server answers" in caplog.text and script.ran == []


def test_the_script_waits_for_an_update_and_gives_up_after_the_wait(monkeypatch):
    states = iter([True, True, False])
    monkeypatch.setattr(warm_cli.updateguard, "is_updating", lambda: next(states))
    monkeypatch.setattr(warm_cli.time, "sleep", lambda seconds: None)
    warm_cli._wait_for_update()
    monkeypatch.setattr(warm_cli.updateguard, "is_updating", lambda: True)
    monkeypatch.setattr(warm_cli, "UPDATE_WAIT_SECONDS", -1)
    with pytest.raises(warm_cli.Failed, match="still running"):
        warm_cli._wait_for_update()


def test_the_warm_up_the_update_runs_does_not_wait_for_the_update(monkeypatch):
    """update.sh runs `planetgen warm` itself, before the reload: waiting for the update it belongs to would never end."""
    monkeypatch.setattr(warm_cli.updateguard, "is_updating", lambda: True)
    monkeypatch.setattr(warm_cli.updateguard, "inside_update", lambda: True)
    monkeypatch.setattr(warm_cli, "UPDATE_WAIT_SECONDS", -1)
    warm_cli._wait_for_update()


def test_the_progress_bar_counts_families():
    """update.sh shows this line while it waits for the warm-up."""
    assert warm_progress.text_bar(1, 4, "lists", width=8) == "[##------] 1/4 lists"
    assert warm_progress.text_bar(4, 4, "search", width=8) == "[########] 4/4 search"


def test_a_sigterm_ends_the_run_as_cancelled_not_stranded(script, redis_server, monkeypatch):
    """SIGTERM (the update guard's cancel, update.sh's timeout) reaches the starter as Ctrl+C does, and the run ends
    as cancelled in Redis, freeing its place."""
    import os
    import signal
    original = signal.getsignal(signal.SIGTERM)
    warm_cli._stop_on_sigterm()
    try:
        with pytest.raises(KeyboardInterrupt):
            os.kill(os.getpid(), signal.SIGTERM)
    finally:
        signal.signal(signal.SIGTERM, original)

    def interrupted(*args, **kwargs):
        raise KeyboardInterrupt

    monkeypatch.setattr(units, "follow_stage", interrupted)
    with pytest.raises(KeyboardInterrupt):
        warm_cli.run_warm(_args("lists"))
    conn = redisqueue.connect(redis_server)
    assert runrecord.read_run(conn, runrecord.recent_runs(conn)[0])["state"] == "cancelled"
    assert runrecord.active_run(conn, runrecord.MAINTENANCE) is None


# --- The Galaxy Map's sector-block nodes (MAP.260) --------------------------------------------------


def _make_block_rows(config):
    """The block aggregates of the seeded galaxy: a few pre-placed stars (the seeded sectors have no grid address, so
    the stored light of the stars is what makes blocks), the totals `generate` makes after a run."""
    from planetgen.db import block_lights, block_stats
    from tests.test_galaxy_nodes import _bright_row
    conn = store.get_connection(config)
    try:
        store.insert_bright_stars(conn, [_bright_row(position, 1e4) for position in
                                         ((900.0, 20.0, 5.0), (-700.0, 300.0, 40.0), (10.0, 10.0, 1.0))])
        block_stats.rebuild(conn)
        block_lights.rebuild(conn)
        conn.commit()
    finally:
        conn.close()


def _run_node_units(galaxy):
    expansion = commands.expand_function("warm-galaxy-nodes")({"database": galaxy.config.database, "cap": None})
    return expansion, [commands.run_alone("warm-galaxy-nodes", args) for args, _cost in expansion.units]


def test_the_nodes_to_warm_are_the_populated_blocks_from_the_top_down_and_a_second_run_builds_none(
        galaxy, node_levels, monkeypatch):
    from planetgen.db import nodes
    from planetgen.galaxy.drill import parse_node_key
    _make_block_rows(galaxy.config)
    conn = store.get_connection(galaxy.config)
    try:
        keys = nodes.populated_node_keys(conn, warmup.NODE_WARM_FLOOR)
    finally:
        conn.close()
    assert keys, "the seeded sector's blocks are stored"
    levels = [parse_node_key(key).m for key in keys]
    assert levels == sorted(levels, reverse=True) and min(levels) >= warmup.NODE_WARM_FLOOR
    assert len(keys) == len(set(keys))
    monkeypatch.setattr(tilecache, "STAMP_TTL_SECONDS", 0)
    expansion, first = _run_node_units(galaxy)
    assert expansion.units and expansion.plan["nodes"] == len(keys)
    assert sum(result["warmed"] for result in first) == len(keys)
    assert sum(result["built"] for result in first) == len(keys) and all(result["left_out"] == 0 for result in first)
    _expansion, again = _run_node_units(galaxy)
    assert sum(result["built"] for result in again) == 0
    assert sum(result["cached"] for result in again) == len(keys)
    assert "already stored" in again[0]["line"]


def test_a_full_cache_leaves_the_nodes_out_and_says_so(galaxy, node_levels, monkeypatch):
    _make_block_rows(galaxy.config)
    monkeypatch.setattr(tilecache, "cache_is_full", lambda *args, **kwargs: True)
    _expansion, results = _run_node_units(galaxy)
    assert sum(result["warmed"] for result in results) == 0
    assert sum(result["left_out"] for result in results) > 0 and "left out" in results[0]["line"]


def test_the_nodes_family_warms_inline_and_belongs_to_the_built_levels_phase(galaxy, node_levels):
    _make_block_rows(galaxy.config)
    assert warmup.phase_of("galaxy-nodes") == "built-levels"
    result = warmup.warm_family(galaxy.config.database, "galaxy-nodes", config=galaxy.config)
    assert result["warmed"] > 0 and result["noun"] == "nodes"


def test_the_nodes_after_a_visit_are_its_children_then_its_parent_and_the_parents_children():
    from planetgen.galaxy.drill import drill_children, drill_parent, format_drill_key, parse_node_key
    from planetgen.web.lib import tileprefetch
    asked = ["243.30.4.2", "243.31.4.2"]
    found = list(tileprefetch.candidates(asked))
    assert len(found) == len(set(found)) and not set(found) & set(asked)
    first = parse_node_key(asked[0])
    children = [format_drill_key(child) for _shell, group in drill_children(first) for child in group]
    assert found[:len(children)] == children, "the children of the first (nearest) node come first"
    parent = format_drill_key(drill_parent(first))
    assert parent in found
    assert list(tileprefetch.candidates(["not-a-key"])) == []
    sector = parse_node_key("3.30.4.2")
    assert all(parse_node_key(key).m in (1, 3, 9) for key in tileprefetch.candidates([format_drill_key(sector)]))




def _plan_and_place(galaxy, count):
    """A planned galaxy and `count` more sectors with a grid address (the sector page shows their Profile)."""
    import math
    from planetgen import tuning
    from planetgen.galaxy.density import build_galaxy_shape
    shape = build_galaxy_shape(disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0,
                               bulge_amplitude=3.11, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4)
    store.save_galaxy_shape(shape, edge_pc=float(tuning.DEFAULT_SECTOR_EDGE_PC), outer_shell_index=1500,
                            expected_system_count_at_density_1=600.0, config=galaxy.config, galaxy_seed=bytes(range(16)))
    edge_pc = ly_to_pc(DEFAULT_SECTOR_EDGE_LY)
    for number in range(count):
        center = (5.0 + 9.0 * number, 5.0, 5.0)
        _place_sector(galaxy.config, f"Addressed {number}", center, address=spherical.sector_address_at(center, edge_pc))


def test_the_warm_up_counts_the_profile_answers_of_addressed_sectors_as_cached_on_a_second_run(galaxy):
    """MAP.289: the Profile each sector page asks for (`/galaxy/profile`) is one more cached answer, so a second run
    over the same sectors builds none."""
    _plan_and_place(galaxy, 3)
    first = _warm(galaxy, "sectors")
    second = _warm(galaxy, "sectors")
    assert first["warmed"] == 5
    assert second["built"] == 0, second["line"]
