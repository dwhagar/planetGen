# tests/test_cache_warmup.py

"""
PERF.104, PERF.127: the warm-up (`planetgen/web/warmup.py`, `planetgen warm`) fills the shared
caches so the first visit answers from them, one command per family, loudly.

The families visit the pages of a real app on a throwaway database; a second, separate app
(a web worker that never saw the warm-up) is then asked for the same pages with a spy on the API
calls, to show nothing is rebuilt. The wrapper's list of the commands runs on the real Redis queue (the
tests' workers run the units here), with the families themselves faked where only the list is under test.
"""

import contextlib
import types
import uuid

import pytest

from planetgen.admin import maintenance
from planetgen.api.config import Config
from planetgen.cli import generate as generate_cli
from planetgen.cli import warm as warm_cli
from planetgen.db import dataversion, store
from planetgen.galaxy import viewport
from planetgen.galaxy.sector import SpaceSector
from planetgen.queue import commands, redisclient, redisqueue, runrecord, units, wrapper
from planetgen.web import warmup
from planetgen.web.app import create_app
from planetgen.web.lib import apiclient, tilecache
from tests.publicids import pid
from tests.test_api import seeded_sector  # noqa: F401 -- the fixture


@pytest.fixture
def galaxy(seeded_sector, redis_server, tmp_path, monkeypatch):  # noqa: F811
    """The seeded sector (two systems) and an empty one, a nebula in the first, and a tile cache of its own."""
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
    assert lists["warmed"] == 6 and search["warmed"] == 1
    web = _web_worker(galaxy.config)
    del spy[:]
    for path in ("/sectors", "/systems", "/phenomena", "/", "/search"):
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
    assert result["noun"] == "tiles" and result["warmed"] > 0
    web = _web_worker(galaxy.config)
    del spy[:]
    assert web.get("/galaxy").status_code == 200
    assert [call for call in spy if call.startswith(("/galaxy/sectors", "/galaxy/shape", "/galaxy/phenomena"))] == []


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
    with pytest.raises(ValueError, match="cannot be cut"):
        commands.expand_function("run-population")


def test_every_cache_family_is_its_own_registered_command_with_one_unit():
    assert list(commands.KIND_OF_WARM_FAMILY) == list(warmup.FAMILIES)
    for family, kind in commands.KIND_OF_WARM_FAMILY.items():
        assert not commands.spec(kind).inline
        if family == "galaxy-tiles":
            continue   # many units, one a chunk of tiles: its own tests below
        expansion = commands.expand_function(kind)({"database": "x", "cap": 4})
        assert expansion.units == [({"database": "x", "cap": 4}, 1.0)], "one unit that warms that family and no other"
    assert commands.spec(commands.END_MAP_REBUILD).summary


# --- The Galaxy Map's tiles (PERF.130) -----------------------------------------------------------

@pytest.fixture
def tile_levels(galaxy, monkeypatch):
    """`warm-galaxy-tiles` on the seeded galaxy with a small octree (levels 0 to 3 in full, a zoom run to level 5)."""
    real = warmup.build_app
    monkeypatch.setattr(warmup, "build_app", lambda config=None: real(galaxy.config))
    levels = types.SimpleNamespace(warm_levels=3, warm_zoom_levels=5)
    monkeypatch.setattr(warmup, "get_settings", lambda: types.SimpleNamespace(tile_cache=levels))
    return levels


def _run_tile_units(galaxy):
    expansion = commands.expand_function("warm-galaxy-tiles")({"database": galaxy.config.database, "cap": None})
    return expansion, [commands.run_alone("warm-galaxy-tiles", args) for args, _cost in expansion.units]


def test_the_first_levels_are_every_tile_coarsest_level_first_and_the_zoom_run_adds_finer_ones():
    bounds = ((-20000.0, -20000.0, -20000.0), (20000.0, 20000.0, 20000.0))
    keys = warmup.octree_tiles(bounds, 3)
    levels = [viewport.parse_tile_key(key)[0] for key in keys]
    assert levels == sorted(levels) and [levels.count(level) for level in range(4)] == [1, 8, 8, 8]
    assert len(set(keys)) == len(keys)
    run = warmup.zoom_run_tiles(bounds, 30000.0, 8)
    assert run and all(viewport.parse_tile_key(key)[0] <= 8 for key in run) and len(set(run)) == len(run)
    assert max(viewport.parse_tile_key(key)[0] for key in run) == 8, "the run goes in to the finest level asked for"
    shape = None
    every = warmup.tiles_to_warm(shape, 3, 8)
    assert every[:len(keys)] == warmup.octree_tiles(warmup.tile_bounds_pc(warmup.ly_to_pc(
        warmup.DEFAULT_SECTOR_EDGE_LY), shape), 3)
    assert warmup.tiles_to_warm(shape, 3, 3) == warmup.tiles_to_warm(shape, 3, 0), "a run no finer adds nothing"


def test_levels_that_would_build_too_many_tiles_are_refused_with_the_count():
    with pytest.raises(warmup.WarmupFailed, match=r"warm_levels \d+ would build \d+ tiles.*lower it"):
        warmup.tiles_to_warm(None, viewport.TILE_MAX_LEVEL, 0)


def test_the_tile_command_is_cut_into_chunks_and_the_map_then_finds_every_tile_stored(galaxy, tile_levels,
                                                                                   monkeypatch):
    expansion, results = _run_tile_units(galaxy)
    keys = [key for args, _cost in expansion.units for key in args["keys"]]
    assert len(keys) == expansion.plan["tiles"] == len(set(keys)) > 0
    assert all(0 < len(args["keys"]) <= warmup.TILE_CHUNK and cost > 0 for args, cost in expansion.units)
    assert sum(result["warmed"] for result in results) == len(keys)
    assert sum(result["built"] for result in results) == len(keys) and results[0]["noun"] == "tiles"
    # A second run finds every tile stored and current: it builds nothing.
    _expansion, again = _run_tile_units(galaxy)
    assert sum(result["built"] for result in again) == 0 and sum(result["warmed"] for result in again) == len(keys)
    # A web worker that never saw the warm-up asks the API for none of them.
    asked = []
    real = tilecache.get_galaxy_tiles
    monkeypatch.setattr(tilecache, "get_galaxy_tiles", lambda db, wanted: asked.append(list(wanted)) or real(db, wanted))
    for start in range(0, len(keys), tilecache.MAX_TILES_PER_REQUEST):
        got = tilecache.fetch_tiles(galaxy.config.database, keys[start:start + tilecache.MAX_TILES_PER_REQUEST])
        assert got["cached"] == len(got["tiles"])
    assert [wanted for wanted in asked if wanted] == []


def test_the_views_the_map_asks_for_on_the_way_in_are_stored(galaxy, tile_levels):
    _run_tile_units(galaxy)
    with warmup.build_app().test_request_context("/"):
        shape = apiclient.get_galaxy_shape(galaxy.config.database)
    stored = set(warmup.tiles_to_warm(shape, tile_levels.warm_levels, tile_levels.warm_zoom_levels))
    opening = [key for key in warmup.opening_request(shape)
               if viewport.parse_tile_key(key)[0] <= tile_levels.warm_zoom_levels]
    assert opening and set(opening) <= stored


# --- Building around a visit (PERF.132) ------------------------------------------------------------

def test_the_tiles_after_a_visit_are_its_neighbours_and_the_next_level_inside_it():
    from planetgen.web.lib import tileprefetch
    asked = ["5/16/16/16", "5/17/16/16"]
    found = list(tileprefetch.candidates(asked))
    assert len(found) == len(set(found)) and not set(found) & set(asked)
    assert found[:8] == [f"6/{32 + dx}/{32 + dy}/{32 + dz}" for dx in (0, 1) for dy in (0, 1) for dz in (0, 1)], \
        "the children of the first (nearest) tile come first"
    assert {"5/15/15/15", "5/17/17/17", "5/18/16/16"} <= set(found)
    corner = list(tileprefetch.candidates(["1/0/0/0"]))
    assert len([key for key in corner if key.startswith("1/")]) == 7, "the neighbours stop at the root cube"
    assert list(tileprefetch.candidates([f"{viewport.TILE_MAX_LEVEL}/0/0/0"])) == [
        f"{viewport.TILE_MAX_LEVEL}/{dx}/{dy}/{dz}" for dx in (0, 1) for dy in (0, 1) for dz in (0, 1)][1:], \
        "the finest tiles have no children"


def test_the_next_zoom_steps_are_the_cuts_for_the_camera_closer_along_its_line():
    from planetgen.web.lib import tileprefetch
    bounds = ((-20000.0, -20000.0, -20000.0), (20000.0, 20000.0, 20000.0))
    camera, look = (3000.0, 1000.0, 500.0), (200.0, 100.0, 0.0)
    found = list(tileprefetch.zoom_candidates(camera, look, bounds, steps=2))
    scale = viewport.DISTANCE_CUT_FADE_BAND ** -1
    first = viewport.fetch_cut_tiles(tuple(t + (c - t) * scale for c, t in zip(camera, look)), look, bounds)
    assert found[:len(first)] == first, "the nearest step comes first"
    assert len(found) == len(set(found)) > len(first)
    asked = ["5/16/16/16"]
    order = list(tileprefetch.wanted(asked, camera + look, bounds))
    assert order[:len(first)] == first and not set(order) & set(asked)
    assert list(tileprefetch.wanted(asked, None, None)) == list(tileprefetch.candidates(asked)), "no camera: the neighbours"


@pytest.mark.parametrize("text,expected", [
    ("1,2,3,4,5,6", (1.0, 2.0, 3.0, 4.0, 5.0, 6.0)), ("-1.5,0,3e3,4,5,6", (-1.5, 0.0, 3000.0, 4.0, 5.0, 6.0)),
    (None, None), ("", None), ("1,2,3", None), ("1,2,3,4,5,x", None), ("1,2,3,4,5,nan", None),
    ("1,2,3,4,5,inf", None), ("1,2,3,4,5,1e9", None)])
def test_the_cameras_of_a_request_are_six_plain_numbers_or_nothing(text, expected):
    from planetgen.web.galaxy_views import _camera_asked
    assert _camera_asked(text) == expected


@pytest.fixture
def prefetch(galaxy, monkeypatch):
    """The background build on the seeded galaxy, with the tiles it queues built in this process on a queue of its
    own, and the limits as told."""
    from planetgen.web.lib import tileprefetch
    real = warmup.build_app
    monkeypatch.setattr(warmup, "build_app", lambda config=None: real(galaxy.config))
    real(galaxy.config)                 # a web process: the app installs the in-process transport the tile cache reads by

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
    done = rq.SimpleWorker([redisqueue.queue(redisqueue.BULK_QUEUE, connection)], connection=connection,
                           serializer=redisqueue.JsonSerializer)
    done.work(burst=True, logging_level="WARNING")


def test_a_visit_queues_the_tiles_around_it_and_they_are_then_stored(galaxy, prefetch):
    with warmup.build_app().test_request_context("/"):
        tilecache.fetch_tiles(prefetch.db, ["5/16/16/16"])          # the visit: the map asked for this tile
    queued = prefetch.schedule(["5/16/16/16"])
    assert queued == 20 and prefetch.ran and prefetch.ran[0].get("count") == 1, "one worker, never the admin's count"
    _run_queued(prefetch)
    around = list(prefetch.module.candidates(["5/16/16/16"]))[:20]
    with warmup.build_app().test_request_context("/"):
        stored = tilecache.fetch_tiles(prefetch.db, around)
    assert stored["cached"] == 20, "a visitor asking for them next is served from the disk"


def test_a_visit_with_a_camera_queues_the_next_zoom_steps_first(galaxy, prefetch):
    camera, look = (9000.0, 3000.0, 800.0), (0.0, 0.0, 0.0)
    with warmup.build_app().test_request_context("/"):
        asked = viewport.fetch_cut_tiles(camera, look, ((-15000.0,) * 3, (15000.0,) * 3))
        queued = prefetch.module.schedule(prefetch.db, asked, camera + look)
        bounds = prefetch.module._bounds(prefetch.db)
    assert queued == 20
    ahead = [key for key in prefetch.module.wanted(asked, camera + look, bounds)][:20]
    _run_queued(prefetch)
    with warmup.build_app().test_request_context("/"):
        assert tilecache.fetch_tiles(prefetch.db, ahead)["cached"] == 20


def test_a_visit_to_a_built_area_queues_nothing_and_a_second_visitor_causes_no_build(galaxy, prefetch):
    assert prefetch.schedule(["5/16/16/16"]) == 20
    assert prefetch.schedule(["5/16/16/16"]) == 0, "the same area is claimed already"
    _run_queued(prefetch)
    assert prefetch.schedule(["5/16/16/16"]) == 14, "the rest of the area's 34 tiles, not the 20 built"
    _run_queued(prefetch)
    every = list(prefetch.module.candidates(["5/16/16/16"]))
    with warmup.build_app().test_request_context("/"):
        assert tilecache.fetch_tiles(prefetch.db, every)["cached"] == len(every) == 34
    assert prefetch.schedule(["5/16/16/16"]) == 0, "built: nothing left to queue"


def test_the_background_build_keeps_to_its_limits(galaxy, prefetch):
    prefetch.prefetch_per_visit, prefetch.prefetch_per_minute = 50, 30
    assert prefetch.schedule(["6/32/32/32"]) == 30, "the minute's cap trims the visit"
    assert prefetch.schedule(["6/10/10/10"]) == 0, "and the minute is spent"
    prefetch.prefetch_per_visit = 0
    assert prefetch.schedule(["6/50/50/50"]) == 0, "off"


def test_the_background_build_waits_while_the_cache_is_full_or_the_queue_busy(galaxy, prefetch):
    root = tilecache.cache_dir()
    prefetch.module.mark_full(root)
    assert prefetch.module.is_full(root) and prefetch.schedule(["7/64/64/64"]) == 0
    import os
    os.unlink(os.path.join(root, prefetch.module.FULL_FILE))
    prefetch.prefetch_max_queued = 1
    assert prefetch.schedule(["7/64/64/64"]) == 20      # one job per 16 tiles: two wait now
    assert prefetch.schedule(["7/90/90/90"]) == 0, "bulk work is waiting"


def test_an_answered_request_for_tiles_queues_the_tiles_around_it(galaxy, prefetch, monkeypatch):
    asked = []
    monkeypatch.setattr(prefetch.module, "schedule", lambda db, keys, at=None: asked.append((db, list(keys), at)) or 0)
    response = _web_worker(galaxy.config).get("/galaxy/tiles?tiles=5/16/16/16,5/17/16/16&at=10,20,30,1,2,3")
    assert response.status_code == 200
    response.close()
    assert asked == [(galaxy.config.database, ["5/16/16/16", "5/17/16/16"], (10.0, 20.0, 30.0, 1.0, 2.0, 3.0))]


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
    with pytest.raises(ValueError, match="no such cache family 'nonsense'"):
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
        return {"family": family, "warmed": 3, "left_out": 0, "seconds": 0.1,
                "line": f"{family}: 3 pages warmed (0 answers built, 3 already cached, 0 left out by the cap) in 0.1 s"}

    monkeypatch.setattr(warm_cli, "_site_database", lambda: "x")
    monkeypatch.setattr(warmup, "warm_family", warm_family)
    # The tile command is many units on a real galaxy; here, like the other families, it is one unit that is faked.
    monkeypatch.setattr(warmup, "expand_warm_galaxy_tiles", lambda params: commands.Expansion(
        [({"database": params["database"], "cap": params.get("cap")}, 1.0)], {}))
    monkeypatch.setattr(warmup, "warm_galaxy_tiles_unit", lambda args, report: warmup._warm_unit("galaxy-tiles", args))
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
    assert "the cache warm-up stopped: unit 0 of stage 4 failed: WarmupFailed: 2 of 9 pages did not answer" in err
    assert [family for family, _cap, _db in script.ran] == ["galaxy-map", "galaxy-tiles", "lists", "sectors", "systems"], \
        "the wrapper's list stops at the failing command"


def test_the_families_and_the_cap_are_passed_on(script):
    warm_cli.run_warm(_args("sectors", "lists", "--cap", "7"))
    assert [(family, cap) for family, cap, _db in script.ran] == [("lists", 7), ("sectors", 7)]


def test_an_unknown_family_is_an_error_before_anything_runs(script):
    with pytest.raises(SystemExit, match="no such cache family 'nonsense'"):
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
