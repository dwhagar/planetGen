# planetgen/web/warmup.py

"""
Builds the site's caches ahead of the first visitor (MAP.134, PERF.104; design in
docs/design/cache-warm-up.md), so the first visit after a release, or after a cache was
cleared, finds the Galaxy Map, any sector, any system, the lists and the search already
there instead of waiting for the database.

The work is split into families (`FAMILIES`), each of which is a command of its own on the work queue
(`warm-galaxy-map`, `warm-lists`, ... `queue/commands.py`; PERF.127) whose unit is the function here, or runs
inline (`warm_family`). `planetgen warm` hands the list of them to the wrapper. A family warms by visiting
the pages the way a visitor does, through the site's own Flask app, so what it caches is
exactly what the page asks for and the two can't drift apart; the answers go to the shared
page cache in Redis (`web/lib/pagecache.py`), the counts to `db/countcache.py`, the Galaxy
Map's tiles to the disk cache (`web/lib/tilecache.py`).

`opening_request` is the one place that says which tiles the Galaxy Map's first frame needs;
the `/galaxy` page asks for exactly them. `planetgen warm` runs it all.
"""

import math
import time

from planetgen.admin import activity_log
from planetgen.api.config import Config
from planetgen.api.limiter import PAGE_LIMITS_OFF
from planetgen.galaxy import viewport
from planetgen.physics.units import ly_to_pc
from planetgen.queue import commands
from planetgen.tuning import DEFAULT_SECTOR_EDGE_LY
from planetgen.util.settings import get_settings
from planetgen.web.lib import apiclient, tilecache
from planetgen.web.lib.tilecache import fetch_stage, fetch_tiles
from planetgen.web.maps.galaxymap3d import (
    CAMERA_START_PHI_DEG, CAMERA_START_THETA_DEG, camera_position_pc, initial_tile_request, tile_bounds_pc,
    view_radius_bounds)

PAGE_SIZE = 100
"""int: Rows asked for at a time when the warm-up lists the objects to visit."""

FAMILIES = {
    "galaxy-map": "the Galaxy Map's opening view, its shape, sectors and phenomena",
    "galaxy-tiles": "the Galaxy Map's tiles: the first levels in full, and the zoom run in to the centre",
    "lists": "the first page of every list, and its counts",
    "sectors": "every sector page",
    "systems": "every system page",
    "phenomena": "every phenomenon page",
    "search": "the search page",
}
"""dict: Every family, in the order they run, and what it warms."""

TILE_CHUNK = 128
"""int: Tiles one unit of `warm-galaxy-tiles` builds (the most a tile request holds, `tilecache.MAX_TILES_PER_REQUEST`)."""

MAX_WARM_TILES = 150000
"""int: Most tiles `warm-galaxy-tiles` will build in one run. `tile_cache.warm_levels` and `warm_zoom_levels` past
what this allows are refused, saying how many the levels would build, rather than filling the disk."""

ZOOM_STEP = viewport.DISTANCE_CUT_FADE_BAND
"""float: The zoom run's camera distance falls by this factor between two cuts; the cut of tiles changes by less than
the fetch of one step, so every tile the map would ask for on the way is in a cut of the run (MAP.152, the fade band)."""

TILE_SECONDS = ((6, 0.5), (7, 0.35), (8, 0.15))
"""tuple: `(level, seconds)`: the seconds one tile up to that level takes to build, measured on the scaled probe galaxy
(MAP.249), for the stage's bar and ETA; finer tiles take `TILE_SECONDS_FINE`."""

TILE_SECONDS_FINE = 0.08


class WarmupFailed(Exception):
    """A family could not warm what it covers; the message says what."""


def opening_request(galaxy_shape):
    """
    The tile keys the Galaxy Map's first frame needs.

    Args:
        galaxy_shape (dict or None): `apiclient.get_galaxy_shape`.

    Returns:
        list: Tile keys, nearest first.
    """
    edge_pc = galaxy_shape["edge_pc"] if galaxy_shape else ly_to_pc(DEFAULT_SECTOR_EDGE_LY)
    _min_radius, max_radius = view_radius_bounds(edge_pc, galaxy_shape)
    return initial_tile_request(max_radius, (0.0, 0.0, 0.0), tile_bounds_pc(edge_pc, galaxy_shape))


def warm_opening_view(db):
    """
    Fetches (and so caches) the opening view's tiles and the galaxy
    stage. Must run inside a Flask request context, so `apiclient` reaches
    the API in process.

    Returns:
        dict: `tiles` (how many), `cached` (how many of them were already
            on disk) and `seconds`.
    """
    started = time.monotonic()
    keys = opening_request(apiclient.get_galaxy_shape(db))
    result = fetch_tiles(db, keys)
    fetch_stage(db)
    return {"tiles": len(result["tiles"]), "cached": result["cached"], "seconds": time.monotonic() - started}


def tile_seconds(level):
    """About how long one tile of `level` takes to build, seconds (the weight of a unit in the stage's bar)."""
    return next((seconds for last, seconds in TILE_SECONDS if level <= last), TILE_SECONDS_FINE)


def _level_ranges(bounds_pc, level):
    """The inclusive index range on each axis of the level's tiles that meet `bounds_pc`."""
    root = viewport.TILE_ROOT_EDGE_PC
    edge = root / 2 ** level
    ranges = []
    for axis in range(3):
        low = max(0, math.floor((bounds_pc[0][axis] + root / 2) / edge))
        high = min(2 ** level - 1, math.floor((bounds_pc[1][axis] + root / 2) / edge))
        ranges.append(range(low, high + 1))
    return ranges


def octree_tiles(bounds_pc, levels):
    """
    Every tile of the first `levels` levels of the tile octree that meets the box `bounds_pc` (the cut's own,
    `tile_bounds_pc`): level 0, the root, down to level `levels` (so `levels` 8 is nine levels), coarsest level
    first and nearest the origin first within a level, so a run that is stopped has built the tiles the most views use.

    Raises:
        WarmupFailed: The levels would build more than `MAX_WARM_TILES` tiles.
    """
    levels = min(max(int(levels), 0), viewport.TILE_MAX_LEVEL)
    count = 0
    for level in range(levels + 1):
        xs, ys, zs = _level_ranges(bounds_pc, level)
        count += len(xs) * len(ys) * len(zs)
    if count > MAX_WARM_TILES:
        raise WarmupFailed(f"tile_cache.warm_levels {levels} would build {count} tiles, more than the {MAX_WARM_TILES} "
                           f"one run builds; lower it")
    keys = []
    for level in range(levels + 1):
        xs, ys, zs = _level_ranges(bounds_pc, level)
        keys += sorted((viewport.tile_key(level, ix, iy, iz) for ix in xs for iy in ys for iz in zs),
                       key=lambda key: (viewport._box_distance_sq(key, (0.0, 0.0, 0.0)), key))
    return keys


def zoom_run_tiles(bounds_pc, start_radius_pc, finest_level):
    """
    The tiles the map fetches on a zoom from the opening view in to the galaxy's centre, along the opening camera's
    line (the page's starting angles round the origin, the central black hole): the union of
    `viewport.fetch_cut_tiles` for the camera at `start_radius_pc`, then each `ZOOM_STEP` closer, down to tiles of
    `finest_level` (a tile of that level is in the cut from about its own edge in), keeping the tiles of that level
    and coarser, nearest first.

    Returns:
        list[str]: Tile keys, each once.
    """
    finest_level = min(max(int(finest_level), 0), viewport.TILE_MAX_LEVEL)
    direction = camera_position_pc(1.0, (0.0, 0.0, 0.0))
    target = (0.0, 0.0, 0.0)
    nearest = viewport.tile_edge_pc(finest_level)
    seen = {}
    distance = start_radius_pc
    while distance >= nearest:
        camera = tuple(distance * axis for axis in direction)
        for key in viewport.fetch_cut_tiles(camera, target, bounds_pc):
            if viewport.parse_tile_key(key)[0] <= finest_level:
                seen[key] = None
        distance /= ZOOM_STEP
    return sorted(seen, key=lambda key: (viewport.parse_tile_key(key)[0], viewport._box_distance_sq(key, target), key))


def tiles_to_warm(galaxy_shape, levels, zoom_levels):
    """
    The tiles `warm-galaxy-tiles` builds (PERF.130): `octree_tiles` of the first `levels` levels, then the zoom run's
    tiles finer than that, down to level `zoom_levels` (nothing when it is not finer than `levels`).

    Args:
        galaxy_shape (dict or None): `apiclient.get_galaxy_shape`.
        levels (int): `tile_cache.warm_levels`.
        zoom_levels (int): `tile_cache.warm_zoom_levels`.

    Returns:
        list[str]: Tile keys, each once, in the order they are built.
    """
    edge_pc = galaxy_shape["edge_pc"] if galaxy_shape else ly_to_pc(DEFAULT_SECTOR_EDGE_LY)
    bounds = tile_bounds_pc(edge_pc, galaxy_shape)
    keys = octree_tiles(bounds, levels)
    if zoom_levels > levels:
        _min_radius, start_radius = view_radius_bounds(edge_pc, galaxy_shape)
        taken = set(keys)
        run = [key for key in zoom_run_tiles(bounds, start_radius, zoom_levels) if key not in taken]
        if len(keys) + len(run) > MAX_WARM_TILES:
            raise WarmupFailed(f"the zoom run to level {zoom_levels} would build {len(keys) + len(run)} tiles in all, "
                               f"more than the {MAX_WARM_TILES} one run builds; lower tile_cache.warm_zoom_levels")
        keys += run
    return keys


def build_app(config=None):
    """The site's Flask app as the warm-up uses it: open to its own requests during maintenance (the sector
    switch keeps the site closed until the maps are rebuilt) and not rate limited.

    Args:
        config (MySQLConfig, optional): The database to show; the site's own when not given.
    """
    from planetgen.web.app import create_app

    class WarmupConfig(Config):
        WARMUP = True
        RATELIMIT_ENABLED = False
        RATELIMIT_PAGES = PAGE_LIMITS_OFF

    if config is not None:
        WarmupConfig.MYSQL_CONFIG = WarmupConfig.WRITE_MYSQL_CONFIG = WarmupConfig.CONTROL_MYSQL_CONFIG = config

    app = create_app(WarmupConfig)
    app.testing = False
    return app


class _Visitor:
    """Visits pages through the app the way a browser does, counting what answered from the cache."""

    def __init__(self, app, db):
        self.app, self.db = app, db
        self.client = app.test_client()
        self.pages = 0
        self.failed = []
        self._hits, self._misses = apiclient.cache_stats["hits"], apiclient.cache_stats["misses"]

    def visit(self, path):
        response = self.client.get(path)
        self.pages += 1
        if response.status_code != 200:
            self.failed.append(f"{path} answered {response.status_code}")
        return response

    @property
    def cached(self):
        return apiclient.cache_stats["hits"] - self._hits

    @property
    def built(self):
        return apiclient.cache_stats["misses"] - self._misses


def _cap(cap):
    return get_settings().page_cache.warm_max_objects if cap is None else int(cap)


def _objects(visitor, fetch, cap):
    """Yields the rows of a paginated listing, at most `cap`, and sets `visitor.total` to how many exist."""
    seen = 0
    visitor.total = 0
    while seen < cap:
        envelope = fetch(limit=min(PAGE_SIZE, cap - seen), offset=seen)
        visitor.total = envelope["total"]
        if not envelope["items"]:
            return
        for row in envelope["items"]:
            yield row
        seen += len(envelope["items"])
        if seen >= envelope["total"]:
            return


def _warm_galaxy_map(visitor, cap):
    opening = warm_opening_view(visitor.db)
    visitor.visit("/galaxy")
    return {"warmed": opening["tiles"], "noun": "tiles", "already": opening["cached"], "left_out": 0}


def _warm_galaxy_tiles(visitor, cap):
    """Warms every tile `warm-galaxy-tiles` builds, in this process (`warm_family`, the command's inline twin)."""
    settings = get_settings().tile_cache
    keys = tiles_to_warm(apiclient.get_galaxy_shape(visitor.db), settings.warm_levels, settings.warm_zoom_levels)
    warmed = already = 0
    for start in range(0, len(keys), TILE_CHUNK):
        result = fetch_tiles(visitor.db, keys[start:start + TILE_CHUNK])
        warmed += len(result["tiles"])
        already += result["cached"]
    return {"warmed": warmed, "noun": "tiles", "already": already, "left_out": 0}


def _warm_lists(visitor, cap):
    paths = ["/", "/sectors", "/systems", "/systems?uncharted=1", "/phenomena", "/classes"]
    status = apiclient.get_population_status(visitor.db)   # the pages are not there before the population pass
    paths += [path for kind, path in (("polities", "/polities"), ("species", "/species")) if status.get(kind)]
    for path in paths:
        visitor.visit(path)
    return {"warmed": visitor.pages, "noun": "pages", "left_out": 0}


def _warm_objects(visitor, cap, fetch, pages, noun):
    warmed = 0
    for row in _objects(visitor, fetch, _cap(cap)):
        for page in pages(row):
            visitor.visit(page)
        warmed += 1
    return {"warmed": warmed, "noun": noun, "left_out": max(0, visitor.total - warmed)}


def _warm_sectors(visitor, cap):
    return _warm_objects(
        visitor, cap, lambda **page: apiclient.get_sectors(visitor.db, **page),
        lambda row: (f"/sector/{row['id']}", f"/sector/{row['id']}/scene"), "sectors")


def _warm_systems(visitor, cap):
    return _warm_objects(
        visitor, cap, lambda **page: apiclient.get_systems(visitor.db, **page),
        lambda row: (f"/system/{row['id']}", f"/system/{row['id']}/scene"), "systems")


def _warm_phenomena(visitor, cap):
    return _warm_objects(
        visitor, cap, lambda **page: apiclient.get_phenomena(visitor.db, **page),
        lambda row: () if row.get("scattered") else (f"/phenomenon/{row['type']}/{row['id']}",), "phenomena")


def _warm_search(visitor, cap):
    visitor.visit("/search")
    return {"warmed": visitor.pages, "noun": "pages", "left_out": 0}


_WARMERS = {"galaxy-map": _warm_galaxy_map, "galaxy-tiles": _warm_galaxy_tiles, "lists": _warm_lists, "sectors": _warm_sectors,
            "systems": _warm_systems, "phenomena": _warm_phenomena, "search": _warm_search}


def warm_family(db, family, cap=None, app=None, config=None):
    """
    Warms one family inline.

    Args:
        db (str): The database the site shows.
        family (str): A key of `FAMILIES`.
        cap (int, optional): Most objects of an object family (default `page_cache.warm_max_objects`).
        app (flask.Flask, optional): The warm-up app (`build_app`).
        config (MySQLConfig, optional): The database to build the app on when `app` is not given.

    Returns:
        dict: `family`, `warmed`, `noun`, `cached` (answers that were already cached), `built`, `left_out` (objects
            past the cap), `seconds` and `line`, the one line that says it all.

    Raises:
        WarmupFailed: A page did not answer 200, or the family is unknown.
    """
    if family not in _WARMERS:
        raise WarmupFailed(f"no such cache family {family!r}; one of {', '.join(FAMILIES)}")
    started = time.monotonic()
    app = app or build_app(config)
    visitor = _Visitor(app, db)
    with app.test_request_context("/"):
        result = _WARMERS[family](visitor, cap)
    seconds = time.monotonic() - started
    result.update(family=family, cached=visitor.cached, built=visitor.built, seconds=seconds)
    result["line"] = (f"{family}: {result['warmed']} {result['noun']} warmed ({visitor.built} answers built, "
                      f"{visitor.cached} already cached, {result['left_out']} left out by the cap) in {seconds:.1f} s")
    if visitor.failed:
        shown = "; ".join(visitor.failed[:5]) + (f"; and {len(visitor.failed) - 5} more" if len(visitor.failed) > 5 else "")
        raise WarmupFailed(f"{family}: {len(visitor.failed)} of {visitor.pages} pages did not answer: {shown}")
    return result


def _warm_unit(family, args):
    """One family's unit of a warm run: warms it and leaves a line in the activity log, the result's line of it."""
    try:
        result = warm_family(args["database"], family, args.get("cap"))
    except Exception as exc:
        activity_log.event("GEN", "warm.failed", db=args["database"], family=family, error=str(exc)[:300])
        raise
    activity_log.event("GEN", "warm.family", db=args["database"], family=family, warmed=result["warmed"],
                       left_out=result["left_out"], seconds=round(result["seconds"], 1))
    return result


def _expand_family(params):
    """The one unit of a family's stage: the site's database and the cap, as the step gave them."""
    return commands.Expansion([({"database": params["database"], "cap": params.get("cap")}, 1.0)], {})


def warm_galaxy_map_unit(args, report):
    """The unit of the `warm-galaxy-map` command (PERF.127): the Galaxy Map's opening view, shape, sectors and phenomena."""
    return _warm_unit("galaxy-map", args)


def warm_lists_unit(args, report):
    """The unit of the `warm-lists` command: the first page of every list, and its counts."""
    return _warm_unit("lists", args)


def warm_sectors_unit(args, report):
    """The unit of the `warm-sectors` command: every sector page, up to `args["cap"]`."""
    return _warm_unit("sectors", args)


def warm_systems_unit(args, report):
    """The unit of the `warm-systems` command: every system page, up to `args["cap"]`."""
    return _warm_unit("systems", args)


def warm_phenomena_unit(args, report):
    """The unit of the `warm-phenomena` command: every phenomenon page, up to `args["cap"]`."""
    return _warm_unit("phenomena", args)


def warm_search_unit(args, report):
    """The unit of the `warm-search` command: the search page."""
    return _warm_unit("search", args)


expand_warm_galaxy_map = expand_warm_lists = expand_warm_sectors = _expand_family
expand_warm_systems = expand_warm_phenomena = expand_warm_search = _expand_family


def warm_galaxy_tiles_unit(args, report):
    """
    The unit of the `warm-galaxy-tiles` command (PERF.130): builds one chunk of the tiles into the disk cache, the
    same call the map's own request makes (`tilecache.fetch_tiles`), so what it stores is what a visitor is served.
    A tile that is already stored and current is not built again.

    Args:
        args (dict): `database` (the site's) and `keys` (the chunk, at most `TILE_CHUNK` tile keys).
        report (callable): Called with `(done, estimate)` when the chunk is built.
    """
    started = time.monotonic()
    count = len(args["keys"])
    if tilecache.cache_is_full():
        # The cache is at its budget (`cache.max_mb`): nothing more is built, the tiles built ahead stay (PERF.133).
        report(count, count)
        return {"family": "galaxy-tiles", "warmed": 0, "noun": "tiles", "cached": 0, "built": 0, "left_out": count,
                "seconds": 0.0,
                "line": f"galaxy-tiles: {count} tiles left out, " + (
                    "the disk cache is off (cache.max_mb is 0)" if tilecache.max_cache_bytes() <= 0 else
                    f"the tile cache is full ({tilecache.max_cache_bytes() / 1048576:.0f} MB, cache.max_mb)")}
    app = build_app()
    with app.test_request_context("/"):
        result = fetch_tiles(args["database"], args["keys"], ahead=True)
    report(count, count)
    seconds = time.monotonic() - started
    built = count - result["cached"]
    return {"family": "galaxy-tiles", "warmed": count, "noun": "tiles", "cached": result["cached"], "built": built,
            "left_out": 0, "seconds": seconds,
            "line": f"galaxy-tiles: {count} tiles warmed ({built} built, {result['cached']} already stored) in "
                    f"{seconds:.1f} s"}


def expand_warm_galaxy_tiles(params):
    """The units of the `warm-galaxy-tiles` stage: `tiles_to_warm` for the galaxy's shape and the named settings
    (`tile_cache.warm_levels`, `warm_zoom_levels`), in chunks of `TILE_CHUNK`, each weighted by about how long its
    tiles take to build."""
    database = params["database"]
    settings = get_settings().tile_cache
    app = build_app()
    with app.test_request_context("/"):
        shape = apiclient.get_galaxy_shape(database)
    keys = tiles_to_warm(shape, settings.warm_levels, settings.warm_zoom_levels)
    chunks = [keys[start:start + TILE_CHUNK] for start in range(0, len(keys), TILE_CHUNK)]
    units = [({"database": database, "keys": chunk},
              sum(tile_seconds(viewport.parse_tile_key(key)[0]) for key in chunk)) for chunk in chunks]
    return commands.Expansion(units, {"tiles": len(keys), "levels": settings.warm_levels,
                                      "zoom_levels": settings.warm_zoom_levels})


def finish_warm_galaxy_tiles(params, plan, results, seconds):
    """Leaves one line in the activity log for the whole `warm-galaxy-tiles` stage: tiles, how many were built and
    how many were already stored."""
    activity_log.event("GEN", "warm.family", db=params["database"], family="galaxy-tiles",
                       warmed=sum(result["warmed"] for result in results),
                       built=sum(result.get("built", 0) for result in results), left_out=0, seconds=round(seconds, 1))
    return {}


def end_map_rebuild_unit(args, report):
    """
    The unit of the `end-map-rebuild-maintenance` command: opens the site if the sector-grid switch closed it
    (`maintenance.SECTOR_SWITCH`, UX.112), the Galaxy Map being back on the new grid. An admin's own maintenance, and one
    with another cause, stays.

    Raises:
        planetgen.admin.maintenance.MaintenanceWriteError: The site could not be opened; the run fails with it.
    """
    from planetgen.admin import maintenance
    if maintenance.turn_off_for(maintenance.SECTOR_SWITCH, by="map rebuild"):
        return {"reopened": True, "line": "Maintenance mode is off: the Galaxy Map is rebuilt."}
    return {"reopened": False, "line": "Maintenance mode was not the map rebuild's, so it is left as it is."}


def expand_end_map_rebuild(params):
    """The one unit of the stage that opens the site."""
    return commands.Expansion([({}, 1.0)], {})
