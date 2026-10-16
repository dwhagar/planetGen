# planetgen/web/warmup.py

"""
Builds the site's caches ahead of the first visitor (MAP.134, PERF.104; design in
docs/design/cache-warm-up.md), so the first visit after a release, or after a cache was
cleared, finds the Galaxy Map, any sector, any system, the lists and the search already
there instead of waiting for the database.

The work is split into families (`FAMILIES`), each of which runs as a job of its own
(`run_family_job`, an RQ job on the interactive queue, so it is a node on the admin runs page
with its state and a line of result), or inline (`warm_family`). A family warms by visiting
the pages the way a visitor does, through the site's own Flask app, so what it caches is
exactly what the page asks for and the two can't drift apart; the answers go to the shared
page cache in Redis (`web/lib/pagecache.py`), the counts to `db/countcache.py`, the Galaxy
Map's tiles to the disk cache (`web/lib/tilecache.py`).

`opening_request` is the one place that says which tiles the Galaxy Map's first frame needs;
the `/galaxy` page asks for exactly them. `python -m planetgen.cli.warm_map` runs it all.
"""

import time

from planetgen.admin import activity_log
from planetgen.api.config import Config
from planetgen.api.limiter import PAGE_LIMITS_OFF
from planetgen.physics.units import ly_to_pc
from planetgen.queue import work as workQueue
from planetgen.tuning import DEFAULT_SECTOR_EDGE_LY
from planetgen.util.settings import get_settings
from planetgen.web.lib import apiclient
from planetgen.web.lib.tilecache import fetch_stage, fetch_tiles
from planetgen.web.maps.galaxymap3d import initial_tile_request, tile_bounds_pc, view_radius_bounds

PAGE_SIZE = 100
"""int: Rows asked for at a time when the warm-up lists the objects to visit."""

FAMILIES = {
    "galaxy-map": "the Galaxy Map's opening view, its shape, sectors and phenomena",
    "lists": "the first page of every list, and its counts",
    "sectors": "every sector page",
    "systems": "every system page",
    "phenomena": "every phenomenon page",
    "search": "the search page",
}
"""dict: Every family, in the order they run, and what it warms."""


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


_WARMERS = {"galaxy-map": _warm_galaxy_map, "lists": _warm_lists, "sectors": _warm_sectors,
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


def run_family_job(db, family, cap=None, config=None):
    """
    The RQ job for one family (`planetgen.queue.api_jobs.submit`): a node on the admin runs page that ends in the
    family's result line, and a line in the activity log. Raises `WarmupFailed` for a family that did not warm,
    which RQ records as a failed job and the node as failed.

    Args:
        db (str): The database the site shows.
        family (str): A key of `FAMILIES`.
        cap (int, optional): See `warm_family`.
        config (MySQLConfig, optional): The database and control database to use; the site's own when not given.

    Returns:
        dict: `warm_family`'s result.
    """
    title = f"Cache warm-up: {family} ({FAMILIES.get(family, 'unknown')})"
    with workQueue.job_node("warm", title, config, database=db) as node:
        try:
            result = warm_family(db, family, cap, config=config)
        except Exception as exc:
            workQueue.record_task(node, "warm-line", f"{family}: {exc}", state="failed", error=str(exc))
            activity_log.event("GEN", "warm.failed", db=db, family=family, error=str(exc)[:300])
            raise
        workQueue.record_task(node, "warm-line", result["line"], seconds=result["seconds"])
        activity_log.event("GEN", "warm.family", db=db, family=family, warmed=result["warmed"],
                           left_out=result["left_out"], seconds=round(result["seconds"], 1))
    return result
