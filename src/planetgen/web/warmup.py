# planetgen/web/warmup.py

"""
Builds the site's caches ahead of the first visitor (MAP.134, PERF.104; design in
docs/design/operations/cache-warm-up.md), so the first visit after a release, or after a cache was
cleared, finds the Galaxy Map, any sector, any system, the lists and the search already
there instead of waiting for the database.

The work is split into families (`FAMILIES`), each of which is a command of its own on the work queue
(`warm-galaxy-map`, `warm-lists`, ... `queue/commands.py`; PERF.127) whose unit is the function here, or runs
inline (`warm_family`). `planetgen warm` hands the list of them to the wrapper. A family warms by visiting
the pages the way a visitor does, through the site's own Flask app, so what it caches is
exactly what the page asks for and the two can't drift apart; the answers go to the shared
page cache in Redis (`web/lib/pagecache.py`), the counts to `db/countcache.py`, the Galaxy
Map's nodes to the disk cache (`web/lib/tilecache.py`).

`opening_request` is the one place that says which nodes the Galaxy Map's first frame needs;
the `/galaxy` page asks for exactly them. `planetgen warm` runs it all.
"""

import time

from planetgen.admin import activity_log
from planetgen.api.config import Config
from planetgen.api.limiter import PAGE_LIMITS_OFF
from planetgen.db import nodes, objectindex, store
from planetgen.galaxy.nodecut import node_cut
from planetgen.physics.units import ly_to_pc
from planetgen.queue import commands
from planetgen.queue import units as units_module
from planetgen.tuning import DEFAULT_SECTOR_EDGE_LY
from planetgen.util.settings import get_settings
from planetgen.web.lib import apiclient, tilecache
from planetgen.web.lib.tilecache import fetch_nodes, fetch_stage
from planetgen.web.maps.galaxymap3d import (
    camera_position_pc, initial_node_request, node_bounds_pc,
    view_radius_bounds)

PAGE_SIZE = 100
"""int: Rows asked for at a time when the warm-up lists the objects to visit."""

PHASES = (
    ("opening", "the galaxy's cluster field, the object index, the opening views and the first page of each list",
     ("cluster-field", "object-index", "galaxy-map", "lists")),
    ("built-levels", "the Galaxy Map's sector-block nodes down to the built level", ("galaxy-nodes",)),
    ("deep-nodes", "the deep nodes along the zoom run in to the centre", ("galaxy-deep-nodes",)),
    ("objects", "the sector, system and phenomenon pages", ("sectors", "systems", "phenomena")),
    ("search", "the search page and its answers", ("search",)),
)
"""tuple: The warm-up's phases in the order they run (PERF.141): `(name, what it warms, the families in it)`. A phase
is the unit the run reports and the admin page shows; each family is a command of its own."""

_FAMILY_TEXT = {
    "cluster-field": "the galaxy's clusters and satellites as one field: each age's mean system mass worked out in small "
                     "parts, then the field kept as a file in the state folder (`galaxy/field_job.py`)",
    "object-index": "the object index: one row per object of every kind the site lists",
    "galaxy-map": "the Galaxy Map's opening view, its shape, sectors and phenomena",
    "lists": "the first page of every list, and its counts",
    "galaxy-nodes": "the Galaxy Map's sector-block nodes: every populated block from the whole galaxy down to level "
                    "`NODE_WARM_FLOOR`",
    "galaxy-deep-nodes": "the Galaxy Map's deep nodes: the zoom run in to the centre, down to `tile_cache.warm_zoom_floor`",
    "sectors": "every sector page",
    "systems": "every system page",
    "phenomena": "every phenomenon page",
    "search": "the search page",
}
FAMILIES = {family: _FAMILY_TEXT[family] for _name, _text, families in PHASES for family in families}
"""dict: Every family, in the order they run (the phases', see `PHASES`), and what it warms."""


class PhaseTally:
    """
    Adds up the results of the families of a warm-up run and says, as each phase ends, what it did (PERF.141): how
    much it built, how much was already stored or cached, how much the budget or a cap left out, and, when it built
    nothing at all, that it was `unchanged`.

    Args:
        families (list[str]): The families of the run, in the order they run.
    """

    def __init__(self, families):
        self.families = [family for family in FAMILIES if family in families]
        self.totals = {family: {"built": 0, "cached": 0, "warmed": 0, "left_out": 0, "seconds": 0.0}
                       for family in self.families}
        self._done = set()
        self._reported = set()
        self.ended = []
        """list[str]: The phases `start` or `finish` reported last, in order."""

    def add(self, family, result):
        """Adds one unit's result (the dict a family's unit returns) to `family`'s totals."""
        total = self.totals[family]
        for key in ("built", "cached", "warmed", "left_out"):
            total[key] += int(result.get(key) or 0)
        total["seconds"] += float(result.get("seconds") or 0.0)

    def result(self, name):
        """What phase `name` did, as the admin page keeps it: `{"state", "built", "cached", "left_out", "seconds",
        "line"}` with `state` one of `built`, `unchanged`, `stopped` and `nothing` (`None` when none of its families ran)."""
        families = [family for family in self.totals if phase_of(family) == name]
        if not families:
            return None
        built, cached, warmed, left_out, seconds = (
            sum(self.totals[family][key] for family in families) for key in ("built", "cached", "warmed", "left_out", "seconds"))
        if not warmed and not left_out:
            state = "nothing"
        elif not built and not left_out:
            state = "unchanged"
        else:
            state = "stopped" if left_out else "built"
        return {"state": state, "built": built, "cached": cached, "left_out": left_out, "seconds": round(seconds, 1),
                "line": self.line(name)}

    def line(self, name):
        """The one line for phase `name`: `Phase 2 of 5, built-levels: ...` (`None` when none of its families ran)."""
        number = next(index for index, (phase, _text, _families) in enumerate(PHASES, 1) if phase == name)
        families = [family for family in self.totals if phase_of(family) == name]
        if not families:
            return None
        built, cached, warmed, left_out, seconds = (
            sum(self.totals[family][key] for family in families) for key in ("built", "cached", "warmed", "left_out", "seconds"))
        head = f"Phase {number} of {len(PHASES)}, {name}"
        if not warmed and not left_out:
            return f"{head}: nothing to warm."
        if not built and not left_out:
            return f"{head}: unchanged, nothing built ({cached:,} already stored) in {seconds:.1f} s."
        return (f"{head}: {built:,} built, {cached:,} already stored, {left_out:,} left out in {seconds:.1f} s"
                f"{'; stopped at the cache budget' if left_out and name in ('built-levels', 'deep-nodes') else ''}.")

    def _ready(self):
        lines, ended = [], []
        for name, _text, _families in PHASES:
            ours = [family for family in self.totals if phase_of(family) == name]
            if ours and name not in self._reported and all(family in self._done for family in ours):
                self._reported.add(name)
                lines.append(self.line(name))
                ended.append(name)
        self.ended = ended
        return lines

    def start(self, step_key):
        """A step is starting: every family before it is done. Returns the lines of the phases that ended."""
        if step_key in self.totals:
            self._done.update(self.families[:self.families.index(step_key)])
        else:
            self._done.update(self.families)
        return self._ready()

    def finish(self):
        """The run ended: returns the lines of the phases not reported yet."""
        self._done.update(self.families)
        return self._ready()


def phase_of(family):
    """The name of the phase `family` belongs to."""
    return next(name for name, _text, families in PHASES if family in families)


def expand_phases(names):
    """
    `names`, family names and phase names mixed, as the families they stand for (a phase is all of its families),
    each once, in the order of `FAMILIES`.

    Raises:
        ValueError: A name that is neither a family nor a phase.
    """
    wanted = set()
    for name in names:
        phase = next((families for phase_name, _text, families in PHASES if phase_name == name), None)
        if phase is not None:
            wanted.update(phase)
        elif name in FAMILIES:
            wanted.add(name)
        else:
            raise ValueError(f"no such cache family or phase {name!r}; families: {', '.join(FAMILIES)}; "
                             f"phases: {', '.join(name for name, _text, _families in PHASES)}")
    return [family for family in FAMILIES if family in wanted]

MAX_WARM_NODES = 100000
"""int: Most nodes the zoom run's `warm-galaxy-deep-nodes` will build in one run. A `tile_cache.warm_zoom_floor` that
reaches more is refused, saying how many nodes it would build, rather than filling the disk."""

ZOOM_STEP = 1.3
"""float: The zoom run's camera distance falls by this factor between two cuts; the cut of nodes changes by less than
the fetch of one step, so every node the map would ask for on the way is in a cut of the run."""

NODE_CHUNK = 64
"""int: Nodes one unit of `warm-galaxy-nodes` and `warm-galaxy-deep-nodes` builds (`tilecache.MAX_NODES_PER_REQUEST` is
the most a request holds)."""

NODE_WARM_FLOOR = 243
"""int: The finest block level `warm-galaxy-nodes` builds ahead (MAP.260): the populated blocks of that level and every
bigger one, the nodes the map opens on its way in; finer blocks are built as a visit asks, then prefetched around it."""

NODE_SECONDS = 0.3
"""float: About how long one node takes to build, seconds (the weight of a unit in the stage's bar)."""

class WarmupFailed(Exception):
    """A family could not warm what it covers; the message says what."""


def opening_request(galaxy_shape):
    """
    The node keys the Galaxy Map's first frame needs.

    Args:
        galaxy_shape (dict or None): `apiclient.get_galaxy_shape`.

    Returns:
        list: Node keys, the coarsest first (`galaxymap3d.initial_node_request`).
    """
    edge_pc = galaxy_shape["edge_pc"] if galaxy_shape else ly_to_pc(DEFAULT_SECTOR_EDGE_LY)
    _min_radius, max_radius = view_radius_bounds(edge_pc, galaxy_shape)
    return initial_node_request(max_radius, (0.0, 0.0, 0.0), edge_pc, node_bounds_pc(edge_pc, galaxy_shape))


def warm_opening_view(db):
    """
    Fetches (and so caches) the opening view's nodes and the galaxy
    stage. Must run inside a Flask request context, so `apiclient` reaches
    the API in process.

    Returns:
        dict: `nodes` (how many), `cached` (how many of them were already
            on disk) and `seconds`.
    """
    started = time.monotonic()
    keys = opening_request(apiclient.get_galaxy_shape(db))
    result = fetch_nodes(db, keys)
    fetch_stage(db)
    return {"nodes": len(result["nodes"]), "cached": result["cached"], "seconds": time.monotonic() - started}


def zoom_run_nodes(edge_pc, bounds_pc, start_radius_pc, finest_block):
    """
    The nodes the map fetches on a zoom from the opening view in to the galaxy's centre, along the opening camera's
    line (the page's starting angles round the origin, the central black hole): the union of `nodecut.node_cut` for the
    camera at `start_radius_pc`, then each `ZOOM_STEP` closer, down to blocks of `finest_block` shells (a block of that
    size is in the cut from about a few of its own sizes in), keeping the nodes finer than `NODE_WARM_FLOOR` and not
    finer than `finest_block`, coarsest first.

    Returns:
        list[str]: Node keys, each once.
    """
    direction = camera_position_pc(1.0, (0.0, 0.0, 0.0))
    target = (0.0, 0.0, 0.0)
    nearest = finest_block * edge_pc
    seen = {}
    distance = start_radius_pc
    while distance >= nearest:
        camera = tuple(distance * axis for axis in direction)
        for key in node_cut(camera, target, edge_pc, bounds_pc=bounds_pc).keys:
            if finest_block <= int(key.split(".", 1)[0]) < NODE_WARM_FLOOR:
                seen[key] = None
        distance /= ZOOM_STEP
    return sorted(seen, key=lambda key: (-int(key.split(".", 1)[0]), key))


def deep_nodes_to_warm(galaxy_shape, finest_block):
    """
    The nodes `warm-galaxy-deep-nodes` builds: `zoom_run_nodes` for the galaxy's shape, down to blocks of
    `finest_block` shells.

    Args:
        galaxy_shape (dict or None): `apiclient.get_galaxy_shape`.
        finest_block (int): `tile_cache.warm_zoom_floor`; at `NODE_WARM_FLOOR` or more there is nothing to build.

    Returns:
        list[str]: Node keys, each once, in the order they are built.

    Raises:
        WarmupFailed: The run would build more than `MAX_WARM_NODES` nodes.
    """
    if finest_block >= NODE_WARM_FLOOR:
        return []
    edge_pc = galaxy_shape["edge_pc"] if galaxy_shape else ly_to_pc(DEFAULT_SECTOR_EDGE_LY)
    _min_radius, start_radius = view_radius_bounds(edge_pc, galaxy_shape)
    keys = zoom_run_nodes(edge_pc, node_bounds_pc(edge_pc, galaxy_shape), start_radius, finest_block)
    if len(keys) > MAX_WARM_NODES:
        raise WarmupFailed(f"the zoom run to blocks of {finest_block} shells would build {len(keys)} nodes, more than "
                           f"the {MAX_WARM_NODES} one run builds; raise tile_cache.warm_zoom_floor")
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


def _objects(visitor, fetch, cap, start=0):
    """Yields the rows of a paginated listing from place `start`, at most `cap` of them, and sets `visitor.total` to how many
    exist."""
    seen = 0
    visitor.total = 0
    while seen < cap:
        envelope = fetch(limit=min(PAGE_SIZE, cap - seen), offset=start + seen)
        visitor.total = envelope["total"]
        if not envelope["items"]:
            return
        for row in envelope["items"]:
            yield row
        seen += len(envelope["items"])
        if start + seen >= envelope["total"]:
            return


def _warm_galaxy_map(visitor, cap):
    opening = warm_opening_view(visitor.db)
    visitor.visit("/galaxy")
    return {"warmed": opening["nodes"], "noun": "nodes", "already": opening["cached"], "left_out": 0}


def _populated_nodes(database):
    """The keys `warm-galaxy-nodes` builds, from the stored block rows (`nodes.populated_node_keys`)."""
    conn = store.get_connection(store.MySQLConfig(database=database), ensure_schema=False)
    try:
        return nodes.populated_node_keys(conn, NODE_WARM_FLOOR)
    finally:
        conn.close()


def _warm_cluster_field(visitor, cap):
    """Makes the cluster field in this process (`warm_family`, the command's inline twin), part after part."""
    from planetgen.galaxy import field_job
    built, cached = field_job.warm_inline(visitor.db)
    return {"warmed": built + cached, "noun": "parts", "already": cached, "left_out": 0}


def _warm_galaxy_nodes(visitor, cap):
    """Warms every node `warm-galaxy-nodes` builds, in this process (`warm_family`, the command's inline twin)."""
    keys = _populated_nodes(visitor.db)
    warmed = already = 0
    for start in range(0, len(keys), NODE_CHUNK):
        result = fetch_nodes(visitor.db, keys[start:start + NODE_CHUNK])
        warmed += len(result["nodes"])
        already += result["cached"]
    return {"warmed": warmed, "noun": "nodes", "already": already, "left_out": 0}


def _warm_lists(visitor, cap):
    from planetgen.web.object_lists import KIND_LABELS
    paths = ["/", "/sectors", "/systems", "/systems?uncharted=1", "/phenomena", "/classes", "/lists"]
    paths += [f"/lists/{kind}" for kind in KIND_LABELS]   # UX.141: a keyset page costs the same at any depth
    status = apiclient.get_population_status(visitor.db)   # the pages are not there before the population pass
    paths += [path for kind, path in (("polities", "/polities"), ("species", "/species")) if status.get(kind)]
    for path in paths:
        visitor.visit(path)
    return {"warmed": visitor.pages, "noun": "pages", "left_out": 0}


def _warm_objects(visitor, cap, fetch, pages, noun, window=None):
    """Warms the pages of a listing's objects, up to `cap`. With a `window` (`{"offset", "count", "left_out"}`, PERF.180) it warms
    that slice of the listing alone, one unit of the family's stage, and `left_out` is the share of the cap's leaving-out that
    unit reports."""
    warmed = 0
    start = 0 if window is None else window["offset"]
    for row in _objects(visitor, fetch, _cap(cap) if window is None else window["count"], start):
        for page in pages(row):
            visitor.visit(page)
        warmed += 1
    left_out = max(0, visitor.total - warmed) if window is None else window["left_out"]
    return {"warmed": warmed, "noun": noun, "left_out": left_out}


def _warm_sectors(visitor, cap, window=None):
    return _warm_objects(
        visitor, cap, lambda **page: apiclient.get_sectors(visitor.db, **page),
        lambda row: (f"/sector/{row['id']}", f"/sector/{row['id']}/scene"), "sectors", window)


def _warm_systems(visitor, cap, window=None):
    return _warm_objects(
        visitor, cap, lambda **page: apiclient.get_systems(visitor.db, **page),
        lambda row: (f"/system/{row['id']}", f"/system/{row['id']}/scene"), "systems", window)


def _warm_phenomena(visitor, cap, window=None):
    return _warm_objects(
        visitor, cap, lambda **page: apiclient.get_phenomena(visitor.db, **page),
        lambda row: () if row.get("scattered") else (f"/phenomenon/{row['type']}/{row['id']}",), "phenomena", window)


def _warm_search(visitor, cap):
    visitor.visit("/search")
    return {"warmed": visitor.pages, "noun": "pages", "left_out": 0}


def _warm_galaxy_deep_nodes(visitor, cap):
    """Warms every node `warm-galaxy-deep-nodes` builds, in this process."""
    keys = deep_nodes_to_warm(apiclient.get_galaxy_shape(visitor.db), get_settings().tile_cache.warm_zoom_floor)
    warmed = already = 0
    for start in range(0, len(keys), NODE_CHUNK):
        result = fetch_nodes(visitor.db, keys[start:start + NODE_CHUNK])
        warmed += len(result["nodes"])
        already += result["cached"]
    return {"warmed": warmed, "noun": "nodes", "already": already, "left_out": 0}


_WARMERS = {"cluster-field": _warm_cluster_field, "galaxy-map": _warm_galaxy_map, "galaxy-nodes": _warm_galaxy_nodes,
            "galaxy-deep-nodes": _warm_galaxy_deep_nodes, "lists": _warm_lists, "sectors": _warm_sectors,
            "systems": _warm_systems, "phenomena": _warm_phenomena, "search": _warm_search}


def warm_family(db, family, cap=None, app=None, config=None, window=None):
    """
    Warms one family inline.

    Args:
        db (str): The database the site shows.
        family (str): A key of `FAMILIES`.
        cap (int, optional): Most objects of an object family (default `page_cache.warm_max_objects`).
        app (flask.Flask, optional): The warm-up app (`build_app`).
        config (MySQLConfig, optional): The database to build the app on when `app` is not given.
        window (dict, optional): For an object family, the slice of its listing to warm (`{"offset", "count", "left_out"}`): one
            unit of the family's stage (PERF.180). Without it the whole listing up to the cap.

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
        result = _WARMERS[family](visitor, cap, **({} if window is None else {"window": window}))
    seconds = time.monotonic() - started
    result.update(family=family, cached=visitor.cached, built=visitor.built, seconds=seconds)
    result["line"] = (f"{family}: {result['warmed']} {result['noun']} warmed ({visitor.built} answers built, "
                      f"{visitor.cached} already cached, {result['left_out']} left out by the cap) in {seconds:.1f} s")
    if visitor.failed:
        shown = "; ".join(visitor.failed[:5]) + (f"; and {len(visitor.failed) - 5} more" if len(visitor.failed) > 5 else "")
        raise WarmupFailed(f"{family}: {len(visitor.failed)} of {visitor.pages} pages did not answer: {shown}")
    return result


def _warm_unit(family, args):
    """One family's unit of a warm run: warms it and leaves a line in the activity log, the result's line of it. A unit that is a
    slice of an object family (`args["offset"]`, PERF.180) leaves no line of its own: the stage's `finish` leaves one for all."""
    window = None if "offset" not in args else {key: args[key] for key in ("offset", "count", "left_out")}
    try:
        result = warm_family(args["database"], family, args.get("cap"), **({} if window is None else {"window": window}))
    except Exception as exc:
        activity_log.event("GEN", "warm.failed", db=args["database"], family=family, error=str(exc)[:300])
        raise
    if window is None:
        activity_log.event("GEN", "warm.family", db=args["database"], family=family, warmed=result["warmed"],
                           left_out=result["left_out"], seconds=round(result["seconds"], 1))
    return result


def _expand_family(params):
    """The one unit of a family's stage: the site's database, the cap and the index mode, as the step gave them."""
    unit = {"database": params["database"], "cap": params.get("cap")}
    if params.get("index_mode"):
        unit["index_mode"] = params["index_mode"]
    return commands.Expansion([(unit, 1.0)], {})


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


def warm_object_index_unit(args, report):
    """
    The unit of the `warm-object-index` command (PERF.140): refreshes the object index (builds the kinds that have no
    record of a run), or with `index_mode` `build` writes every object again. Reports like a page family: the rows
    are what was warmed, the rows it wrote are what it built.
    """
    started = time.monotonic()
    database = args["database"]
    try:
        conn = store.get_connection(store.MySQLConfig(database=database), ensure_schema=False)
        try:
            build = args.get("index_mode") == "build"
            found = (objectindex.build if build else objectindex.refresh)(conn, database=database)
        finally:
            conn.close()
    except Exception as exc:
        activity_log.event("GEN", "warm.failed", db=database, family="object-index", error=str(exc)[:300])
        raise
    rows = sum(kind["rows"] for kind in found.values())
    written = sum(kind["changed"] + kind["removed"] for kind in found.values())
    seconds = time.monotonic() - started
    activity_log.event("GEN", "warm.family", db=database, family="object-index", warmed=rows, left_out=0,
                       seconds=round(seconds, 1))
    return {"family": "object-index", "warmed": rows, "noun": "objects", "cached": max(rows - written, 0),
            "built": written, "left_out": 0, "seconds": seconds,
            "line": (f"object-index: {rows:,} objects indexed ({written:,} rows written, {max(rows - written, 0):,} "
                     f"already current, 0 left out by the cap) in {seconds:.1f} s")}


OBJECT_CHUNK = 250
"""int: How many objects of a sector, system or phenomenon family one unit warms (PERF.180): the page and its scene, two answers
each. The stage is cut into units of this many (a whole multiple of it past `units.MAX_UNITS` units, `units.coarsen`)."""

OBJECT_SECONDS = 0.3
"""float: About how long one object takes to warm (its page and scene built, none cached), the weight of a unit and, times
`OBJECT_CHUNK`, the longest unit the command declares (`queue.commands`)."""

_OBJECT_TOTALS = {"sectors": "get_sectors", "systems": "get_systems", "phenomena": "get_phenomena"}


def _expand_objects(params, family):
    """
    The units of an object family's stage (PERF.180): the listing up to the cap, in slices of `OBJECT_CHUNK` objects, each unit the
    `offset` and `count` of its slice. The first unit also carries the objects the cap leaves out. A listing with nothing in it is one
    unit that says so.
    """
    database = params["database"]
    cap = _cap(params.get("cap"))
    app = build_app()
    with app.test_request_context("/"):
        total = int(getattr(apiclient, _OBJECT_TOTALS[family])(database, limit=1)["total"])
    wanted = min(total, cap)
    size = units_module.coarsen(wanted, OBJECT_CHUNK)
    units = [({"database": database, "cap": params.get("cap"), "offset": start, "count": min(size, wanted - start),
               "left_out": total - wanted if start == 0 else 0}, OBJECT_SECONDS * min(size, wanted - start))
             for start in range(0, wanted, size)]
    if not units:
        units = [({"database": database, "cap": params.get("cap"), "offset": 0, "count": 0, "left_out": 0}, OBJECT_SECONDS)]
    return commands.Expansion(units, {family: total, "cap": cap})


def expand_warm_sectors(params):
    """The units of the `warm-sectors` stage: slices of the sector listing (`_expand_objects`)."""
    return _expand_objects(params, "sectors")


def expand_warm_systems(params):
    """The units of the `warm-systems` stage: slices of the system listing (`_expand_objects`)."""
    return _expand_objects(params, "systems")


def expand_warm_phenomena(params):
    """The units of the `warm-phenomena` stage: slices of the phenomenon listing (`_expand_objects`)."""
    return _expand_objects(params, "phenomena")


def finish_warm_sectors(params, plan, results, seconds):
    """Leaves one line in the activity log for the whole `warm-sectors` stage."""
    return _finish_nodes(params, results, seconds, "sectors")


def finish_warm_systems(params, plan, results, seconds):
    """Leaves one line in the activity log for the whole `warm-systems` stage."""
    return _finish_nodes(params, results, seconds, "systems")


def finish_warm_phenomena(params, plan, results, seconds):
    """Leaves one line in the activity log for the whole `warm-phenomena` stage."""
    return _finish_nodes(params, results, seconds, "phenomena")


expand_warm_object_index = _expand_family
expand_warm_galaxy_map = expand_warm_lists = _expand_family
expand_warm_search = _expand_family


def _fetch_in_requests(app, args, report):
    """
    Builds the unit's `args["keys"]` with `fetch_nodes`, `NODE_CHUNK` keys per call, the most one request holds. A unit
    holds more keys than that on a big galaxy, when the stage's cut is made coarser to stay within `units.MAX_UNITS`
    (`units.coarsen`). The cache is looked at before each call after the first: once it is full the rest of the unit is
    left out, as a unit that finds it full at the start leaves out all of its keys.

    Returns:
        tuple: `(cached, left_out)`: how many keys were already stored and how many were not built.
    """
    keys = args["keys"]
    cached = 0
    for start in range(0, len(keys), NODE_CHUNK):
        if start and tilecache.cache_is_full():
            report(start, len(keys))
            return cached, len(keys) - start
        with app.test_request_context("/"):
            cached += fetch_nodes(args["database"], keys[start:start + NODE_CHUNK], ahead=True)["cached"]
        report(min(start + NODE_CHUNK, len(keys)), len(keys))
    return cached, 0


def _warm_node_chunk(args, report, family):
    """
    The unit of the `warm-galaxy-nodes` and `warm-galaxy-deep-nodes` commands (MAP.260): builds one chunk of the nodes
    into the disk cache, the same call the map's own request makes (`tilecache.fetch_nodes`), so what it stores is what
    a visitor is served. A node that is already stored and current is not built again.

    Args:
        args (dict): `database` (the site's) and `keys` (the chunk: `NODE_CHUNK` node keys, or a whole multiple of that on a galaxy too big for
            the stage's unit cap, see `units.coarsen`).
        report (callable): Called with `(done, estimate)` when the chunk is built.
        family (str): The family the chunk belongs to, for its line.
    """
    started = time.monotonic()
    count = len(args["keys"])
    if not count:
        report(0, 0)
        return {"family": family, "warmed": 0, "noun": "nodes", "cached": 0, "built": 0, "left_out": 0,
                "seconds": 0.0, "line": f"{family}: no node to warm" + (
                    ", no block is stored yet" if family == "galaxy-nodes" else "")}
    if tilecache.cache_is_full():
        report(count, count)
        return {"family": family, "warmed": 0, "noun": "nodes", "cached": 0, "built": 0, "left_out": count,
                "seconds": 0.0,
                "line": f"{family}: {count} nodes left out, " + (
                    "the disk cache is off (cache.max_mb is 0)" if tilecache.max_cache_bytes() <= 0 else
                    f"the cache is full ({tilecache.max_cache_bytes() / 1048576:.0f} MB, cache.max_mb)")}
    app = build_app()
    cached, left_out = _fetch_in_requests(app, args, report)
    seconds = time.monotonic() - started
    warmed = count - left_out
    built = warmed - cached
    return {"family": family, "warmed": warmed, "noun": "nodes", "cached": cached, "built": built,
            "left_out": left_out, "seconds": seconds,
            "line": f"{family}: {warmed} nodes warmed ({built} built, {cached} already stored) in {seconds:.1f} s"
                    + (f", {left_out} left out as the cache filled" if left_out else "")}


def warm_galaxy_nodes_unit(args, report):
    """The unit of the `warm-galaxy-nodes` command (MAP.260): one chunk of the built levels' nodes."""
    return _warm_node_chunk(args, report, "galaxy-nodes")


def warm_galaxy_deep_nodes_unit(args, report):
    """The unit of the `warm-galaxy-deep-nodes` command: one chunk of the zoom run's deep nodes."""
    return _warm_node_chunk(args, report, "galaxy-deep-nodes")


def expand_warm_galaxy_nodes(params):
    """The units of the `warm-galaxy-nodes` stage: the populated blocks of `NODE_WARM_FLOOR` and above, in chunks of
    `NODE_CHUNK` (a whole multiple of it past `units.MAX_UNITS` chunks), coarsest first."""
    database = params["database"]
    keys = _populated_nodes(database)
    size = units_module.coarsen(len(keys), NODE_CHUNK)
    units = [({"database": database, "keys": keys[start:start + size]},
              NODE_SECONDS * len(keys[start:start + size])) for start in range(0, len(keys), size)]
    if not units:       # nothing is stored yet: one unit with no nodes, so the stage still runs and says so
        units = [({"database": database, "keys": []}, NODE_SECONDS)]
    return commands.Expansion(units, {"nodes": len(keys), "floor": NODE_WARM_FLOOR})


def finish_warm_galaxy_nodes(params, plan, results, seconds):
    """Leaves one line in the activity log for the whole `warm-galaxy-nodes` stage."""
    return _finish_nodes(params, results, seconds, "galaxy-nodes")


def expand_warm_galaxy_deep_nodes(params):
    """The units of the `warm-galaxy-deep-nodes` stage: `deep_nodes_to_warm` for the galaxy's shape and
    `tile_cache.warm_zoom_floor`, in chunks of `NODE_CHUNK` (a whole multiple of it past `units.MAX_UNITS` chunks), coarsest first."""
    database = params["database"]
    floor = get_settings().tile_cache.warm_zoom_floor
    app = build_app()
    with app.test_request_context("/"):
        shape = apiclient.get_galaxy_shape(database)
    keys = deep_nodes_to_warm(shape, floor)
    size = units_module.coarsen(len(keys), NODE_CHUNK)
    units = [({"database": database, "keys": keys[start:start + size]},
              NODE_SECONDS * len(keys[start:start + size])) for start in range(0, len(keys), size)]
    if not units:
        units = [({"database": database, "keys": []}, NODE_SECONDS)]
    return commands.Expansion(units, {"nodes": len(keys), "zoom_floor": floor})


def _finish_nodes(params, results, seconds, family):
    activity_log.event("GEN", "warm.family", db=params["database"], family=family,
                       warmed=sum(result["warmed"] for result in results),
                       built=sum(result.get("built", 0) for result in results),
                       left_out=sum(result.get("left_out", 0) for result in results), seconds=round(seconds, 1))
    return {}


def finish_warm_galaxy_deep_nodes(params, plan, results, seconds):
    """Leaves one line in the activity log for the whole `warm-galaxy-deep-nodes` stage."""
    return _finish_nodes(params, results, seconds, "galaxy-deep-nodes")


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
