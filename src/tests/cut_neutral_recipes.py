# tests/cut_neutral_recipes.py

"""
The recipes of the cut-neutral check (PERF.260, `tests/cut_neutral_support.py`): for each command whose stage is cut into units, how to put a small database in the state the stage works
on and which constant cuts it. They are kept by command kind (`RECIPES`); a command flagged "dynamic split is byte identical" with no recipe here fails the merge gate.
"""

from planetgen.queue import wrapper
from tests.cut_neutral_support import Recipe

LONG_AGO = "2000-01-01 00:00:00"


def _database(context):
    return wrapper.database_args(context.mysql_config)


# --- link-neighbors and settle-sector-paths: runs of sectors (PERF.113, PERF.115) -------------------------------

def _save_sectors(context):
    from tests import test_sector_paths_db as helper

    for wedge in (700, 701, 702, 703):
        context.monkeypatch.setattr(helper, "ADDRESS", (1500, 2, wedge))
        helper._save(context.mysql_config, systems=2)


def _prepare_link(context):
    from planetgen.db import store

    _save_sectors(context)
    conn = store.get_connection(context.mysql_config)
    try:
        conn.execute("DELETE FROM nearest_systems")      # saving a sector links it: clear the links so the stage has them to make
        conn.commit()
    finally:
        conn.close()


LINK = Recipe(
    "link-neighbors", (("planetgen.generation.run_galaxy", "LINK_UNIT_SECTORS"),), ((1,), (10 ** 6,)),
    prepare=_prepare_link, run_params=lambda context: {"database": _database(context), "since": LONG_AGO, "workers": 4},
    ignore_columns=(("nearest_systems", "id"),),
    note="`nearest_systems.id` is the row's surrogate key, taken in insertion order, and nothing refers to it. Four saved sectors of two systems each, their nearest-system links cleared; the stage links them, a chunk of sectors to a unit.")


def _prepare_settle(context):
    from planetgen.db import store

    _save_sectors(context)
    conn = store.get_connection(context.mysql_config)
    try:
        conn.execute("DELETE FROM sector_paths")      # saving a sector settles it: clear the paths so the stage has them to make
        conn.commit()
    finally:
        conn.close()


SETTLE = Recipe(
    "settle-sector-paths", (("planetgen.generation.run_galaxy", "SETTLE_UNIT_SECTORS"),), ((1,), (10 ** 6,)),
    prepare=_prepare_settle, run_params=lambda context: {"database": _database(context), "since": LONG_AGO, "workers": 4},
    ignore_columns=(("sector_paths", "id"), ("sector_path_knots", "id"), ("sector_path_knots", "path_id")),
    views=(("paths with their knots",
            "SELECT p.sector_id, p.object_table, p.object_id, p.exited, p.duration_years, k.position, k.t_years, k.x_pc, k.y_pc, k.z_pc, k.vx_kms, k.vy_kms, k.vz_kms"
            " FROM sector_paths p JOIN sector_path_knots k ON k.path_id = p.id"),),
    note="`sector_paths.id` and `sector_path_knots.id` are surrogate keys taken in insertion order (which sector's path is saved first), referred to only by `sector_path_knots.path_id`: the pair is compared through the join "
         "of each knot with its path. Four saved sectors of two systems each, their paths cleared; the stage integrates them, a chunk of sectors to a unit.")


# --- fill-map-luminosity: ranges of ids of the scattered phenomena (PERF.106) --------------------------------------

def _prepare_map_luminosity(context):
    from planetgen.db import store
    from tests.test_block_light_parts import _point

    conn = store.get_connection(context.mysql_config)
    try:
        cells = ((10, 1, 2), (12, 0, 1), (13, 3, 3), (14, 2, 5))
        points = [_point(cells[index % len(cells)], 1.0e6 + 1.0e5 * index, -2.0e6 + 7.0e4 * index, 3.0e4 * index, None, index + 1) for index in range(9)]
        store.insert_phenomenon_scatter(conn, points)      # every one has a mass and no map luminosity: the stage fills them
        conn.commit()
    finally:
        conn.close()


MAP_LUMINOSITY = Recipe(
    "fill-map-luminosity", (("planetgen.generation.run_plan", "MAP_LUMINOSITY_BATCH_IDS"),), ((2,), (10 ** 9,)),
    prepare=_prepare_map_luminosity, run_params=lambda context: {"database": _database(context)},
    note="Nine scattered phenomena in four cells with a mass and no map luminosity; the stage fills them, and the cells' light parts with them, a range of ids to a unit.")


# --- the warm-up stages: slices of a listing, chunks of nodes (PERF.180, PERF.257) ----------------------------

def _prepare_warm(context):
    import tempfile

    from planetgen.web import warmup
    from planetgen.web.lib import tilecache

    root = context.extras.get("tile_root") or tempfile.mkdtemp(prefix="cut-neutral-tiles-")
    context.extras["tile_root"] = root
    context.monkeypatch.setenv("PLANETGEN_TILE_CACHE_DIR", root)
    context.monkeypatch.setattr(tilecache, "PRUNE_PROBABILITY", 0.0)
    real = warmup.build_app
    context.monkeypatch.setattr(warmup, "build_app", lambda config=None: real(context.mysql_config))
    _save_sectors(context)


def _reset_warm(context):
    import shutil

    from planetgen.queue import redisclient

    client = redisclient.client(decode_responses=True)
    for key in list(client.scan_iter(f"planetgen:*:{context.mysql_config.database}*")):
        client.delete(key)
    root = context.extras["tile_root"]
    shutil.rmtree(root, ignore_errors=True)
    import os

    os.makedirs(root, exist_ok=True)


def _without_clocks(value):
    """`value` (parsed JSON) with every key that ends `_at` (the time a row was written) left out, at any depth."""
    if isinstance(value, dict):
        return {name: _without_clocks(inner) for name, inner in value.items() if not name.endswith("_at")}
    if isinstance(value, list):
        return [_without_clocks(inner) for inner in value]
    return value


def _capture_warm(context):
    import hashlib
    import json
    import os

    from planetgen.queue import redisclient
    from planetgen.web.lib import tilecache

    root = context.extras["tile_root"]
    files = sorted(f"{os.path.basename(path)} {hashlib.sha256(open(path, 'rb').read()).hexdigest()}" for path, _size, _mtime, _ahead in tilecache.cache_files(root))
    client = redisclient.client(decode_responses=True)
    pages, windows = [], []
    for key in sorted(client.scan_iter(f"planetgen:page:{context.mysql_config.database}:*")):      # the answers; the other keys are clocks, stamps and markers of the run
        fields = client.hgetall(key)
        body = fields.get("b")
        try:
            parsed = json.loads(body)
        except (TypeError, ValueError):
            parsed = None
        shown = json.dumps(_without_clocks(parsed), sort_keys=True) if parsed is not None else repr(body)
        # `c` is the seconds the page took to build and `i` the id of the stamp it was stored under, new each run; `v` is the counters it was stored at
        entry = f"{key} {shown} {fields.get('v')}"
        (windows if isinstance(parsed, dict) and "offset" in parsed and "limit" in parsed else pages).append(entry)
    return {"tile files": files, "redis object answers": pages, "redis listing windows": windows}


def _warm_recipe(kind, family, note):
    return Recipe(
        kind, (("planetgen.web.warmup", "OBJECT_CHUNK"),), ((1,), (10 ** 6,)),
        prepare=_prepare_warm, run_params=lambda context: {"database": context.mysql_config.database, "cap": None},
        reset=_reset_warm, capture=_capture_warm,
        note=f"Four saved sectors of two systems each; the stage builds the {family} pages and scenes into the cache (the tile files and the answers in Redis, read back whole), a slice of the listing to a unit. {note}")


WARM_SECTORS = _warm_recipe("warm-sectors", "sector", "")
WARM_SYSTEMS = _warm_recipe("warm-systems", "system", "")
WARM_PHENOMENA = _warm_recipe("warm-phenomena", "phenomenon", "")


DEEP_SLICE = 40
"""int: How many of the zoom run's nodes the deep-node recipe warms (the whole run is thousands)."""


def _prepare_warm_nodes(context):
    import types

    from planetgen.db import block_stats, store
    from planetgen.web import warmup

    _prepare_warm(context)
    context.monkeypatch.setattr(warmup, "get_settings", lambda: types.SimpleNamespace(tile_cache=types.SimpleNamespace(warm_zoom_floor=81)))
    whole = warmup.zoom_run_nodes
    context.monkeypatch.setattr(warmup, "zoom_run_nodes", lambda *args: whole(*args)[:DEEP_SLICE])
    conn = store.get_connection(context.mysql_config)
    try:
        block_stats.rebuild(conn)           # the blocks over the saved sectors, which `warm-galaxy-nodes` lists
        conn.commit()
    finally:
        conn.close()


def _node_recipe(kind, what):
    return Recipe(
        kind, (("planetgen.web.warmup", "NODE_CHUNK"),), ((4,), (10 ** 6,)),
        prepare=_prepare_warm_nodes, run_params=lambda context: {"database": context.mysql_config.database, "cap": None},
        reset=_reset_warm, capture=_capture_warm,
        note=f"Four saved sectors of two systems each, the block totals rolled up; the stage builds {what} into the tile cache, a chunk of nodes to a unit.")


WARM_NODES = _node_recipe("warm-galaxy-nodes", "the populated nodes of the built levels")
WARM_DEEP_NODES = _node_recipe("warm-galaxy-deep-nodes", f"the first {DEEP_SLICE} nodes of the zoom run")


RECIPES = {recipe.kind: recipe for recipe in (LINK, SETTLE, MAP_LUMINOSITY, WARM_SECTORS, WARM_SYSTEMS, WARM_PHENOMENA, WARM_NODES, WARM_DEEP_NODES)}
"""dict: Command kind -> `Recipe`."""
