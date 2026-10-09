# planetgen/web/warmup.py

"""
Builds the Galaxy Map's opening view ahead of time (MAP.134): the tiles
and the galaxy stage the first frame of `/galaxy` asks for, written to the
tile cache (`web/lib/tilecache.py`) so the first visitor after a release,
or after the cache was cleared, finds them already there instead of
waiting for the database (the bright-star sample, MAP.47, is read once
per cold tile request).

`opening_request` is the one place that says which tiles those are; the
`/galaxy` page asks for exactly them, so the warm-up and the page can't
drift apart. `python -m planetgen.cli.warm_map` runs it (update.sh does,
in the background).
"""

import time

from planetgen.physics.units import ly_to_pc
from planetgen.tuning import DEFAULT_SECTOR_EDGE_LY
from planetgen.web.lib import apiclient
from planetgen.web.lib.tilecache import fetch_stage, fetch_tiles
from planetgen.web.maps.galaxymap3d import initial_tile_request, view_radius_bounds


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
    return initial_tile_request(max_radius)


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
