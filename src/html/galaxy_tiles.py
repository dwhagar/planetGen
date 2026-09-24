#!/usr/bin/env python3
# html/galaxy_tiles.py

"""
Interactive 3D Galaxy Map tile endpoint -- the JSON `static/
galaxymap3d.js` fetches as its camera moves, replacing `galaxy_view.py`'s
"everything within R of this point" query with fixed cubes of space (see
`stellarObjects.galaxyViewport`'s "Cube tiles" section).

`tiles` is a comma-separated list of `level/ix/iy/iz` keys and `density`
an optional single key to anchor the illustrative density cloud on. Both
come from the GET query string for the same reason `galaxy_view.py`'s
parameters do (see that script's docstring): this is only ever fetched by
the page's own script, never navigated to.

Every tile goes through `lib/tilecache.py`'s disk cache first; only tiles
not already cached reach the API. The response carries the database's
content `stamp`, which the browser keys its own tile cache by.
"""

import os
import sys

_HTML_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HTML_DIR, "lib"))

from page import query_params, render_json_error, run_json
from tilecache import TileRequestError, fetch_tiles


def handler():
    params = query_params()
    db_name = params.get("db", "")
    tile_keys = [key for key in (params.get("tiles") or "").split(",") if key]
    density_key = params.get("density") or None

    try:
        return fetch_tiles(db_name, tile_keys, density_key)
    except TileRequestError as exc:
        # The page's own script sent a bad request -- a 400, not the 502
        # `run_json` gives an API failure (same distinction
        # `galaxy_view.py` draws).
        render_json_error(str(exc))


run_json(handler)
