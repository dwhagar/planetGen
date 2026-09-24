#!/usr/bin/env python3
# html/galaxy_tiles.py

"""
Interactive 3D Galaxy Map tile endpoint -- the JSON `static/
galaxymap3d.js` fetches as its camera moves, replacing the old `galaxy_view.py`'s
"everything within R of this point" query with fixed cubes of space (see
`stellarObjects.galaxyViewport`'s "Cube tiles" section).

`tiles` is a comma-separated list of `level/ix/iy/iz` keys and `density`
an optional single key to anchor the illustrative density cloud on. Both
come from the GET query string rather than `nav_params`'s POST
convention: this is only ever fetched by the page's own script, never
navigated to, so a `fetch()` URL never reaches the address bar, history
or a `Referer` header.

Every tile goes through `lib/tilecache.py`'s disk cache first; only tiles
not already cached reach the API. The response carries the database's
content `stamp` and cache `generation`; `stamp` is the one the browser's
own tile cache is at, so when that's out of date the response also lists
which tiles changed since (`history`), and the browser drops only those.
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
    known_stamp = params.get("stamp") or None

    try:
        return fetch_tiles(db_name, tile_keys, density_key, known_stamp)
    except TileRequestError as exc:
        # The page's own script sent a bad request -- a 400, not the 502
        # `run_json` gives an API failure.
        render_json_error(str(exc))


run_json(handler)
