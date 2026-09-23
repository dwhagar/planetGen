#!/usr/bin/env python3
# html/galaxy_view.py

"""
Interactive 3D Galaxy Map live-viewport JSON endpoint.

Every other script under `html/` is a page: server-rendered once, its own
data baked in at render time, never fetched again by the browser's own
script. This one is the single exception -- `static/galaxymap3d.js`
(loaded by `galaxy.py`) calls this directly via `fetch()`, debounced,
every time its 3D camera moves, since the whole galaxy's own placed/
planned/density content can never be baked into one page load the way a
sector's handful of systems can (`lib/starmap.py`'s own one-shot JSON
block). This script is the browser-reachable, same-origin proxy for
`GET /api/galaxy/view` (`queryDb.galaxy_view`) -- the browser never calls
the API directly (same reasoning as every other page here: the API's own
reachable address, `PLANETGEN_API_BASE_URL`, is a server-side deployment
setting that may not even be reachable from a visitor's own browser).

Returns raw JSON (`page.run_json`), not an HTML page -- `cx`/`cy`/`cz`
(the view center, galaxy-frame parsecs) and `radius_pc` (the view radius)
come from the GET query string (`page.query_params`), not `nav_params`'s
own POST-preferred convention: this is never navigated to or bookmarked
by a person, only ever fetched by this page's own script, so the
privacy/history reasoning that makes every *navigational* link in `html/`
post its params instead (see `lib/page.py`'s own module docstring)
doesn't apply here -- a `fetch()` URL never appears in the browser's
address bar, history, or an outgoing `Referer` header regardless of
method.
"""

import os
import sys

_HTML_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HTML_DIR, "lib"))

from apiclient import get_galaxy_view
from page import query_params, render_json_error, run_json


def handler():
    params = query_params()
    db_name = params.get("db", "")

    try:
        cx = float(params["cx"])
        cy = float(params["cy"])
        cz = float(params["cz"])
        radius_pc = float(params["radius_pc"])
    except KeyError:
        # A malformed request from this script's own client-side JS, not a
        # planetGen API failure -- reported directly (400) rather than
        # raised as `apiclient.ApiError` (which `page.run_json` maps to
        # 502, meaning "the API itself failed/is unreachable"), same
        # distinction `html/api/routes.py`'s own query-param validation
        # draws for the server-side `/api/galaxy/view` route this proxies.
        render_json_error("cx/cy/cz/radius_pc query parameters are all required.")
    except ValueError:
        render_json_error("cx/cy/cz/radius_pc must all be numbers.")

    return get_galaxy_view(db_name, cx, cy, cz, radius_pc)


run_json(handler)
