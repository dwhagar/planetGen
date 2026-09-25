#!/usr/bin/env python3
# html/search.py

"""
Moved: the search page is now served by the Flask app at `/search`
(`web/views.py`'s `search`, `web/searchpage.py`), where every filter is a
GET parameter. This shim 301-redirects an old `search.py` link, form post
or bookmark there, carrying over every search parameter the old page took
(name fields, size ranges, repeated tag facets and per-panel page
numbers) from the GET query or POST body. `db` is dropped: the Flask
pages take the database from config. Removed in the cleanup PR, once
nothing links here any more.
"""

import os
import sys
from urllib.parse import urlencode

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))

from page import nav_multi_params, redirect, start_request_log  # noqa: E402

CARRIED = (
    ("q", "sector_q", "system_q", "star_q", "planet_q", "moon_q")
    + tuple(f"{entity}_{bound}_radius_km" for entity in ("star", "planet", "moon") for bound in ("min", "max"))
    + ("type", "spectral", "luminosity", "class", "body", "life", "moon_class", "moon_body", "moon_life", "density")
    + tuple(f"{panel}_page" for panel in ("sectors", "systems", "stars", "planets", "moons", "belts"))
)
"""tuple: The old page's parameters, in the order `/search` lists them."""


def new_url(params):
    """`/search?...` for the old request's parameters (`parse_qs`
    shape), keeping only non-empty values of `CARRIED` names."""
    pairs = [(name, value.strip()) for name in CARRIED for value in params.get(name, ()) if value.strip()]
    return f"/search?{urlencode(pairs)}" if pairs else "/search"


if __name__ == "__main__":
    start_request_log()
    redirect(new_url(nav_multi_params()), status="301 Moved Permanently")
