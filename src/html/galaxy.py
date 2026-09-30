#!/usr/bin/env python3
# html/galaxy.py

"""
Moved: the Galaxy Map is now served by the Flask app at `/galaxy`
(`web/galaxy_views.py`). This shim 301-redirects an old `galaxy.py` link
or bookmark there, keeping `quadrant` (`I`-`IV`) and `page` from the old
GET query or POST body. The `db` parameter is dropped: the Flask pages
take the database from config. Removed in the cleanup PR, once nothing
links here any more.
"""

import os
import sys
from urllib.parse import urlencode

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))

from page import nav_params, redirect, start_request_log  # noqa: E402

QUADRANTS = ("I", "II", "III", "IV")


def new_url(params):
    """`/galaxy` plus whichever old parameters still mean something."""
    carried = []
    quadrant = (params.get("quadrant") or "").strip().upper()
    if quadrant in QUADRANTS:
        carried.append(("quadrant", quadrant))
        page = params.get("page") or ""
        if page.isdigit() and int(page) > 1:
            carried.append(("page", page))
    return f"/galaxy?{urlencode(carried)}" if carried else "/galaxy"


if __name__ == "__main__":
    start_request_log()
    redirect(new_url(nav_params()), status="301 Moved Permanently")
