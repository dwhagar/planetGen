#!/usr/bin/env python3
# html/galaxy_tiles.py

"""
Moved: the Galaxy Map's tile JSON is now served by the Flask app at
`/galaxy/tiles` (`web/galaxy_views.py`). This shim 301-redirects a fetch
from a Galaxy Map page still open from before the move there, keeping
`tiles`, `density` and `stamp` and dropping `db` (the Flask app takes the
database from config). Removed in the cleanup PR.
"""

import os
import sys
from urllib.parse import urlencode

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))

from page import query_params, redirect, start_request_log  # noqa: E402

KEEP = ("tiles", "density", "stamp")


def new_url(params):
    """`/galaxy/tiles` with the tile request's own parameters."""
    carried = [(name, params[name]) for name in KEEP if params.get(name)]
    return f"/galaxy/tiles?{urlencode(carried)}" if carried else "/galaxy/tiles"


if __name__ == "__main__":
    start_request_log()
    redirect(new_url(query_params()), status="301 Moved Permanently")
