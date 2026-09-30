#!/usr/bin/env python3
# html/phenomenon.py

"""
Moved: the phenomenon page is now served by the Flask app at
`/phenomenon/<type>/<id>` (`web/system_pages.py`'s `phenomenon`). This
shim 301-redirects an old `phenomenon.py` link, bookmark, `post_link`
form or map marker there, reading `type` and `id` from the old GET query
or POST body. The `db` parameter is dropped: the Flask pages take the
database from config. A missing or malformed type or id goes to
`/phenomena`. Removed in the cleanup PR, once nothing links here any
more.
"""

import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))

from page import moved_permanently, nav_params  # noqa: E402


def new_path(params):
    """The Flask URL for the old request's parameters."""
    phenomenon_type = params.get("type", "")
    phenomenon_id = params.get("id", "")
    if not re.fullmatch(r"[a-z_]{1,40}", phenomenon_type) or not phenomenon_id.isdigit():
        return "/phenomena"
    return f"/phenomenon/{phenomenon_type}/{int(phenomenon_id)}"


if __name__ == "__main__":
    moved_permanently(new_path(nav_params()))
