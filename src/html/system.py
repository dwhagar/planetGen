#!/usr/bin/env python3
# html/system.py

"""
Moved: the system page is now served by the Flask app at `/system/<id>`
(`web/system_pages.py`'s `system`). This shim 301-redirects an old
`system.py` link, bookmark or `post_link` form there, reading `id` (and
the `code=wikitext|markdown` view) from the old GET query or POST body.
The `db` parameter is dropped: the Flask pages take the database from
config. A missing or non-numeric id goes to `/systems`. The old admin
"Upload to Wiki" POST is not replayed (the form now posts to the new
page). Removed in the cleanup PR, once nothing links here any more.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))

from page import moved_permanently, nav_params  # noqa: E402


def new_path(params):
    """The Flask URL for the old request's parameters."""
    system_id = params.get("id", "")
    if not system_id.isdigit():
        return "/systems"
    code = params.get("code")
    query = f"?code={code}" if code in ("wikitext", "markdown") else ""
    return f"/system/{int(system_id)}{query}"


if __name__ == "__main__":
    moved_permanently(new_path(nav_params()))
