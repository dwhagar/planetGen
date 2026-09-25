#!/usr/bin/env python3
# html/sector.py

"""
Moved: the sector page is now served by the Flask app at `/sector/<id>`
(`web/sector_page.py`). This shim 301-redirects an old `sector.py` link,
bookmark or form post there, keeping `contents_page` from the old GET
query or POST body. The `db` parameter is dropped (the Flask pages take
the database from config), and an old admin form post is not replayed:
the visitor lands on the new page and submits it from there. Without a
valid `id` it goes to `/sectors`. Removed in the cleanup PR.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))

from page import moved_permanently, nav_params  # noqa: E402

if __name__ == "__main__":
    sector_id = nav_params().get("id", "").strip()
    if sector_id.isdigit() and int(sector_id) > 0:
        moved_permanently(f"/sector/{int(sector_id)}", ("contents_page",))
    else:
        moved_permanently("/sectors")
