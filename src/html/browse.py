#!/usr/bin/env python3
# html/browse.py

"""
Moved: the sector and standalone-system lists are now served by the
Flask app at `/` (both tables, each paged with `?sectors_page=N`/
`?standalone_page=N`), `/sectors` and `/systems` (`web/views.py`). This
shim 301-redirects an old `browse.py` link, form post or bookmark to `/`,
keeping both page numbers; a link that targeted `browse.py#sectors` or
`#standalone-systems` keeps its fragment, since browsers carry a
fragment across a redirect whose `Location` has none. The `db` parameter
is dropped: the Flask pages take the database from config. Removed in
the cleanup PR, once nothing links here any more.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))

from page import moved_permanently  # noqa: E402

if __name__ == "__main__":
    moved_permanently("/", ("sectors_page", "standalone_page"))
