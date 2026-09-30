#!/usr/bin/env python3
# html/phenomena.py

"""
Moved: the phenomena list is now served by the Flask app at `/phenomena`
(`web/system_pages.py`'s `phenomena`). This shim 301-redirects an old
`phenomena.py` link or bookmark there, keeping the page number (`page`)
from the old GET query or POST body. The `db` parameter is dropped: the
Flask pages take the database from config. Removed in the cleanup PR,
once nothing links here any more.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))

from page import moved_permanently  # noqa: E402

if __name__ == "__main__":
    moved_permanently("/phenomena", ("page",))
