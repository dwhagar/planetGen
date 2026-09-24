#!/usr/bin/env python3
# html/admin.py

"""
Moved: the admin page (API keys, a sector's manual wiki link) is now
`/admin` on the Flask app (`web/admin_pages.py`), keeping `keys_page`.
This shim 301-redirects an old link, form post or bookmark there. A form
POSTed here is not replayed: the visitor lands on the new page's form.
Removed in the cleanup PR, once nothing links here any more.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))

from page import moved_permanently  # noqa: E402

if __name__ == "__main__":
    moved_permanently("/admin", ("keys_page",))
