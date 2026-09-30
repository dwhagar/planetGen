#!/usr/bin/env python3
# html/logout.py

"""
Moved: logging out is now `/logout` on the Flask app (`web/admin_pages.py`):
`GET /logout` shows a confirm button and changes nothing; the button
POSTs (with a CSRF token) to end the session. This shim no longer logs
anyone out, it only sends an old link or bookmark to that page.
This shim 301-redirects an old link, form post or bookmark there. A form
POSTed here is not replayed: the visitor lands on the new page's form.
Removed in the cleanup PR, once nothing links here any more.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))

from page import moved_permanently  # noqa: E402

if __name__ == "__main__":
    moved_permanently("/logout")
