#!/usr/bin/env python3
# html/logout.py

"""
Ends the current admin session -- `POST /api/auth/logout` (via
`apiclient.auth_logout`, forwarding the browser's own session cookie
server-side) and relays the cleared cookie back. A plain link
(`page.py`'s sidenav "Logout" item) rather than a form -- there's nothing
to confirm and no request body needed, so the extra step of a POST-only
form button buys nothing here (unlike the destructive/credential-
changing actions on `admin.py`/`changecreds.py`, which do use POST forms).
"""

import os
import sys

_HTML_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HTML_DIR, "lib"))

from apiclient import ApiError, auth_logout
from page import incoming_cookie_header, redirect, render_error

try:
    set_cookie_headers = auth_logout(incoming_cookie_header())
except ApiError as exc:
    render_error(f"Could not reach the planetGen API ({exc}).", status="502 Bad Gateway")

redirect("index.py", set_cookie_headers=set_cookie_headers)
