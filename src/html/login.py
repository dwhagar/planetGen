#!/usr/bin/env python3
# html/login.py

"""
Admin login form -- `POST /api/auth/login` via `apiclient.auth_login`,
relaying the session cookie it sets back to the browser (see
`page.py`'s `send_headers`/`render`'s `set_cookie_headers` and
`apiclient._auth_request`'s docstring for why this module never
constructs a cookie value itself, only relays the API's).

Not built on `page.run()` (every other page in `html/` is) -- a
successful login needs a redirect (to `changecreds.py` if the admin is
still on the seeded default credentials, else `admin.py`), not a
rendered page, and a wrong username/password needs the *same* form
re-rendered with an inline error rather than `run()`'s generic 502/500
error page. Both cases are handled directly here instead.
"""

import os
import sys

_HTML_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HTML_DIR, "lib"))

from apiclient import ApiError, auth_login
from fmt import esc
from page import form_params, redirect, render, render_error


def _form_html(error=None):
    error_html = f'<p class="error">{esc(error)}</p>' if error else ""
    return f"""
<section class="panel">
<h2>Admin Login</h2>
{error_html}
<form method="post" action="login.py" class="search-form">
  <div class="search-fields">
    <label class="search-field">Username
      <input type="text" name="username" autocomplete="username" required autofocus>
    </label>
    <label class="search-field">Password
      <input type="password" name="password" autocomplete="current-password" required>
    </label>
  </div>
  <div class="search-actions">
    <button type="submit" class="btn">Log in</button>
  </div>
</form>
</section>
"""


method = os.environ.get("REQUEST_METHOD", "GET").upper()

if method != "POST":
    render("Admin Login", _form_html())
else:
    fields = form_params()
    username = fields.get("username", "")
    password = fields.get("password", "")

    try:
        result, set_cookie_headers = auth_login(username, password)
    except ApiError as exc:
        if exc.status_code == 401:
            render("Admin Login", _form_html(error="Invalid username or password."))
            sys.exit(0)
        render_error(f"Could not reach the planetGen API ({exc}).", status="502 Bad Gateway")

    destination = "changecreds.py" if result["must_change_credentials"] else "admin.py"
    redirect(destination, set_cookie_headers=set_cookie_headers)
