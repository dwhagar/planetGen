#!/usr/bin/env python3
# html/changecreds.py

"""
Change the logged-in admin's username/password -- `POST
/api/auth/change-credentials` via `apiclient.auth_change_credentials`.

Reached two ways: forced (`login.py` redirects here straight after a
successful login while `must_change_credentials` is still set on the
seeded `admin`/`password` bootstrap row -- see
`stellarObjects/adminAuth.py`'s module docstring) and voluntary (an
already-"fresh" admin can still come here to rotate credentials; the
sidenav has no direct link for this today, but the URL works either way
-- both cases render the exact same form).

Not built on `page.run()`, same reasoning as `login.py`: a successful
change needs a redirect carrying the freshly re-issued session cookie,
and a policy/wrong-password failure needs the same form re-rendered with
an inline error, not `run()`'s generic error page.
"""

import os
import sys

_HTML_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HTML_DIR, "lib"))

from apiclient import ApiError, auth_change_credentials, auth_me
from fmt import esc
from page import form_params, incoming_cookie_header, redirect, render, render_error


def _form_html(current_username, error=None):
    error_html = f'<p class="error">{esc(error)}</p>' if error else ""
    return f"""
<section class="panel">
<h2>Change Credentials</h2>
<p>Logged in as <strong>{esc(current_username)}</strong>.</p>
{error_html}
<form method="post" action="changecreds.py" class="search-form">
  <div class="search-fields">
    <label class="search-field">Current password
      <input type="password" name="current_password" autocomplete="current-password" required>
    </label>
    <label class="search-field">New username
      <input type="text" name="new_username" autocomplete="username" value="{esc(current_username)}" required>
    </label>
    <label class="search-field">New password
      <input type="password" name="new_password" autocomplete="new-password" required minlength="12">
    </label>
  </div>
  <div class="search-actions">
    <button type="submit" class="btn">Change credentials</button>
  </div>
</form>
</section>
"""


cookie_header = incoming_cookie_header()
try:
    identity = auth_me(cookie_header)
except ApiError as exc:
    render_error(f"Could not reach the planetGen API ({exc}).", status="502 Bad Gateway")
if identity is None:
    redirect("login.py")
    sys.exit(0)

method = os.environ.get("REQUEST_METHOD", "GET").upper()

if method != "POST":
    render("Change Credentials", _form_html(identity["username"]))
else:
    fields = form_params()
    try:
        _result, set_cookie_headers = auth_change_credentials(
            cookie_header,
            current_password=fields.get("current_password", ""),
            new_username=fields.get("new_username", ""),
            new_password=fields.get("new_password", ""),
        )
    except ApiError as exc:
        if exc.status_code in (400, 401):
            render("Change Credentials", _form_html(identity["username"], error=str(exc)))
            sys.exit(0)
        render_error(f"Could not reach the planetGen API ({exc}).", status="502 Bad Gateway")

    redirect("admin.py", set_cookie_headers=set_cookie_headers)
