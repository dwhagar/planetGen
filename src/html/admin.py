#!/usr/bin/env python3
# html/admin.py

"""
Admin landing page: API key management (list/create/revoke -- `GET`/
`POST /api/auth/api-keys`, `DELETE /api/auth/api-keys/<id>`). Sector/
system creation, modification, and deletion are done directly against the
JSON API (`docs/api.md`'s "Write endpoints"), authenticated with one of
these keys -- this page doesn't duplicate that as a web form, only
manages the keys themselves and a caller's own admin identity.

Requires an authenticated session with fresh (non-default) credentials --
redirects to `login.py` (not logged in) or `changecreds.py` (still on the
seeded `admin`/`password` bootstrap row) otherwise, same gate
`authz.require_admin(fresh=True)` enforces server-side for every write
route this page's keys would be used against.
"""

import os
import sys

_HTML_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HTML_DIR, "lib"))

from apiclient import ApiError, auth_create_api_key, auth_list_api_keys, auth_me, auth_revoke_api_key
from fmt import esc
from page import form_params, incoming_cookie_header, redirect, render, render_error


def _keys_table_html(keys):
    if not keys:
        return '<p class="hint">No API keys yet.</p>'

    rows = []
    for key in keys:
        if key["revoked_at"]:
            status = f'<span class="badge">revoked {esc(key["revoked_at"])}</span>'
            action = ""
        else:
            status = '<span class="badge">active</span>'
            action = f"""
<form method="post" action="admin.py" class="table-form">
  <input type="hidden" name="action" value="revoke_key">
  <input type="hidden" name="key_id" value="{key['id']}">
  <button type="submit" class="btn">Revoke</button>
</form>
"""
        rows.append(
            "<tr>"
            f"<td>{esc(key['label'])}</td>"
            f"<td>{esc(key['created_at'] or '')}</td>"
            f"<td>{esc(key['last_used_at'] or 'never')}</td>"
            f"<td>{status}</td>"
            f"<td>{action}</td>"
            "</tr>"
        )
    return f"""
<div class="table-scroll"><table>
  <thead><tr><th>Label</th><th>Created</th><th>Last used</th><th>Status</th><th></th></tr></thead>
  <tbody>{''.join(rows)}</tbody>
</table></div>
"""


def _page_html(identity, keys, new_key=None, error=None):
    error_html = f'<p class="error">{esc(error)}</p>' if error else ""
    new_key_html = ""
    if new_key is not None:
        new_key_html = f"""
<section class="panel">
<h2>New API Key: {esc(new_key['label'])}</h2>
<p class="hint">Copy this now -- it will not be shown again.</p>
<p class="api-key-value">{esc(new_key['key'])}</p>
</section>
"""
    return f"""
<section class="panel">
<h2>Signed in</h2>
<p><strong>{esc(identity['username'])}</strong> &mdash;
<a href="changecreds.py">change username/password</a> &mdash;
<a href="logout.py">log out</a></p>
</section>
{new_key_html}
{error_html}
<section class="panel">
<h2>API Keys</h2>
<p class="hint">Used as an <code>Authorization: Bearer &lt;key&gt;</code>
header against the write endpoints documented in the API reference
(sector/system create, update, delete).</p>
{_keys_table_html(keys)}
<form method="post" action="admin.py" class="search-form">
  <input type="hidden" name="action" value="create_key">
  <div class="search-fields">
    <label class="search-field">Label
      <input type="text" name="label" placeholder="e.g. galaxy-gen script" required>
    </label>
  </div>
  <div class="search-actions">
    <button type="submit" class="btn">Create key</button>
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
if identity["must_change_credentials"]:
    redirect("changecreds.py")
    sys.exit(0)

new_key = None
error = None

if os.environ.get("REQUEST_METHOD", "GET").upper() == "POST":
    fields = form_params()
    action = fields.get("action")
    try:
        if action == "create_key":
            label = fields.get("label", "").strip()
            if not label:
                error = "Label is required."
            else:
                new_key = auth_create_api_key(cookie_header, label)
        elif action == "revoke_key":
            auth_revoke_api_key(cookie_header, int(fields.get("key_id", "0")))
        else:
            error = "Unrecognized form action."
    except ApiError as exc:
        error = str(exc)
    except ValueError:
        error = "Invalid key id."

try:
    keys = auth_list_api_keys(cookie_header)
except ApiError as exc:
    render_error(f"Could not reach the planetGen API ({exc}).", status="502 Bad Gateway")

render("Admin", _page_html(identity, keys, new_key=new_key, error=error))
