#!/usr/bin/env python3
# html/admin.py

"""
Admin landing page (server health and database stats have their own
page, `adminstats.py`, linked from here): API key management (list/create/revoke -- `GET`/
`POST /api/auth/api-keys`, `DELETE /api/auth/api-keys/<id>`). Sector/
system creation, modification, and deletion are done directly against the
JSON API (`docs/api.md`'s "Write endpoints"), authenticated with one of
these keys -- this page doesn't duplicate that as a web form, only
manages the keys themselves and a caller's own admin identity.

The one exception is a sector's wiki link (`sectors.wiki_url`, see
`schema.sql`'s "v22" header note): a small form below lets an admin set
or clear it directly (`PATCH /api/sectors/<id>` `{"wiki_url": ...}`, via
the authenticated session cookie already established here, not an API
key) for a sector that already has a hand-written wiki page from outside
this app, or to correct/remove a link `html/sector.py`'s own "Upload to
Wiki" button set. Since this page isn't scoped to one content database
the way `html/sector.py`/`html/system.py` are, the form asks for the
target database by name alongside the sector id.

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

from apiclient import (
    ApiError, admin_set_sector_wiki_url, auth_create_api_key, auth_list_api_keys, auth_me, auth_revoke_api_key,
)
from fmt import esc
from page import form_params, incoming_cookie_header, redirect, render, render_error
from pagination import page_slice, parse_page, render_pagination


def _keys_table_html(keys, keys_page=1):
    if not keys:
        return '<p class="hint">No API keys yet.</p>'

    page_keys, keys_page = page_slice(keys, keys_page)
    rows = []
    for key in page_keys:
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
{render_pagination("admin.py", {}, "keys_page", keys_page, len(keys), anchor="api-keys", label="API key pages")}
"""


def _page_html(identity, keys, new_key=None, error=None, wiki_message=None, wiki_error=None, keys_page=1):
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
    wiki_message_html = f'<p class="hint">{esc(wiki_message)}</p>' if wiki_message else ""
    wiki_error_html = f'<p class="error">{esc(wiki_error)}</p>' if wiki_error else ""
    return f"""
<section class="panel">
<h2>Signed in</h2>
<p><strong>{esc(identity['username'])}</strong> &mdash;
<a href="changecreds.py">change username/password</a> &mdash;
<a href="adminstats.py">server and database stats</a> &mdash;
<a href="logout.py">log out</a></p>
</section>
{new_key_html}
{error_html}
<section class="panel" id="api-keys">
<h2>API Keys</h2>
<p class="hint">Used as an <code>Authorization: Bearer &lt;key&gt;</code>
header against the write endpoints documented in the API reference
(sector/system create, update, delete).</p>
{_keys_table_html(keys, keys_page)}
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
<section class="panel">
<h2>Sector Wiki Link</h2>
<p class="hint">Manually set or clear a sector's wiki link (leave the URL
blank to clear it) -- for a sector with a hand-written page from outside
this app, or to fix a link <a href="sector.py">a sector's own "Upload to
Wiki" button</a> set.</p>
{wiki_message_html}{wiki_error_html}
<form method="post" action="admin.py" class="search-form">
  <input type="hidden" name="action" value="set_sector_wiki_url">
  <div class="search-fields">
    <label class="search-field">Database
      <input type="text" name="db" placeholder="e.g. planetgen" required>
    </label>
    <label class="search-field">Sector ID
      <input type="number" name="sector_id" min="1" required>
    </label>
    <label class="search-field">Wiki URL
      <input type="text" name="wiki_url" placeholder="https://wiki.example.com/Sector_Name">
    </label>
  </div>
  <div class="search-actions">
    <button type="submit" class="btn">Save link</button>
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
wiki_message = None
wiki_error = None
keys_page = 1

if os.environ.get("REQUEST_METHOD", "GET").upper() == "POST":
    fields = form_params()
    action = fields.get("action")
    # The key table's pager posts just `keys_page`, no action.
    keys_page = parse_page(fields.get("keys_page"))
    try:
        if action == "create_key":
            label = fields.get("label", "").strip()
            if not label:
                error = "Label is required."
            else:
                new_key = auth_create_api_key(cookie_header, label)
        elif action == "revoke_key":
            auth_revoke_api_key(cookie_header, int(fields.get("key_id", "0")))
        elif action == "set_sector_wiki_url":
            db_name = fields.get("db", "").strip()
            sector_id_raw = fields.get("sector_id", "").strip()
            wiki_url = fields.get("wiki_url", "").strip() or None
            if not db_name or not sector_id_raw:
                wiki_error = "Database and Sector ID are required."
            else:
                try:
                    sector_id = int(sector_id_raw)
                except ValueError:
                    wiki_error = "Invalid sector id."
                else:
                    admin_set_sector_wiki_url(cookie_header, db_name, sector_id, wiki_url)
                    wiki_message = (
                        f"Wiki link for sector {sector_id} in {db_name!r} "
                        f"{'cleared' if wiki_url is None else f'set to {wiki_url}'}."
                    )
        elif action:
            error = "Unrecognized form action."
    except ApiError as exc:
        if action == "set_sector_wiki_url":
            wiki_error = str(exc)
        else:
            error = str(exc)
    except ValueError:
        error = "Invalid key id."

try:
    keys = auth_list_api_keys(cookie_header)
except ApiError as exc:
    render_error(f"Could not reach the planetGen API ({exc}).", status="502 Bad Gateway")

render("Admin", _page_html(identity, keys, new_key=new_key, error=error, wiki_message=wiki_message,
                           wiki_error=wiki_error, keys_page=keys_page))
