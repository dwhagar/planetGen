# html/web/old_urls.py

"""
Old `/<name>.py` URLs (the CGI pages before they moved to Flask) send an
old bookmark or link to the page that replaced it with a 301, so it still
lands somewhere sensible. GET only: an old form POST gets the CSRF check's
400 like any other stale form.

The query string is carried over as it is, minus `db` (the pages take the
database from config), since the new pages kept the old parameter names
(`sectors_page`, `quadrant`, the search filters, ...); empty values are
dropped. The detail pages
turn their `id` (and a phenomenon's `type`) into the new path, falling back
to the bare list page without a valid one. An unknown name gets the normal 404 page.
"""

import re
from urllib.parse import urlencode

from flask import abort, redirect, request, url_for

from . import bp

OLD_PAGES = {
    "index": "web.index",
    "browse": "web.index",
    "galaxy": "web.galaxy",
    "galaxy3d": "web.galaxy",
    "galaxy_view": "web.galaxy",
    "galaxy_tiles": "web.galaxy_tiles",
    "sector": "web.sectors",
    "system": "web.systems",
    "phenomena": "web.phenomena",
    "phenomenon": "web.phenomena",
    "nav": "web.nav",
    "search": "web.search",
    "login": "web.login",
    "logout": "web.logout",
    "changecreds": "web.account",
    "admin": "web.admin",
    "adminstats": "web.admin_stats",
}
"""dict: Old script name (without `.py`) -> the endpoint that replaced it."""

_TYPE = re.compile(r"^[a-z_]{1,40}$")


def _nav_query(args):
    """The old NAV picker's `from=12&from_kind=phenomenon&from_type=nebula`
    becomes `/nav`'s `from=nebula:12` (a bare id is a system)."""
    query = []
    for prefix in ("from", "to"):
        raw = args.get(prefix, "").strip()
        if not raw.isdigit():
            continue
        kind = args.get(f"{prefix}_type", "")
        if args.get(f"{prefix}_kind") != "phenomenon" or not _TYPE.match(kind):
            kind = "system"
        query.append((prefix, f"{kind}:{int(raw)}"))
    query += [(name, args[name]) for name in ("from_sector", "to_sector") if args.get(name, "").isdigit()]
    return query


def new_url(name, args):
    """The new URL for old script `name` with query `args` (a MultiDict),
    or `None` for a name that never was a page."""
    endpoint = OLD_PAGES.get(name)
    if endpoint is None:
        return None
    raw_id = args.get("id", "").strip()
    query = [(key, value) for key, values in args.lists() if key not in ("db", "id")
             for value in values if value.strip()]
    if name == "nav":
        query = _nav_query(args)
    elif name in ("sector", "system") and raw_id.isdigit() and int(raw_id) > 0:
        endpoint = f"web.{name}"
        query = [(key, value) for key, value in query if key in ("contents_page", "code")]
        return url_for(endpoint, **{f"{name}_id": int(raw_id)}) + (f"?{urlencode(query)}" if query else "")
    elif name == "phenomenon" and raw_id.isdigit() and _TYPE.match(args.get("type", "")):
        return url_for("web.phenomenon", phenomenon_type=args["type"], phenomenon_id=int(raw_id))
    elif name in ("sector", "system", "phenomenon"):
        query = []  # no valid id: the bare list page
    url = url_for(endpoint)
    return f"{url}?{urlencode(query, safe=':')}" if query else url


@bp.route("/<name>.py")
def old_page(name):
    """301 from an old CGI page's URL to the page that replaced it."""
    url = new_url(name, request.args)
    if url is None:
        abort(404)
    return redirect(url, code=301)
