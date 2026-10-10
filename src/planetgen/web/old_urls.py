# planetgen/web/old_urls.py

"""
Old `/<name>.py` URLs (the CGI pages before they moved to Flask) send an
old bookmark or link to the page that replaced it with a 301, so it still
lands somewhere sensible. GET only: an old form POST gets the CSRF check's
400 like any other stale form.

The query string is carried over as it is, minus `db` (the pages take the
database from config), since the new pages kept the old parameter names
(`sectors_page`, `quadrant`, the search filters, ...); empty values are
dropped. The detail pages and the NAV page go to their bare list page: the
row numbers their old URLs carried name nothing since objects have IDs
(API.23). An unknown name gets the normal 404 page.
"""

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


def new_url(name, args):
    """The new URL for old script `name` with query `args` (a MultiDict),
    or `None` for a name that never was a page."""
    endpoint = OLD_PAGES.get(name)
    if endpoint is None:
        return None
    query = [(key, value) for key, values in args.lists() if key not in ("db", "id")
             for value in values if value.strip()]
    if name in ("nav", "sector", "system", "phenomenon"):
        query = []  # the old row numbers name nothing now (API.23): the page's bare list
    url = url_for(endpoint)
    return f"{url}?{urlencode(query, safe=':')}" if query else url


@bp.route("/<name>.py")
def old_page(name):
    """301 from an old CGI page's URL to the page that replaced it."""
    url = new_url(name, request.args)
    if url is None:
        abort(404)
    return redirect(url, code=301)
