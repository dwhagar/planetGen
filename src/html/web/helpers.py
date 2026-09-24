# html/web/helpers.py

"""
The helpers every Flask-served page uses. A page view normally needs only
these:

    from web.helpers import crumb, db_name, page_url, render_page, trusted_html

    @bp.route("/sector/<int:sector_id>")
    def sector(sector_id):
        detail = apiclient.get_sector(db_name(), sector_id)
        return render_page(
            "sector.html",
            title=detail["name"],
            section="sectors",
            breadcrumbs=[crumb("Sectors", "sectors"), crumb(detail["name"])],
            detail=detail,
            map_html=trusted_html(starmap.render_map_panel(...)),
        )

- `db_name()`: the one database this site shows, from config. Never read
  it from the request and never put it in a URL.
- `page_url(name, **params)`: the URL of another page by its endpoint
  name. Moved pages resolve through `url_for("web.<name>")`; pages still
  served by CGI resolve through `LEGACY_PAGES` below. A page PR that moves
  a page adds its route under the same endpoint name and deletes its
  `LEGACY_PAGES` entry; every link to it then flips on its own.
- `crumb(label, name=None, **params)`: one breadcrumb (the last one, the
  current page, takes no `name`).
- `render_page(template, *, title, section=None, breadcrumbs=(), **ctx)`:
  renders a template that extends `base.html`.
- `trusted_html(html)`: marks HTML built by a `lib/` renderer (pager, star
  map, system map, ...) safe to insert unescaped. Those renderers escape
  their own inputs with `fmt.esc`; never pass request or database text
  through this directly.
- `current_admin()`: the logged-in admin (or `None`), looked up once per
  request.
- `pager(...)`: `lib/pagination.render_pagination` in GET mode, as Markup.
"""

from urllib.parse import urlencode

from flask import current_app, g, render_template, request, url_for
from markupsafe import Markup

import apiclient
from pagination import render_pagination

# ---------------------------------------------------------------------
# Pages still served by CGI
# ---------------------------------------------------------------------

LEGACY_PAGES = {
    # endpoint name: (CGI script, {page_url keyword: CGI parameter})
    "galaxy": ("galaxy.py", {"quadrant": "quadrant", "page": "page"}),
    "sector": ("sector.py", {"sector_id": "id", "contents_page": "contents_page"}),
    "system": ("system.py", {"system_id": "id"}),
    "phenomena": ("phenomena.py", {"page": "page"}),
    "phenomenon": ("phenomenon.py", {"phenomenon_type": "type", "phenomenon_id": "id"}),
    "nav": ("nav.py", {"from_id": "from", "to_id": "to"}),
    "login": ("login.py", {}),
    "logout": ("logout.py", {}),
    "account": ("changecreds.py", {}),
    "admin": ("admin.py", {}),
    "admin_stats": ("adminstats.py", {"names_page": "names_page"}),
}
"""dict: Every page the Flask app links to that is still a CGI script.
Links to them are plain GET links (`/<script>?db=<db>&...`): the CGI
pages read their parameters with `page.nav_params()`, which falls back to
the query string for a GET request. They are the one place the database
name still appears in a URL, until each page moves. The page PRs delete
entries here as they add the matching `web.<name>` route."""

def db_name():
    """
    The one database this site shows: `WEB_DATABASE` if a config sets
    it, else the API's own configured database (`config.json`'s
    `mysql.database`, or `PLANETGEN_MYSQL_DATABASE`).
    """
    return current_app.config.get("WEB_DATABASE") or current_app.config["MYSQL_CONFIG"].database


def page_url(name, _anchor=None, **params):
    """
    URL of another page, by endpoint name (`"index"`, `"sectors"`,
    `"sector"`, ...), whether it has moved to Flask yet or not.

    Args:
        name (str): The endpoint name without the `web.` prefix.
        _anchor (str, optional): A fragment to append (`#...`).
        **params: The route's own parameters (`sector_id=5`) plus any
            query parameters. `None` values are dropped.

    Raises:
        KeyError: For a name that is neither a route nor in
            `LEGACY_PAGES` -- a typo should fail loudly in tests.
    """
    params = {key: value for key, value in params.items() if value is not None}
    endpoint = f"web.{name}"
    if endpoint in current_app.view_functions:
        return url_for(endpoint, _anchor=_anchor, **params)
    script, translation = LEGACY_PAGES[name]
    query = [("db", db_name())]
    for key, value in params.items():
        query.append((translation.get(key, key), value))
    url = f"/{script}?{urlencode(query)}"
    return f"{url}#{_anchor}" if _anchor else url


def crumb(label, name=None, _anchor=None, **params):
    """One breadcrumb: `{"label", "url"}`. Leave `name` out for the
    current page (the last crumb), which is shown as text."""
    return {"label": label, "url": page_url(name, _anchor=_anchor, **params) if name else None}


def trusted_html(html):
    """
    Marks HTML from a `lib/` renderer safe for a template (`{{ x }}`
    would otherwise escape it). Only for HTML those renderers built
    themselves -- they escape every value they interpolate.
    """
    return Markup(html or "")


def pager(page_param, page, total, anchor=None, label="Pages", keep=None):
    """
    The shared pager (`lib/pagination.py`) for one table on the current
    page, with plain `?<page_param>=N` GET links back to the current
    URL.

    Args:
        page_param (str): e.g. `"sectors_page"`.
        page (int): Current (clamped) page number.
        total (int): Total rows.
        anchor (str, optional): Element id to land on.
        label (str): Accessible name of the pager's `<nav>`.
        keep (dict, optional): Other query parameters to carry through
            unchanged (e.g. another table's page number).
    """
    keep = {key: value for key, value in (keep or {}).items() if value not in (None, "", 1)}
    return Markup(render_pagination(
        request.path, keep, page_param, page, total, anchor=anchor, label=label, method="get",
    ))


_UNSET = object()


def current_admin():
    """
    The logged-in admin (`{"username", "must_change_credentials"}`) or
    `None`, from the request's session cookie. Looked up once per request
    (in-process `GET /api/auth/me`) and cached on `g`; a lookup failure
    counts as "not logged in" rather than breaking the page.
    """
    admin = g.get("web_admin", _UNSET)
    if admin is _UNSET:
        try:
            admin = apiclient.auth_me(request.headers.get("Cookie"))
        except (apiclient.ApiError, apiclient.NotFoundError):
            admin = None
        g.web_admin = admin
    return admin


SECTIONS = (
    # (endpoint name, label)
    ("galaxy", "Galaxy"),
    ("sectors", "Sectors"),
    ("systems", "Systems"),
    ("phenomena", "Phenomena"),
    ("nav", "Nav"),
)
"""tuple: The header's main sections, in order. A view marks one active
with `render_page(section=...)`."""


def render_page(template, *, title, section=None, breadcrumbs=(), description=None, status=200, **context):
    """
    Renders `template` (which extends `base.html`).

    Args:
        template (str): Template file under `web/templates/`.
        title (str): Page title: the `<h1>` and the `<title>` prefix.
        section (str, optional): The active `SECTIONS` name, marked with
            `aria-current="page"` in the header.
        breadcrumbs (list[dict]): `crumb(...)` items, outermost first.
            The home crumb is added in front automatically.
        description (str, optional): `<meta name="description">`.
        status (int): HTTP status.
        **context: Anything else the template needs.

    Returns:
        tuple: `(html, status)`, ready to return from a view.
    """
    return render_template(
        template,
        title=title,
        section=section,
        breadcrumbs=list(breadcrumbs),
        description=description,
        **context,
    ), status
