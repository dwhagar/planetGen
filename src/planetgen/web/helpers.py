# planetgen/web/helpers.py

"""
The helpers every Flask-served page uses. A page view normally needs only
these:

    from planetgen.web.helpers import crumb, db_name, page_url, render_page, trusted_html

    @bp.route("/sector/<uid:sector_id>")
    def sector(sector_id):
        detail = apiclient.get_sector(db_name(), sector_id)
        return render_page(
            "sector.html",
            title=detail["name"],
            section="sectors",
            breadcrumbs=[crumb("Sectors", "sectors"), crumb(detail["name"])],
            detail=detail,
            map_html=trusted_html(starmap.map_scene_data(...)),
        )

- `db_name()`: the one database this site shows, from config. Never read
  it from the request and never put it in a URL.
- `page_url(name, **params)`: the URL of another page by its endpoint
  name (`url_for("web.<name>")`).
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
- `pager(...)`: `lib/pagination.render_pagination`, as Markup.
- `bookmark(kind, value, name, url, sector_id=None)`: what a page's
  ☆ Bookmark button saves (`templates/partials/bookmark.html`).
"""

from flask import current_app, g, render_template, request, url_for
from markupsafe import Markup

from planetgen.web.lib import apiclient
from planetgen.api.authz import SESSION_COOKIE_NAME
from planetgen.web.lib.pagination import render_pagination

from . import csrf


def db_name():
    """
    The one database this site shows: `WEB_DATABASE` if a config sets
    it, else the API's own configured database (`config.json`'s
    `mysql.database`, or `PLANETGEN_MYSQL_DATABASE`).
    """
    return current_app.config.get("WEB_DATABASE") or current_app.config["MYSQL_CONFIG"].database


def bookmark(kind, value, name, url, sector_id=None):
    """
    The entry a page's ☆ Bookmark button saves in the browser
    (`static/bookmarks.js`, MAP.23; rendered by
    `templates/partials/bookmark.html`): `kind` is `"system"`, a
    phenomenon type or `"sector"`; `value` is the NAV endpoint
    (`nav_page.endpoint`) or the sector's designation; `url` is the page
    it opens; `sector_id` lets the NAV page open a sector's system picker.
    `db` names the browser's list, one per database.
    """
    return {"db": db_name(), "kind": kind, "value": value, "name": name, "url": url, "sector_id": sector_id}


def page_url(name, _anchor=None, **params):
    """
    URL of another page, by endpoint name (`"index"`, `"sectors"`,
    `"sector"`, ...).

    Args:
        name (str): The endpoint name without the `web.` prefix.
        _anchor (str, optional): A fragment to append (`#...`).
        **params: The route's own parameters (`sector_id=5`) plus any
            query parameters. `None` values are dropped.

    Raises:
        werkzeug.routing.BuildError: For a name that is not a route -- a
            typo should fail loudly in tests.
    """
    params = {key: value for key, value in params.items() if value is not None}
    return url_for(f"web.{name}", _anchor=_anchor, **params)


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
    The shared pager (`planetgen/web/lib/pagination.py`) for one table on the current
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
        request.path, keep, page_param, page, total, anchor=anchor, label=label,
    ))


_UNSET = object()


def generate_target(admin):
    """Where the map pages' Generate buttons post (`static/
    generatebuttons.js`), for an admin who can use the Generate page;
    `None` for everyone else."""
    if admin is None or admin.get("must_change_credentials"):
        return None
    return {"url": url_for("web.generate"), "csrfField": csrf.FIELD_NAME, "csrfToken": csrf.csrf_token()}


def current_admin():
    """
    The logged-in admin (`{"username", "must_change_credentials"}`) or
    `None`, from the request's session cookie. Looked up once per request
    (in-process `GET /api/auth/me`) and cached on `g`; a lookup failure
    counts as "not logged in" rather than breaking the page. A request with
    no session cookie at all (nearly every visitor) skips the lookup.
    """
    admin = g.get("web_admin", _UNSET)
    if admin is _UNSET:
        if not request.cookies.get(SESSION_COOKIE_NAME):
            g.web_admin = None
            return None
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
    ("nearby", "Nearby"),
    ("classes", "Classes"),
    ("species", "Species"),
)
"""tuple: The header's main sections, in order. A view marks one active
with `render_page(section=...)`."""

POPULATION_SECTIONS = {"species"}
"""set: Sections shown only once population data exists
(`population_status()["species"]`)."""


def population_status():
    """
    What population data this site's database has (`apiclient.
    get_population_status`: `generated`/`species`/`polities`/
    `territories`), looked up once per request. Until a population pass
    has made species, the Species pages, "Dominant species" and
    "Territory of ..." stay hidden.
    """
    status = g.get("web_population")
    if status is None:
        try:
            status = apiclient.get_population_status(db_name())
        except Exception:  # noqa: BLE001 -- an unreachable API hides the pages, never breaks one
            status = dict(apiclient.POPULATION_NONE)
        g.web_population = status
    return status


def visible_sections():
    """`SECTIONS` less the ones with nothing to show yet."""
    hidden = set() if population_status()["species"] else POPULATION_SECTIONS
    return tuple((name, text) for name, text in SECTIONS if name not in hidden)


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
