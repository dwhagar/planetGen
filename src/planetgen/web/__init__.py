# planetgen/web/__init__.py

"""
The HTML pages of the site, served by the same Flask app as the JSON API
(`web/app.py`'s `create_app` calls `init_app` below).

Every page has a real GET route here (`/`, `/sectors`, `/sector/<id>`,
...) and a Jinja2 template in `web/templates/` extending `base.html`.
See `docs/reference/html-interface.md` ("Flask pages") for the
full guide; in short, a new page is:

    # web/views.py
    @bp.route("/phenomena")
    def phenomena():
        data = apiclient.get_phenomena(db_name(), limit=..., offset=...)
        return render_page("phenomena.html", title="Phenomena", section="phenomena",
                           breadcrumbs=[crumb("Phenomena")], rows=data["items"])

    {# web/templates/phenomena.html #}
    {% extends "base.html" %}
    {% block content %} ... {% endblock %}

Pieces:

- `helpers.py`: `db_name`, `page_url`, `crumb`, `render_page`,
  `trusted_html`, `pager`, `current_admin`.
- `transport.py`: lets `apiclient` run in-process here (no HTTP hop).
- `csrf.py`: `csrf_field()` for POST forms; checked on every unsafe
  request outside `/api`.
- `errors.py`: HTML 404/502/500 pages that never show a traceback.
"""

import os
import re

from flask import Blueprint, current_app, g, has_app_context, request
from markupsafe import Markup, escape

# src/html/static/, which Apache serves directly.
_STATIC_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                           "html", "static")

from planetgen.admin import maintenance  # noqa: E402
from planetgen.queue import stopall  # noqa: E402
from planetgen.web.lib import apiclient  # noqa: E402
from planetgen.web.lib import classref  # noqa: E402
from planetgen.web.lib import habitability_symbols, poi_icons  # noqa: E402
from planetgen.web.lib import pagecache  # noqa: E402
from planetgen.web.lib import pagetags  # noqa: E402
from planetgen.web.lib.slowpages import take_notice as slow_notice  # noqa: E402
from planetgen.db import dataversion  # noqa: E402
from planetgen.web.lib.fmt import STATIC_VERSION, format_number, utc_time_html  # noqa: E402
from planetgen.web.lib.fmt import static_url as fmt_static_url  # noqa: E402
from planetgen.web.lib.modulegraph import module_graph  # noqa: E402
from planetgen.api.limiter import page_limit  # noqa: E402
from planetgen.api.version import API_VERSION, API_VERSION_HEADER  # noqa: E402
from planetgen.util.settings import get_settings  # noqa: E402
from planetgen.util import chem, settings_view  # noqa: E402

from . import adminregistry, csrf, errors, seo, transport  # noqa: E402
from .helpers import HELP_MENU, HELP_SECTIONS, NAV_MENU, current_admin, page_url, reference_link, visible_sections  # noqa: E402

# `default-src 'self'` covers scripts, styles, images, fonts and fetch():
# the pages only load same-origin `static/` files, the map pages build
# their canvases/WebGL textures in memory (no data:/blob: URLs), and
# nothing inline -- no `<script>` blocks, no `style="..."` attributes (JS
# setting `element.style` is not affected by CSP). The rest are the
# directives `default-src` does not fall back to: no `<base>` hijacking,
# forms only post back to this site, no framing (the modern form of
# X-Frame-Options, kept below for old browsers), and no plugins.
CONTENT_SECURITY_POLICY = ("default-src 'self'; base-uri 'self'; form-action 'self'; "
                           "frame-ancestors 'none'; object-src 'none'")

JOB_LOG_CONTENT_SECURITY_POLICY = CONTENT_SECURITY_POLICY + "; style-src 'self' 'unsafe-inline'"
"""str: The policy of the admin pages that show a job's output (`/admin/generate` and
`/admin/generate/jobs/<id>`): the same, plus inline styles, which Xterm.js
needs (it writes `<style>` elements for its colours and cell sizes). Scripts
stay `'self'` only, and every other page keeps `CONTENT_SECURITY_POLICY`."""

SECURITY_HEADERS = (
    ("X-Content-Type-Options", "nosniff"),
    ("X-Frame-Options", "DENY"),
    ("Referrer-Policy", "no-referrer"),
    ("Content-Security-Policy", CONTENT_SECURITY_POLICY),
)
"""tuple: The one place the HTML pages' security headers are decided
(`web/app.py` adds them to every `text/html` response). The JSON API sets
its own stricter set there."""

STRICT_TRANSPORT_SECURITY = ("Strict-Transport-Security", "max-age=31536000")
"""tuple: Added to every response (pages and API alike) to a request that
came in over HTTPS (`request.is_secure`), so a browser that has seen the
site once never falls back to plain HTTP for a year. Never on a plain-HTTP
response, where browsers ignore it anyway and a local `python
src/html/wsgi.py` must keep working. No `includeSubDomains`: the site
can't speak for other hosts under its domain."""

STATIC_DIR = _STATIC_DIR
"""str: `src/html/static/`. `create_app` makes it the app's static folder
so `/static/...` works under `python src/html/wsgi.py` and in tests; in
production Apache serves `/static/` itself (examples/apache/)."""

bp = Blueprint("web", __name__, template_folder="templates")
# Every page shares one generous per-IP limit; /search, /galaxy and
# /galaxy/nodes have their own (see api/limiter.py's page limits).
page_limit("other")(bp)


def static_url(filename):
    """
    `/static/<filename>?v=<version>` -- `fmt.static_url` made absolute, since Flask pages live at nested paths
    (`/sector/5`). A release changes every static URL, so Apache can let
    browsers cache `static/` for a year.
    """
    return f"{request.script_root}/{fmt_static_url(filename)}"


ICON_NAMES = frozenset({"edit", "delete", "show-on-map", "filter", "navigate", "menu", "settings",
                        "back", "forward", "up", "reset", "bookmark"})
"""frozenset[str]: The symbols in `static/icons.svg` (UX.28)."""


def icon(name):
    """
    One icon from the site's sprite (`static/icons.svg`, UX.28), as inline
    `<svg>` markup hidden from screen readers: the button or link around it
    carries the words, as an `aria-label` and a `title` tooltip.

    Raises:
        ValueError: `name` isn't in the sprite.
    """
    if name not in ICON_NAMES:
        raise ValueError(f"no icon named {name!r} in static/icons.svg")
    return Markup(f'<svg class="icon" aria-hidden="true" focusable="false">'
                  f'<use href="{escape(static_url("icons.svg"))}#{name}"></use></svg>')


def poi_icon(icon_id, tier=None):
    """One point-of-interest symbol (MAP.320) as markup, `icon_id` one of `poi_icons.ICONS`; `tier` (0 Blue to 3 Red)
    colours the habitable and the inhabited world.

    Raises:
        ValueError: `icon_id` or `tier` is not one the symbols have.
    """
    return Markup(poi_icons.symbol_html(icon_id, static_url, tier))


def _habitability_help_url(link):
    """The habitability page's URL for a symbol that links there (UX.181), `None` for one that does not (`link`
    false: the page's own legend)."""
    return page_url("habitability_levels") if link else None


def habitability_symbol(name, tier, detail=None, link=True):
    """
    One PHI-6 symbol (UX.101) as markup: `name` is a pressure, temperature, air, water, radiation, suitability or
    overall symbol (`habitability_symbols.SYMBOL_NAMES`), `tier` the stored rating (0 Blue to 3 Red) that colours it
    and names it in its words, and `detail` the figure that follows the rating in them ("Pressure: Blue, 101 kPa").
    The symbol links to its section of the habitability page unless `link` is false.

    Raises:
        ValueError: `name` or `tier` is not one the symbols have.
    """
    return Markup(habitability_symbols.symbol_html(name, tier, static_url("icons.svg"), detail, _habitability_help_url(link)))


def habitability_symbol_row(body, with_overall=True, with_life=True, page=None):
    """The symbols of one body's row (`habitability_symbols.symbol_row`), overall first and the life cluster last,
    each linked to its section of the habitability page (UX.181), or with `page` (the body's own page) to that page's
    Habitability and Life groups; empty for an unscored body."""
    return Markup(habitability_symbols.symbol_row(body, static_url("icons.svg"), with_overall, with_life,
                                                  _habitability_help_url(True), page))


def habitability_life(name, score, link=True):
    """One life score (UX.101) as markup: `name` is microbial, flora, fauna or sentient, `score` a whole number 0 to 9.
    It links to its row of the habitability page unless `link` is false.

    Raises:
        ValueError: `name` or `score` is not one the life symbols have.
    """
    return Markup(habitability_symbols.life_html(name, score, static_url("icons.svg"), _habitability_help_url(link)))


def math_check_failures():
    """The names of the math checks that failed when this process started
    (`init_app`), shown to admins on every page; empty when all passed."""
    from planetgen.physics import mathcheck
    return [result.name for result in mathcheck.startup_failures()]


def galaxy_status_now():
    """What `lib/galaxystatus` says about a run changing the galaxy, for the banner on every page (UX.3)."""
    from .helpers import db_name
    from .lib import galaxystatus
    return galaxystatus.status(db_name())


def edit_flashes():
    """What the last admin Delete or Regenerate flashed (`edit_actions`)."""
    from .edit_actions import messages
    return messages()


@bp.app_context_processor
def _template_globals():
    return {
        "site_name": get_settings().site_name,
        "api_version": API_VERSION,
        "api_version_header": API_VERSION_HEADER,
        "site_version": STATIC_VERSION,
        "seo_head": seo.head,
        "site_tagline": get_settings().site.tagline,
        "privacy_note": get_settings().seo.privacy_note,
        "sections": visible_sections,
        "nav_menu": NAV_MENU,
        "help_menu": HELP_MENU,
        "help_sections": HELP_SECTIONS,
        "static_url": static_url,
        "module_graph": module_graph,
        "icon": icon,
        "poi_icon": poi_icon,
        "poi_icons": poi_icons.ICONS,
        "habitability_symbol": habitability_symbol,
        "habitability_symbol_row": habitability_symbol_row,
        "habitability_life": habitability_life,
        "page_url": page_url,
        "reference_link": reference_link,
        "current_admin": current_admin,
        "admin_shell": lambda section: adminregistry.shell(section, page_url, current_admin()),
        "admin_job": lambda: generate_page.strip_view(),
        "maintenance_state": maintenance.read,
        "stop_state": stopall.state,
        "galaxy_status": galaxy_status_now,
        "pending_restart_options": settings_view.pending_restart,
        "csrf_field": csrf.csrf_field,
        "csrf_token": csrf.csrf_token,
        "edit_flashes": edit_flashes,
        "slow_notice": slow_notice,
        "math_check_failures": math_check_failures,
        "utc_time": lambda value: Markup(utc_time_html(value)),
        # UX.20: every count or measurement a template shows, `num(x)` or
        # `num(x, ",.2f")` -- scientific past 4 whole digits.
        "num": format_number,
    }


@bp.app_template_filter("formula")
def _formula(value):
    """UX.149: one chemical formula as HTML, `H<sub>2</sub>O` (`planetgen.util.chem.html`)."""
    return Markup(chem.html(value))


@bp.app_template_filter("chem")
def _chem(value):
    """UX.149: a text that names chemicals (a nebula's contents, a limit's meaning) as HTML, each formula
    `chem.split` finds there written out; text already marked up is left as it is."""
    if not value or isinstance(value, Markup):
        return value
    return Markup(chem.mark_html(value))


errors.register(bp)

from . import views  # noqa: E402,F401 -- registers the routes on bp
from . import system_pages  # noqa: E402,F401 -- /system, /phenomena, /phenomenon
from . import admin_pages  # noqa: E402,F401 -- /login, /logout, /account, /admin, /admin/stats
from . import galaxy_views  # noqa: E402,F401 -- /galaxy, /galaxy/nodes
from . import nav_page, nearby_page, sector_page  # noqa: E402,F401 -- /nav, /nearby, /sector/<id>
from . import generate_page  # noqa: E402,F401 -- /admin/generate
from . import queue_page  # noqa: E402,F401 -- /admin/queue
from . import admin_search_page  # noqa: E402,F401 -- /admin/search, /admin/palette.json
from . import settings_page, backups_page, images_page  # noqa: E402,F401 -- /admin/settings, /backups, /images
from . import seo_pages  # noqa: E402,F401 -- robots.txt, the sitemaps, security.txt, humans.txt, the manifest, /site-assets
from . import facility_type_pages  # noqa: E402,F401 -- /admin/facility-types
from . import system_page  # noqa: E402,F401 -- /admin/generate/system
from . import population_pages  # noqa: E402,F401 -- /species, /polities
from . import class_pages  # noqa: E402,F401 -- /classes, /classes/<type>, /classes/<type>/<code>
from . import references_pages  # noqa: E402,F401 -- /references (DOC.88), the one section that speaks out of universe
from . import maintenance_page  # noqa: E402,F401 -- /maintenance (UX.112)
from . import galaxy_status_page  # noqa: E402,F401 -- /galaxy-status (UX.3)
from . import object_lists  # noqa: E402,F401 -- /lists, /lists/<kind> (UX.141)
from . import preferences_page  # noqa: E402,F401 -- /preferences (UX.236)


PAGE_CACHE_EXTENSION = "planetgen_page_cache"
"""str: Where `_install_page_cache` keeps the app's
`pagecache.ResponseCache` (`app.extensions`)."""

_WRITE_METHODS = frozenset(("POST", "PUT", "PATCH", "DELETE"))


def _app_page_cache():
    """`apiclient.set_response_cache`'s provider: the current app's
    cache, or `None` outside an app (a CLI, a test calling `apiclient`
    directly)."""
    if not has_app_context():
        return None
    return current_app.extensions.get(PAGE_CACHE_EXTENSION)


def _written_database():
    """The database a write through the API changes: `?db=`, else the site's own."""
    return request.args.get("db") or current_app.config.get("WEB_DATABASE") or current_app.config["MYSQL_CONFIG"].database


def _sector_of_system(system_id):
    """The printed ID of the sector holding `system_id`, or `None` (standalone, or not there)."""
    try:
        return apiclient.get_system(_written_database(), system_id).get("sector_id")
    except Exception:  # noqa: BLE001 -- not found, or no database: the change is then placed more widely
        return None


def _record_write(response, sector_before):
    """Notes the data change a successful write through the API made (PERF.104), before the response goes out."""
    path = request.path[len("/api"):]
    kind, ids = pagetags.change_for_write(request.method, path)
    if kind in ("system-edited", "system-set"):
        sector = sector_before
        if sector is None and kind == "system-set" and request.method == "POST":
            sector = (request.get_json(silent=True) or {}).get("sector_id")
        if sector is not None:
            ids["sector_id"] = sector
        elif kind == "system-edited":
            kind, ids = "system-any", {}
    try:
        dataversion.record(_written_database(), kind, **ids)
    except ValueError:
        dataversion.record(_written_database(), "everything")


def _install_page_cache(app):
    """
    Gives `app` its own `pagecache.ResponseCache` (PERF.2) unless
    `page_cache.enabled` is off, and notes the data change of every successful
    write under `/api` -- a page's admin form or an API client alike -- by bumping the
    tags it touches (`pagetags.change_for_write`, `dataversion.record`), so only the answers that read it
    are rebuilt.
    """
    settings = get_settings().page_cache.model_dump()
    if not settings["enabled"]:
        return
    cache = pagecache.ResponseCache(apiclient.get_galaxy_changes, settings)
    app.extensions[PAGE_CACHE_EXTENSION] = cache
    apiclient.set_response_cache(_app_page_cache)

    @app.before_request
    def _remember_sector_of_system():
        g.page_cache_sector = None
        if request.method in _WRITE_METHODS and re.match(r"/api/systems/[^/]+", request.path):
            g.page_cache_sector = _sector_of_system(request.path.split("/")[3])

    @app.after_request
    def _record_api_writes(response):
        if (request.method in _WRITE_METHODS and request.path.startswith("/api/")
                and response.status_code < 400):
            _record_write(response, g.get("page_cache_sector"))
        return response


def init_app(app, limiter=None):
    """
    Registers the pages on `app`. Called by `planetgen.web.app.create_app`.

    Args:
        app (flask.Flask): The app.
        limiter (flask_limiter.Limiter, optional): Unused; kept for
            callers. The pages' per-IP limits are declared on `bp` and
            its routes (`api.limiter.page_limit`, configured by
            `RATELIMIT_PAGES`). The API calls they make in-process skip
            only the API's default limits; see
            `api.limiter.is_in_process_call`.
    """
    app.register_blueprint(bp)
    # The math check (TEST.63) runs once per process, now: a failure is
    # logged and shown to admins on every page (`math_check_failures`),
    # but the site still serves.
    math_check_failures()
    # The class reference catalog is built from the generator's tables
    # now, once, so the /classes pages never build it on a request (a dataset not downloaded yet does not stop the site).
    classref.preload()
    app.before_request(csrf.protect)
    app.after_request(csrf.set_cookie)
    transport.install()
    _install_page_cache(app)
    if app.config.get("SECRET_KEY_IS_EPHEMERAL"):
        app.logger.warning(
            "No secret_key in config.json (or PLANETGEN_SECRET_KEY): using a random key for this "
            "process. Forms still work, but ones left open across a restart must be resubmitted."
        )
