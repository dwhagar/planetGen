# html/web/__init__.py

"""
The HTML pages of the site, served by the same Flask app as the JSON API
(`api/app.py`'s `create_app` calls `init_app` below).

Every page has a real GET route here (`/`, `/sectors`, `/sector/<id>`,
...) and a Jinja2 template in `web/templates/` extending `base.html`.
The old CGI URLs (`/<name>.py`) 301 to their replacements
(`old_urls.py`). See `docs/html-interface.md` ("Flask pages") for the
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
import sys

from flask import Blueprint, current_app, has_app_context, request
from markupsafe import Markup

_HTML_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_LIB_DIR = os.path.join(_HTML_DIR, "lib")
_STATIC_DIR = os.path.join(_HTML_DIR, "static")
# lib/ holds the modules the pages share (apiclient, fmt, pagination, the
# map renderers).
if _LIB_DIR not in sys.path:
    sys.path.insert(0, _LIB_DIR)

import apiclient  # noqa: E402
import classref  # noqa: E402
import pagecache  # noqa: E402
from fmt import STATIC_VERSION, utc_time_html  # noqa: E402
from fmt import static_url as fmt_static_url  # noqa: E402
from api.limiter import page_limit  # noqa: E402
from stellarObjects.appconfig import load_config  # noqa: E402

from . import csrf, errors, transport  # noqa: E402
from .helpers import current_admin, page_url, visible_sections  # noqa: E402

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

SECURITY_HEADERS = (
    ("X-Content-Type-Options", "nosniff"),
    ("X-Frame-Options", "DENY"),
    ("Referrer-Policy", "no-referrer"),
    ("Content-Security-Policy", CONTENT_SECURITY_POLICY),
)
"""tuple: The one place the HTML pages' security headers are decided
(`api/app.py` adds them to every `text/html` response). The JSON API sets
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
# /galaxy/tiles have their own (see api/limiter.py's page limits).
page_limit("other")(bp)


def static_url(filename):
    """
    `/static/<filename>?v=<version>` -- `fmt.static_url` made absolute, since Flask pages live at nested paths
    (`/sector/5`). A release changes every static URL, so Apache can let
    browsers cache `static/` for a year.
    """
    return f"{request.script_root}/{fmt_static_url(filename)}"


@bp.app_context_processor
def _template_globals():
    return {
        "site_name": load_config()["site_name"],
        "site_version": STATIC_VERSION,
        "sections": visible_sections,
        "static_url": static_url,
        "page_url": page_url,
        "current_admin": current_admin,
        "csrf_field": csrf.csrf_field,
        "utc_time": lambda value: Markup(utc_time_html(value)),
    }


errors.register(bp)

from . import views  # noqa: E402,F401 -- registers the routes on bp
from . import system_pages  # noqa: E402,F401 -- /system, /phenomena, /phenomenon
from . import admin_pages  # noqa: E402,F401 -- /login, /logout, /account, /admin, /admin/stats
from . import galaxy_views  # noqa: E402,F401 -- /galaxy, /galaxy/tiles
from . import nav_page, sector_page  # noqa: E402,F401 -- /nav, /sector/<id>
from . import generate_page  # noqa: E402,F401 -- /admin/generate
from . import system_page  # noqa: E402,F401 -- /admin/generate/system
from . import population_pages  # noqa: E402,F401 -- /species, /polities
from . import class_pages  # noqa: E402,F401 -- /classes, /classes/<type>, /classes/<type>/<code>
from . import old_urls  # noqa: E402,F401 -- /<name>.py -> 301 to the page that replaced it


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


def _install_page_cache(app):
    """
    Gives `app` its own `pagecache.ResponseCache` (TODO 8) unless
    `page_cache.enabled` is off, and clears it after every successful
    write under `/api` -- a page's admin form or an API client alike.
    """
    settings = pagecache.settings_from(load_config())
    if not settings["enabled"]:
        return
    cache = pagecache.ResponseCache(lambda db: apiclient.get_galaxy_changes(db)["stamp"], settings)
    app.extensions[PAGE_CACHE_EXTENSION] = cache
    apiclient.set_response_cache(_app_page_cache)

    @app.after_request
    def _clear_page_cache_after_writes(response):
        if (request.method in _WRITE_METHODS and request.path.startswith("/api/")
                and response.status_code < 400):
            cache.clear()
        return response


def init_app(app, limiter=None):
    """
    Registers the pages on `app`. Called by `api.app.create_app`.

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
    # The class reference catalog is built from the generator's tables
    # now, once, so the /classes pages never build it on a request.
    classref.catalog()
    app.before_request(csrf.protect)
    app.after_request(csrf.set_cookie)
    transport.install()
    _install_page_cache(app)
    if app.config.get("SECRET_KEY_IS_EPHEMERAL"):
        app.logger.warning(
            "No secret_key in config.json (or PLANETGEN_SECRET_KEY): using a random key for this "
            "process. Forms still work, but ones left open across a restart must be resubmitted."
        )
