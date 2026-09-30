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

from flask import Blueprint, request

_HTML_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_LIB_DIR = os.path.join(_HTML_DIR, "lib")
_STATIC_DIR = os.path.join(_HTML_DIR, "static")
# lib/ holds the modules the pages share (apiclient, fmt, pagination, the
# map renderers).
if _LIB_DIR not in sys.path:
    sys.path.insert(0, _LIB_DIR)

from fmt import STATIC_VERSION  # noqa: E402
from fmt import static_url as fmt_static_url  # noqa: E402
from stellarObjects.appconfig import load_config  # noqa: E402

from . import csrf, errors, transport  # noqa: E402
from .helpers import SECTIONS, current_admin, page_url  # noqa: E402

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

STATIC_DIR = _STATIC_DIR
"""str: `src/html/static/`. `create_app` makes it the app's static folder
so `/static/...` works under `python src/html/wsgi.py` and in tests; in
production Apache serves `/static/` itself (examples/apache/)."""

bp = Blueprint("web", __name__, template_folder="templates")


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
        "sections": SECTIONS,
        "static_url": static_url,
        "page_url": page_url,
        "current_admin": current_admin,
        "csrf_field": csrf.csrf_field,
    }


errors.register(bp)

from . import views  # noqa: E402,F401 -- registers the routes on bp
from . import system_pages  # noqa: E402,F401 -- /system, /phenomena, /phenomenon
from . import admin_pages  # noqa: E402,F401 -- /login, /logout, /account, /admin, /admin/stats
from . import galaxy_views  # noqa: E402,F401 -- /galaxy, /galaxy/tiles
from . import nav_page, sector_page  # noqa: E402,F401 -- /nav, /sector/<id>
from . import generate_page  # noqa: E402,F401 -- /admin/generate
from . import system_page  # noqa: E402,F401 -- /admin/generate/system
from . import old_urls  # noqa: E402,F401 -- /<name>.py -> 301 to the page that replaced it


def init_app(app, limiter=None):
    """
    Registers the pages on `app`. Called by `api.app.create_app`.

    Args:
        app (flask.Flask): The app.
        limiter (flask_limiter.Limiter, optional): When given, the pages
            themselves are exempt from rate limiting (the old CGI pages never
            are). The API calls they make in-process skip only the
            default limits; see `api.limiter.is_in_process_call`.
    """
    app.register_blueprint(bp)
    app.before_request(csrf.protect)
    app.after_request(csrf.set_cookie)
    if limiter is not None:
        limiter.exempt(bp)
    transport.install()
    if app.config.get("SECRET_KEY_IS_EPHEMERAL"):
        app.logger.warning(
            "No secret_key in config.json (or PLANETGEN_SECRET_KEY): using a random key for this "
            "process. Forms still work, but ones left open across a restart must be resubmitted."
        )
