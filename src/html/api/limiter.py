# html/api/limiter.py

"""
Rate limiting for the planetGen API (Flask-Limiter).

One shared `Limiter` instance, created here uninitialized (no `app` bound
yet) and wired up in `app.py`'s `create_app` via `limiter.init_app(app)`
-- the standard Flask-extension pattern, needed because `routes.py`'s
individual view functions want to decorate themselves with
`@limiter.limit(...)` for a stricter per-route limit, and importing a
single shared instance is how both modules end up talking about the same
limiter rather than each creating their own.

The global default limits and storage backend are configured through
`app.config` (`RATELIMIT_DEFAULT`/`RATELIMIT_STORAGE_URI`, see
`config.py`) rather than passed to the constructor here, since `Limiter`
reads those from the app at `init_app` time -- this module doesn't need
its own copy of that configuration.
"""

from flask import current_app, request
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address

from stellarObjects.appconfig import DEFAULT_CONFIG

IN_PROCESS_ENVIRON_KEY = "planetgen.in_process"
"""str: WSGI environ key `web/transport.py` sets on the API requests the
Flask-served pages (`html/web/`) dispatch to this app in-process. Not an
`HTTP_*` key, so no client can set it from outside."""


def is_in_process_call():
    """
    True for an API request made in-process by one of the Flask-served
    pages (`web/transport.py`) rather than by a real client.

    Such a call is exempt from the *default* limits only: a page view
    already fetches its data this way (the CGI pages did the same over
    loopback HTTP), and counting each of those against the viewer's
    hourly budget would throttle ordinary browsing. A route's own
    explicit `@limiter.limit(...)` (login, every write) still applies,
    keyed by the real viewer's address, which `web/transport.py` passes
    through as `REMOTE_ADDR`.
    """
    return bool(request.environ.get(IN_PROCESS_ENVIRON_KEY))


limiter = Limiter(key_func=get_remote_address, default_limits_exempt_when=is_in_process_call)


# ---------------------------------------------------------------------
# Per-page limits (the HTML pages and /api/health)
# ---------------------------------------------------------------------
#
# Each page view does real work (API calls, database queries, template
# rendering) on one of a fixed number of WSGI threads, so without a limit
# one client could tie up every thread. The expensive public pages get
# their own per-IP limit and every other page shares one generous limit
# (`web/__init__.py` applies `page_limit("other")` to the pages'
# blueprint). The values come from `app.config["RATELIMIT_PAGES"]`
# (`config.json`'s `ratelimit.pages`, see `config.py` and
# `docs/config.md`), read on each request:
#
# - `search`: `/search` (default 30 per minute)
# - `galaxy`: `/galaxy`, the Galaxy Map page (default 60 per minute)
# - `galaxy_tiles`: `/galaxy/tiles` and `/galaxy/stage`, fetched as the map's camera moves
#   (default 600 per minute)
# - `health`: `/api/health` (default 60 per minute)
# - `other`: every other page, all together (default 300 per minute)
#
# An empty value turns that limit off. Exceeding one is a 429: the HTML
# error page for a page, JSON under `/api` and for `/galaxy/tiles`
# (`app.py`'s 429 handler).

DEFAULT_PAGE_LIMITS = DEFAULT_CONFIG["ratelimit"]["pages"]
"""dict: The built-in values, used for any name missing from the app's
`RATELIMIT_PAGES` (a test config that doesn't set it at all)."""


PAGE_LIMITS_OFF = {name: "" for name in DEFAULT_PAGE_LIMITS}
"""dict: A `RATELIMIT_PAGES` that turns every page limit off (for a test
harness that loads many pages from one address)."""


def configured_page_limit(name):
    """The limit string for `name`, or `""` when it is turned off."""
    limits = current_app.config.get("RATELIMIT_PAGES")
    if not isinstance(limits, dict):
        limits = {}
    value = limits.get(name, DEFAULT_PAGE_LIMITS.get(name, ""))
    return value.strip() if isinstance(value, str) else ""


def page_limit(name):
    """
    A `limiter.limit(...)` decorator for the page limit `name`, applied
    to a view function or to a blueprint (every route in it without its
    own limit then shares the one counter). Replaces the app-wide default
    limits for what it covers.
    """
    return limiter.shared_limit(
        # Not consulted while the limit is off (exempt_when), but always a
        # valid limit string.
        lambda: configured_page_limit(name) or "1000 per second",
        scope=f"page:{name}",
        exempt_when=lambda: not configured_page_limit(name),
    )
