# planetgen/web/errors.py

"""
HTML error pages for the Flask-served pages, with these safety rules:

- an `apiclient.NotFoundError` (a bad id, no such database) is a 404,
- an `apiclient.ApiError` (the API failed) is a 502,
- anything else is a 500 whose page says only "An unexpected error
  occurred." The traceback goes to the log (and the debug log, when
  debug is on), never into the page.

`api/app.py`'s app-wide handlers call `render_error` for any request
outside `/api`, so an unknown URL or a 405 on a page gets this HTML page,
while `/api/...` keeps its JSON errors.
"""

import os
import time

from flask import current_app, request

from planetgen.web.lib import apiclient
from planetgen.util import log
from planetgen.util.logpaths import debug_enabled

from .helpers import render_page

_TITLES = {
    400: "Bad request",
    403: "Forbidden",
    404: "Not found",
    405: "Method not allowed",
    413: "Too large",
    429: "Too many requests",
    500: "Something went wrong",
    502: "Data unavailable",
    504: "Took too long",
}

TIMEOUT_MESSAGE = ("The database is busy (a generation run may be using it) and could not answer this page in "
                   "time, so it was stopped. The page asks again by itself in a few seconds; a narrower search "
                   "or a smaller page is quicker.")
"""str: What a page whose query hit the statement time limit (PERF.17) says."""

RETRY_AFTER_SECONDS = 10
"""int: How long a page that timed out waits before asking again (PERF.64), and what `Retry-After` says."""


def timeout_page():
    """The 504 page of a statement that hit the time limit: it reloads itself after `RETRY_AFTER_SECONDS`, so a
    page opened during a fill fills in once the database has caught up (PERF.64)."""
    html, status = render_page("error.html", title=_TITLES[504], message=TIMEOUT_MESSAGE, status=504,
                               retry_after=RETRY_AFTER_SECONDS)
    return html, status, {"Retry-After": str(RETRY_AFTER_SECONDS)}


def render_error(status, message):
    """
    Renders `error.html` with `status`. `message` is plain text (escaped
    by the template).
    """
    return render_page("error.html", title=_TITLES.get(status, "Error"), message=message, status=status)


def unexpected_error_message():
    """The one message a 500 page shows."""
    if debug_enabled():
        return (f"An unexpected error occurred. The traceback is in the debug log "
                f"(process {os.getpid()}, {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}).")
    return "An unexpected error occurred."


def register(bp):
    """Error handlers for exceptions raised inside `web` views."""

    @bp.app_errorhandler(apiclient.NotFoundError)
    def _not_found(exc):
        log.debug(f"Not found: {exc}")
        return render_error(404, str(exc) or "Not found.")

    @bp.app_errorhandler(apiclient.ApiError)
    def _api_error(exc):
        if exc.status_code == 504:
            log.debug(f"Query time limit on {request.path}: {exc}")
            return timeout_page()
        current_app.logger.error(f"API error while building {request.path}: {exc}")
        log.exception(f"API error while building the page: {exc}")
        return render_error(502, "The data for this page could not be loaded. Please try again shortly.")
