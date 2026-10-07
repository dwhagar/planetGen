# html/api/common.py

"""
Small pieces shared by `routes.py` (content-schema reads/writes) and
`auth.py`/`authz.py` (control-schema auth) -- split out here specifically
so neither of those two needs to import the other (both need `ApiError`
and a control-schema connection; only one direction would create an
import cycle otherwise).
"""

from urllib.parse import urlsplit

import pymysql
from flask import current_app, g, request

from planetgen.db.store import get_control_connection

CONTROL_DATABASE_UNAVAILABLE = "database unavailable"
"""str: What a route needing the control database says when it can't be
reached -- the same words as `routes.DATABASE_UNAVAILABLE`."""


class ApiError(Exception):
    """
    Raised by a route to end the request with a JSON `{"error": ...}` body
    and a specific status code, handled by `app.py`'s error handler. Beats
    each route hand-rolling its own `return jsonify(...), status` for
    validation failures, so every 4xx in this API is worded and shaped the
    same way.
    """

    def __init__(self, message, status_code=400):
        super().__init__(message)
        self.message = message
        self.status_code = status_code


def require_json_body():
    """
    Parses the request body as JSON, for every route that needs one.

    `request.get_json(silent=True)` returns `None` for a missing/empty
    body, a body that isn't valid JSON, *or* a `Content-Type` other than
    `application/json` -- this API doesn't need to tell those apart for a
    caller, they're all just "you didn't send a JSON object".

    Returns:
        dict: The parsed body.

    Raises:
        ApiError: If the body is missing, isn't valid JSON, or isn't a
            JSON *object* (a bare list/string/number is valid JSON but
            not a usable request body here).
    """
    body = request.get_json(silent=True)
    if not isinstance(body, dict):
        raise ApiError("request body must be a JSON object")
    return body


def get_control_db():
    """
    Returns the request-scoped control-schema connection, opening one on
    first use against `current_app.config["CONTROL_MYSQL_CONFIG"]`.
    Reused for the lifetime of the request, then closed by
    `close_control_db` in the app's teardown handler -- exactly the
    `routes.get_db`/`close_db` pattern, mirrored here for the control
    schema instead of a content database.

    If the server can't be reached, the request ends with a 503 saying
    only "database unavailable" (the server's address goes to the log),
    the same answer the content routes give, not a 500.
    """
    if "control_db" not in g:
        try:
            g.control_db = get_control_connection(current_app.config["CONTROL_MYSQL_CONFIG"], ensure_schema=False)
        except pymysql.err.OperationalError as exc:
            current_app.logger.error("Control database unavailable on %s %s: %s", request.method, request.path, exc)
            raise ApiError(CONTROL_DATABASE_UNAVAILABLE, status_code=503) from exc
    return g.control_db


def close_control_db(exception=None):
    db = g.pop("control_db", None)
    if db is not None:
        db.close()


def is_http_url(value):
    """
    True when `value` is an absolute `http://` or `https://` URL with a
    host -- the only kind of link a sector's `wiki_url` may hold, since
    the sector page renders it as an `<a href>` (a `javascript:` or
    `data:` URL there would run in the viewer's session).
    """
    if not isinstance(value, str) or value != value.strip():
        return False
    if any(ord(ch) < 0x21 or ord(ch) == 0x7f for ch in value):
        # Whitespace or control characters inside a URL, which browsers
        # strip or reinterpret (e.g. "java\tscript:").
        return False
    try:
        parts = urlsplit(value)
        port = parts.port  # noqa: F841 -- raises ValueError on a bad port
    except ValueError:
        return False
    return parts.scheme.lower() in ("http", "https") and bool(parts.hostname)
