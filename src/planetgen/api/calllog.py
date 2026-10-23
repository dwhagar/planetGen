# planetgen/api/calllog.py

"""
Every API call in the activity log (API.15)
===========================================

One `API call` line for each request to `/api/...`, written when the answer is ready (`docs/design/api-design-standards.md`
section 8): the time and client address of every activity-log line, then the `user` (the key's owner for an API key, the
signed-in user for the web site, `god` for the console, `-` when nobody), `via` (`key`, `session`, `anonymous` or `console`),
the HTTP `method`, the matched `route` (`-` when none matched) and the raw `path` (cut to 200 characters), the HTTP `status`
code, the time taken in milliseconds (`ms`), the bytes in and out, the request id (`rid`, also sent back in the
`X-Request-Id` header), the user agent (120 characters) and the API version the client said it speaks (`client_api`, from
`X-PlanetGen-API-Version`, API.31). An API key adds `key_id` and `key_prefix`, never the key. `internal=1` marks a call a web
page made to the API in-process (`web/transport.py`), about six for each page view; the setting `api_log.internal` set to `off`
leaves those out. `ratelimited=1` marks a 429. A body, a token or an `Authorization` header is never written.

The status codes are HTTP's, for parity with the web server's own logs. A request the web server refuses before this program sees
it (a 413 from `LimitRequestBody`, a failed TLS handshake) appears only in the web server's log.

The console: `activity_log.console_call` writes the same line for an API-shaped action a command-line tool runs with no request, with
`user=god`, `via=console` and a status of 200 for success or 500 for failure.

Writing a line never stops a request (`activity_log.event` never raises, and this module catches the rest).
"""

import re
import time
import uuid

from flask import g, request

from planetgen.admin import activity_log
from planetgen.api.limiter import is_in_process_call
from planetgen.api.version import API_VERSION_HEADER
from planetgen.util import log
from planetgen.util.settings import get_settings

REQUEST_ID_HEADER = "X-Request-Id"
"""str: Carried by every `/api` answer; a client's own (a short token) is kept."""

AGENT_LENGTH = 120
PATH_LENGTH = 200

_CLIENT_REQUEST_ID = re.compile(r"^[A-Za-z0-9._-]{8,64}$")
_BEARER = "Bearer "


def _is_api_path():
    return request.path == "/api" or request.path.startswith("/api/")


def _begin():
    g.api_call_started = time.perf_counter()
    if not _is_api_path():
        return
    offered = request.headers.get(REQUEST_ID_HEADER, "")
    g.api_request_id = offered if _CLIENT_REQUEST_ID.match(offered) else uuid.uuid4().hex


def _caller():
    """`(user, via, key_id, key_prefix)` of the request: from the route's own check when it ran one, else from the
    credential the request carries, looked up once (a refused or absent credential is anonymous)."""
    from planetgen.admin import auth as adminAuth
    from planetgen.api import authz
    from planetgen.api.common import get_control_db

    admin = g.get("admin_user")
    bearer = request.headers.get("Authorization", "").startswith(_BEARER)
    if admin is None and bearer:
        admin = authz.api_key_admin(request.headers["Authorization"][len(_BEARER):].strip())
    elif admin is None and request.cookies.get(authz.SESSION_COOKIE_NAME):
        admin = adminAuth.validate_session(get_control_db(), request.cookies.get(authz.SESSION_COOKIE_NAME))
    if admin is None:
        return None, "anonymous", None, None
    if bearer:
        return admin["username"], "key", admin.get("api_key_id"), admin.get("api_key_prefix")
    return admin["username"], "session", None, None


def _safe_caller():
    try:
        return _caller()
    except Exception as exc:  # noqa: BLE001 -- the log never decides the answer (a control database that is down)
        log.debug(f"API call log: could not tell who made {request.path}: {exc}")
        return None, "anonymous", None, None


def _end(response):
    if not _is_api_path():
        return response
    request_id = g.get("api_request_id")
    if request_id:
        response.headers[REQUEST_ID_HEADER] = request_id
    try:
        internal = is_in_process_call()
        if internal and get_settings().api_log.internal == "off":
            return response
        user, via, key_id, key_prefix = _safe_caller()
        started = g.get("api_call_started")
        rule = request.url_rule.rule if request.url_rule is not None else "-"
        activity_log.event(
            "API", "call", user=user, method=request.method, route=rule, path=request.path[:PATH_LENGTH],
            status=response.status_code, ms=(time.perf_counter() - started) * 1000 if started is not None else None,
            bytes_in=request.content_length or 0, bytes_out=response.calculate_content_length() or 0,
            via=via, key_id=key_id, key_prefix=key_prefix, internal=1 if internal else None,
            ratelimited=1 if response.status_code == 429 else None, rid=request_id,
            client_api=request.headers.get(API_VERSION_HEADER, "")[:8] or None,
            ua=request.headers.get("User-Agent", "")[:AGENT_LENGTH] or None)
    except Exception as exc:  # noqa: BLE001
        log.debug(f"API call log: line not written for {request.path}: {exc}")
    return response


def install(app):
    """Registers the request hooks: first among the `before_request` ones, so every answer has a start time."""
    app.before_request(_begin)
    app.after_request(_end)
