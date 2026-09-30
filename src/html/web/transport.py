# html/web/transport.py

"""
In-process transport for `lib/apiclient.py`.

The Flask-served pages call the same `apiclient` functions the CGI pages
use (`get_sectors(db, ...)`, `auth_me(cookie_header)`, ...). Inside the
Flask app those calls must not go out over HTTP to the very process
making them, so `install()` (called once from `web.init_app`) gives
`apiclient` this transport: each call is dispatched straight through the
app's own WSGI callable, so it runs the same route, validation, auth
check and JSON serialization an HTTP call would, without a socket.

Details that matter:

- Only used while a Flask request is active. Outside one (a CLI, a test
  calling `apiclient` directly, a CGI process) it declines and
  `apiclient` falls back to HTTP, so CGI pages are unaffected.
- Each call runs in its own application context, so the API's
  request-scoped database connections (`g.db`, `g.control_db`) are
  opened and closed per call exactly as they are for a real request, and
  never leak into the calling page's own `g`.
- The viewer's address is passed through as `REMOTE_ADDR`, so a route's
  own rate limit (login, writes) applies to the real visitor. The
  in-process marker (`api.limiter.IN_PROCESS_ENVIRON_KEY`) exempts the
  call from the app-wide *default* limits only -- see
  `api.limiter.is_in_process_call`.
- Cookies are forwarded only when the `apiclient` call passes a
  `cookie_header` (the `auth_*` functions), exactly as over HTTP.
"""

from flask import current_app, has_request_context, request
from werkzeug.test import EnvironBuilder, run_wsgi_app

import apiclient
from api.limiter import IN_PROCESS_ENVIRON_KEY

API_PREFIX = "/api"


def in_process_transport(method, target, data, headers, timeout):
    """`apiclient.set_transport` callable -- see this module's docstring.
    `timeout` is ignored: nothing here can hang on a network read."""
    if not has_request_context():
        return None
    app = current_app._get_current_object()
    path, _, query = target.partition("?")
    builder = EnvironBuilder(
        path=API_PREFIX + path,
        method=method,
        query_string=query,
        data=data,
        headers=headers,
        environ_base={
            "REMOTE_ADDR": request.remote_addr or "127.0.0.1",
            IN_PROCESS_ENVIRON_KEY: True,
        },
    )
    try:
        environ = builder.get_environ()
    finally:
        builder.close()

    with app.app_context():
        app_iter, status, response_headers = run_wsgi_app(app.wsgi_app, environ, buffered=True)
        try:
            body = b"".join(app_iter).decode("utf-8", errors="replace")
        finally:
            close = getattr(app_iter, "close", None)
            if close is not None:
                close()

    code, _, reason = status.partition(" ")
    return apiclient.TransportResponse(
        int(code), body, response_headers.getlist("Set-Cookie"), reason,
    )


def install():
    """Makes `apiclient` use `in_process_transport` inside Flask requests."""
    apiclient.set_transport(in_process_transport)
