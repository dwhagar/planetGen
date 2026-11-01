# planetgen/web/rebuild_gate.py

"""
The request side of `admin/cachebuild.py`: while the caches are being rebuilt, a request that the cache cannot answer
and that would build one itself gets a 503 with `Retry-After` (HTML page or JSON) instead of running into the
daemon's request limit, which would cut off every request in flight.

Two places ask: `check`, called by `apiclient` just before a page's in-process call goes to the API (so a page the
cache can answer never gets here), and `gate`, a `before_request` hook for a client calling a heavy `/api` route
directly. Not refused: the warm-up (it fills the caches) and work outside a request (a queue worker building a page job). A
signed-in admin is refused too: the daemon's request limit does not care who asked.

`gate` also marks the request a visitor's (`cachebuild.visitor`), so that code about to build something huge in the
process (`cachebuild.refuse_expensive_build`, the galaxy's cluster field) stops with the same 503 whichever route
reached it, in-process calls included.
"""

from flask import current_app, g, has_request_context, jsonify, request

from planetgen.admin import cachebuild
from planetgen.api.limiter import is_in_process_call


def _exempt():
    return not has_request_context() or bool(current_app.config.get("WARMUP"))


def check(path):
    """Raises `cachebuild.CacheRebuilding` for a visitor's heavy `path` while the caches are rebuilt."""
    if cachebuild.is_heavy(path) and cachebuild.state() is not None and not _exempt():
        raise cachebuild.CacheRebuilding(cachebuild.state())


def api_response(exc):
    """The JSON 503 for `exc`."""
    g.retry_after_fixed = True
    response = jsonify({"error": cachebuild.MESSAGE, "rebuilding": True, "cause": exc.state["cause"],
                        "since": exc.state.get("since")})
    response.status_code = 503
    response.headers["Retry-After"] = str(cachebuild.RETRY_AFTER_SECONDS)
    response.headers["Cache-Control"] = "no-store"
    return response


def gate():
    """The `before_request` hook: marks the request a visitor's, and refuses a direct heavy `/api` call."""
    cachebuild.visitor.set(not current_app.config.get("WARMUP"))
    if is_in_process_call() or not request.path.startswith("/api/") or request.method != "GET":
        return None
    try:
        check(request.path)
    except cachebuild.CacheRebuilding as exc:
        return api_response(exc)
    return None


def install(app):
    """Registers `gate` on `app` (after the maintenance gate)."""
    app.before_request(gate)
    app.teardown_request(lambda _exc: cachebuild.visitor.set(False))
