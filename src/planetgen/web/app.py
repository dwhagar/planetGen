# planetgen/api/app.py

"""
Flask application factory for the planetGen API.

See `docs/api.md` for how to run this in development and how it deploys
behind the project's existing Apache2 vhost (`examples/apache/`).
"""

import os
import time

import pymysql
from flask import Flask, abort, current_app, g, jsonify, make_response, request

from planetgen.db import store
from planetgen.admin import activity_log
from planetgen.util import log

from planetgen.api.admin import bp as admin_bp
from planetgen.api.auth import bp as auth_bp
from planetgen.api.common import ApiError, close_control_db
from planetgen.api.config import Config
from planetgen.api.edits import bp as edits_bp
from planetgen.api.limiter import limiter
from planetgen.api.population import bp as population_bp
from planetgen.api.routes import bp, close_db
from planetgen.api.workqueue import bp as workqueue_bp


def _is_api_request():
    """True for a request to the JSON API (`/api/...`); everything else
    is one of the HTML pages (`planetgen/web/`)."""
    return request.path == "/api" or request.path.startswith("/api/")


def create_app(config_object=Config):
    # The HTML pages (planetgen/web/), imported here rather than at the top so
    # `import planetgen.web.app` stays cheap for callers that never build an app.
    from planetgen import web

    # /static/ is src/html/static/ -- Apache serves it directly in
    # production (examples/apache/); this only matters for the dev server
    # (`python src/html/wsgi.py`) and tests.
    # Named outside the "planetgen" logger tree (which doesn't propagate):
    # the app's errors keep going to Flask's own stderr handler, so to
    # Apache's error log, and through the root logger to the debug log.
    app = Flask("planetgen_app", root_path=os.path.dirname(os.path.abspath(__file__)),
                static_folder=web.STATIC_DIR)
    app.config.from_object(config_object)
    app.before_request(_reject_undecodable_query_string)
    app.before_request(_reject_oversized_body)
    limiter.init_app(app)
    app.register_blueprint(bp)
    app.register_blueprint(auth_bp)
    app.register_blueprint(admin_bp)
    app.register_blueprint(population_bp)
    app.register_blueprint(edits_bp)
    app.register_blueprint(workqueue_bp)
    web.init_app(app, limiter=limiter)
    app.teardown_appcontext(close_db)
    app.teardown_appcontext(close_control_db)
    _register_error_handlers(app)
    _register_security_headers(app)
    _register_request_logging(app)
    _warn_if_unshared_ratelimit_storage(app)
    _apply_proxy_fix(app)
    return app


def _apply_proxy_fix(app):
    """
    Behind a reverse proxy (nginx, Caddy, IIS, Apache `mod_proxy`), takes
    the client address, scheme and host from the proxy's `X-Forwarded-*`
    headers, trusting as many proxies per header as `PROXY_FIX` says
    (`config.json`'s `proxy_fix`, see `config._proxy_fix`). Rate limits
    are keyed by `request.remote_addr` and HSTS depends on
    `request.is_secure`, so both need this behind a proxy. Nothing
    changes when every count is 0 (the default, right for Apache +
    mod_wsgi).

    The in-process API calls the pages make (`web/transport.py`) pass
    through this too; they carry no `X-Forwarded-*` headers, so they keep
    the `REMOTE_ADDR` the page gave them.
    """
    counts = app.config.get("PROXY_FIX") or {}
    if not any(counts.get(name, 0) for name in ("x_for", "x_proto", "x_host")):
        return
    from werkzeug.middleware.proxy_fix import ProxyFix

    app.wsgi_app = ProxyFix(
        app.wsgi_app,
        x_for=counts.get("x_for", 0),
        x_proto=counts.get("x_proto", 0),
        x_host=counts.get("x_host", 0),
    )
    log.debug(f"Trusting X-Forwarded-* headers from proxies: {counts}")


def _reject_undecodable_query_string():
    """
    First `before_request` hook: a query string holding raw (not
    %-escaped) bytes that aren't UTF-8 is a 400. Werkzeug decodes
    `request.args` strictly, so without this every route that reads a
    parameter would fail with a UnicodeDecodeError (a 500).
    """
    try:
        request.query_string.decode("utf-8")
    except UnicodeDecodeError:
        abort(400, description="The query string is not valid UTF-8.")


def _reject_oversized_body():
    """
    A request whose `Content-Length` is over `MAX_CONTENT_LENGTH` is a 413
    before anything else runs (TEST.47) -- before the login checks and
    the database, not only once a route reads the body. A body sent
    without a length is still cut off at the limit when it's read.
    """
    limit = current_app.config.get("MAX_CONTENT_LENGTH")
    if limit is not None and (request.content_length or 0) > limit:
        abort(413)


_SIGN_IN_PATHS = ("/api/auth/login", "/api/auth/change-credentials")
"""tuple: The routes where a 429 from the rate limiter means password
guessing, logged as `AUTH login.ratelimited`."""

_SECRET_WORDS = ("password", "token", "secret", "key")


def _register_request_logging(app):
    """
    With the debug log on (see `planetgen.util.log`), records every API
    request as it arrives and as it's answered -- method, path, query,
    caller, which credential it carried (never the credential itself),
    status, size and time taken. Request bodies are summarized by their
    field names only; a login's password never reaches the log.
    """
    @app.before_request
    def _log_request_start():
        g.log_request_started = time.perf_counter()
        if not log.debug_log_active():
            return
        args = {k: ("<withheld>" if any(w in k.lower() for w in _SECRET_WORDS) else v)
                for k, v in request.args.lists()}
        credential = ("API key" if request.headers.get("Authorization", "").startswith("Bearer ")
                      else "session cookie" if request.cookies else "none")
        body = ""
        if request.is_json:
            payload = request.get_json(silent=True)
            if isinstance(payload, dict):
                body = f", JSON body fields {sorted(payload)}"
        log.debug(f"API request: {request.method} {request.path} args={args} from {request.remote_addr} "
                  f"(credential: {credential}, user agent {request.headers.get('User-Agent', '?')!r}){body}")

    @app.after_request
    def _log_request_end(response):
        if log.debug_log_active():
            started = g.get("log_request_started")
            elapsed = f" in {(time.perf_counter() - started) * 1000:.1f}ms" if started is not None else ""
            log.debug(f"API response: {request.method} {request.path} -> {response.status} "
                      f"({response.calculate_content_length() or 0} bytes){elapsed}")
        return response


def _register_security_headers(app):
    """
    Adds defense-in-depth response headers to every response. The HTML
    pages (`planetgen/web/`) get `web.SECURITY_HEADERS` (the one source of
    the pages' CSP and the rest). Everything else
    (API JSON, static files under the dev server) gets the same three
    basic headers with `Content-Security-Policy: default-src 'none'` --
    a JSON body never needs to load anything -- except JavaScript, which
    gets the pages' own CSP: a script started as a Web Worker (the Galaxy
    Map's `static/galaxyblocks.js`) runs under its own response's policy,
    and `'none'` would stop it importing its sibling modules. Either way a
    request that came in over HTTPS also gets
    `web.STRICT_TRANSPORT_SECURITY`.
    """
    from planetgen.web import CONTENT_SECURITY_POLICY, SECURITY_HEADERS, STRICT_TRANSPORT_SECURITY

    @app.after_request
    def _add_headers(response):
        if request.is_secure:
            response.headers.setdefault(*STRICT_TRANSPORT_SECURITY)
        if response.mimetype == "text/html":
            for name, value in SECURITY_HEADERS:
                response.headers.setdefault(name, value)
            return response
        response.headers.setdefault("X-Content-Type-Options", "nosniff")
        response.headers.setdefault("X-Frame-Options", "DENY")
        response.headers.setdefault("Referrer-Policy", "no-referrer")
        script = response.mimetype in ("text/javascript", "application/javascript")
        response.headers.setdefault("Content-Security-Policy", CONTENT_SECURITY_POLICY if script else "default-src 'none'")
        return response


def _warn_if_unshared_ratelimit_storage(app):
    """
    The rate limits and login lockouts count on Redis unless
    `RATELIMIT_STORAGE_URI` says `memory://` (SEC.30). This surfaces in the
    running app's own logs that such a choice is per worker process, so a
    deployment that later grows past one worker notices.
    """
    if app.config.get("RATELIMIT_STORAGE_URI") == "memory://":
        app.logger.warning(
            "Flask-Limiter and the login lockouts are counting in memory: "
            "rate limits are tracked per worker process, not shared across "
            "them. This is correct for a single-process deployment (Flask's "
            "dev server, or mod_wsgi/gunicorn with exactly one worker) but "
            "silently under-enforces limits by roughly a factor of the "
            "worker count otherwise -- leave ratelimit.storage_uri empty in "
            "config.json to count on the Redis server (redis.url)."
        )


def _register_error_handlers(app):
    """
    Forces every error response for `/api/...` -- not just the ones
    routes.py already handles explicitly -- through the same `{"error":
    "..."}` JSON shape. Any other path is one of the HTML pages
    (`planetgen/web/`), which get an HTML error page instead
    (`web.errors.render_error`) -- equally free of tracebacks.
    Without this, an unmatched URL or an uncaught exception falls through
    to Flask's default HTML error page, which is the wrong content type
    for a JSON-only API and leaks a stack trace to the client in
    production if `DEBUG` is ever left off by omission rather than intent.
    """

    @app.errorhandler(ApiError)
    def _handle_api_error(exc):
        log.debug(f"API error {exc.status_code} on {request.method} {request.path}: {exc.message}")
        return jsonify({"error": exc.message}), exc.status_code

    from planetgen.web.errors import TIMEOUT_MESSAGE, render_error, unexpected_error_message

    @app.errorhandler(pymysql.err.OperationalError)
    def _handle_statement_timeout(exc):
        # PERF.17: a read that ran past `mysql.statement_timeout_seconds`
        # is a 504 with a plain message; any other database error is
        # still the generic 500.
        if not (exc.args and exc.args[0] in store.STATEMENT_TIMEOUT_ERRORS):
            return _handle_internal_error(exc)
        log.debug(f"Statement time limit hit on {request.method} {request.path}: {exc}")
        if _is_api_request():
            return jsonify({"error": "QUERY_TIMEOUT: the database took too long to answer"}), 504
        return render_error(504, TIMEOUT_MESSAGE)

    @app.errorhandler(400)
    def _handle_bad_request(exc):
        if _is_api_request():
            return jsonify({"error": exc.description or "bad request"}), 400
        return render_error(400, exc.description or "Bad request.")

    @app.errorhandler(404)
    def _handle_not_found(exc):
        if _is_api_request():
            return jsonify({"error": "not found"}), 404
        return render_error(404, "There is no page at this address.")

    @app.errorhandler(405)
    def _handle_method_not_allowed(exc):
        if _is_api_request():
            response = make_response(jsonify({"error": "method not allowed"}), 405)
        else:
            response = make_response(render_error(405, "This page can't be used that way."))
        # A 405 names the methods the URL does take (RFC 9110).
        if getattr(exc, "valid_methods", None):
            response.headers["Allow"] = ", ".join(exc.valid_methods)
        return response

    @app.errorhandler(500)
    def _handle_internal_error(exc):
        # No exception detail in the body -- this fires for genuinely
        # unexpected failures (a bug, a lost DB connection), and echoing
        # str(exc) back to the client risks leaking internals (file
        # paths, query text) that aren't this API's business to expose.
        # The real detail still reaches Flask's own logger.
        # With the debug log on, this also lands there (app.logger
        # propagates to the root logger, which the debug log listens on).
        if _is_api_request():
            app.logger.exception("Unhandled exception in API request")
            return jsonify({"error": "internal server error"}), 500
        app.logger.exception("Unhandled exception while building a page")
        log.exception(f"Unhandled exception while building {request.path}")
        return render_error(500, unexpected_error_message())

    @app.errorhandler(413)
    def _handle_too_large(exc):
        # A body over `MAX_CONTENT_LENGTH` (config.py).
        log.debug(f"Request body too large on {request.method} {request.path}: {request.content_length} bytes")
        if _is_api_request():
            return jsonify({"error": "request body too large"}), 413
        return render_error(413, "That request was too large.")

    @app.errorhandler(429)
    def _handle_rate_limit_exceeded(exc):
        # Flask-Limiter raises flask_limiter.errors.RateLimitExceeded, a
        # werkzeug HTTPException subclass for 429 -- registering by the
        # plain integer code catches it (and anything else that ever
        # raises a bare 429) the same way `404`/`405` above catch Flask's
        # own routing exceptions, without importing Flask-Limiter's
        # exception type here just to reference it once.
        # /galaxy/tiles is fetched by the map's script, which wants JSON
        # like the API; every other page gets the HTML error page.
        log.debug(f"Rate limit exceeded by {request.remote_addr} on {request.path}: {exc.description}")
        if request.path in _SIGN_IN_PATHS:
            # A password guesser past the per-address limit (SEC.20);
            # the form never got as far as reading a username.
            activity_log.event("AUTH", "login.ratelimited", path=request.path)
        if not _is_api_request() and request.endpoint != "web.galaxy_tiles":
            return render_error(429, "Too many requests. Please wait a minute and try again.")
        return jsonify({"error": "rate limit exceeded", "detail": exc.description}), 429


if __name__ == "__main__":
    create_app().run(debug=True)
