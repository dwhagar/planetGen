# html/api/app.py

"""
Flask application factory for the planetGen API.

See `docs/api.md` for how to run this in development and how it deploys
behind the project's existing Apache2 vhost (`examples/apache/`).
"""

import time

from flask import Flask, g, jsonify, request

from stellarObjects import log

from .admin import bp as admin_bp
from .auth import bp as auth_bp
from .common import ApiError, close_control_db
from .config import Config
from .limiter import limiter
from .routes import bp, close_db


def create_app(config_object=Config):
    app = Flask(__name__)
    app.config.from_object(config_object)
    limiter.init_app(app)
    app.register_blueprint(bp)
    app.register_blueprint(auth_bp)
    app.register_blueprint(admin_bp)
    app.teardown_appcontext(close_db)
    app.teardown_appcontext(close_control_db)
    _register_error_handlers(app)
    _register_security_headers(app)
    _register_request_logging(app)
    _warn_if_unshared_ratelimit_storage(app)
    return app


_SECRET_WORDS = ("password", "token", "secret", "key")


def _register_request_logging(app):
    """
    With the debug log on (see `stellarObjects.log`), records every API
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
    Adds a handful of defense-in-depth response headers this JSON-only API
    has no legitimate reason to omit -- `default-src 'none'` in particular
    is safe here specifically because every response is `application/
    json` (see `routes.py`), never HTML/JS that would need to load
    anything of its own.
    """

    @app.after_request
    def _add_headers(response):
        response.headers.setdefault("X-Content-Type-Options", "nosniff")
        response.headers.setdefault("X-Frame-Options", "DENY")
        response.headers.setdefault("Referrer-Policy", "no-referrer")
        response.headers.setdefault("Content-Security-Policy", "default-src 'none'")
        return response


def _warn_if_unshared_ratelimit_storage(app):
    """
    `config.py`'s own docstring already explains why `RATELIMIT_STORAGE_URI`
    defaults to `memory://` and when that stops being accurate (more than
    one worker process); this surfaces the same warning in the running
    app's own logs, at the moment it's actually decided, so a deployment
    that later grows past one worker doesn't have to go re-read that
    docstring to notice the mismatch.
    """
    if app.config.get("RATELIMIT_STORAGE_URI") == "memory://":
        app.logger.warning(
            "Flask-Limiter is using in-memory storage: rate limits are "
            "tracked per worker process, not shared across them. This is "
            "correct for a single-process deployment (Flask's dev server, "
            "or mod_wsgi/gunicorn with exactly one worker) but silently "
            "under-enforces limits by roughly a factor of the worker count "
            "otherwise -- set PLANETGEN_RATELIMIT_STORAGE_URI to a shared "
            "backend (e.g. redis://...) before scaling past one worker."
        )


def _register_error_handlers(app):
    """
    Forces every error response -- not just the ones routes.py already
    handles explicitly -- through the same `{"error": "..."}` JSON shape.
    Without this, an unmatched URL or an uncaught exception falls through
    to Flask's default HTML error page, which is the wrong content type
    for a JSON-only API and leaks a stack trace to the client in
    production if `DEBUG` is ever left off by omission rather than intent.
    """

    @app.errorhandler(ApiError)
    def _handle_api_error(exc):
        log.debug(f"API error {exc.status_code} on {request.method} {request.path}: {exc.message}")
        return jsonify({"error": exc.message}), exc.status_code

    @app.errorhandler(404)
    def _handle_not_found(exc):
        return jsonify({"error": "not found"}), 404

    @app.errorhandler(405)
    def _handle_method_not_allowed(exc):
        return jsonify({"error": "method not allowed"}), 405

    @app.errorhandler(500)
    def _handle_internal_error(exc):
        # No exception detail in the body -- this fires for genuinely
        # unexpected failures (a bug, a lost DB connection), and echoing
        # str(exc) back to the client risks leaking internals (file
        # paths, query text) that aren't this API's business to expose.
        # The real detail still reaches Flask's own logger.
        # With the debug log on, this also lands there (app.logger
        # propagates to the root logger, which the debug log listens on).
        app.logger.exception("Unhandled exception in API request")
        return jsonify({"error": "internal server error"}), 500

    @app.errorhandler(429)
    def _handle_rate_limit_exceeded(exc):
        # Flask-Limiter raises flask_limiter.errors.RateLimitExceeded, a
        # werkzeug HTTPException subclass for 429 -- registering by the
        # plain integer code catches it (and anything else that ever
        # raises a bare 429) the same way `404`/`405` above catch Flask's
        # own routing exceptions, without importing Flask-Limiter's
        # exception type here just to reference it once.
        log.debug(f"Rate limit exceeded by {request.remote_addr} on {request.path}: {exc.description}")
        return jsonify({"error": "rate limit exceeded", "detail": exc.description}), 429


if __name__ == "__main__":
    create_app().run(debug=True)
