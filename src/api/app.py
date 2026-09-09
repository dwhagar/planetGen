# api/app.py

"""
Flask application factory for the planetGen API.

See `docs/api.md` for how to run this in development and how it deploys
behind the project's existing Apache2 vhost (`examples/apache/`).
"""

from flask import Flask, jsonify

from .config import Config
from .limiter import limiter
from .routes import ApiError, bp, close_db


def create_app(config_object=Config):
    app = Flask(__name__)
    app.config.from_object(config_object)
    limiter.init_app(app)
    app.register_blueprint(bp)
    app.teardown_appcontext(close_db)
    _register_error_handlers(app)
    return app


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
        return jsonify({"error": "rate limit exceeded", "detail": exc.description}), 429


if __name__ == "__main__":
    create_app().run(debug=True)
