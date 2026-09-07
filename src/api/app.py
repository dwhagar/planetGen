# api/app.py

"""
Flask application factory for the read-only planetGen API.

See `api/README.md` for how to run this in development and how it deploys
behind the project's existing Apache2 vhost (`apache/`).
"""

from flask import Flask

from .config import Config
from .routes import bp, close_db


def create_app(config_object=Config):
    app = Flask(__name__)
    app.config.from_object(config_object)
    app.register_blueprint(bp)
    app.teardown_appcontext(close_db)
    return app


if __name__ == "__main__":
    create_app().run(debug=True)
