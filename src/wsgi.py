# src/wsgi.py

"""
mod_wsgi/gunicorn entry point for the read-only planetGen API.

Points Apache's `WSGIScriptAlias` (or a `gunicorn wsgi:application`
invocation) at this file's `application` object. See `docs/api.md` for
the full deployment story alongside the existing `examples/apache/` vhost
that serves the `html/` CGI browser. Also runnable directly
(`python src/wsgi.py`) to start Flask's own dev server locally.

This file lives alongside `api/` under `src/`, so Python's own
sys.path[0] (the running script's directory) already makes `api`
importable -- no sys.path shim needed.
"""

from api.app import create_app

application = create_app()

if __name__ == "__main__":
    application.run(debug=True)
