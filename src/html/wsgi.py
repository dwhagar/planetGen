# html/wsgi.py

"""
mod_wsgi/gunicorn entry point for the read-only planetGen API.

Points Apache's `WSGIScriptAlias` (or a `gunicorn wsgi:application`
invocation) at this file's `application` object. See `docs/api.md` for
the full deployment story alongside the existing `examples/apache/` vhost
that serves the `html/` CGI browser this file now lives alongside. Also
runnable directly (`python src/html/wsgi.py`) to start Flask's own dev
server locally.

This file lives alongside `api/` under `html/`, so Python's own
sys.path[0] (the running script's directory) already makes `api`
importable -- no shim needed for that. `api/routes.py`/`api/config.py`
also import top-level `queryDb`/`stellarObjects`, which live one level up
at `src/` (src layout, same as every CGI script under this directory
reaches for them -- see e.g. `nav.py`), so that directory is added
explicitly below.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api.app import create_app

application = create_app()

if __name__ == "__main__":
    application.run(debug=True)
