# html/wsgi.py

"""
mod_wsgi/gunicorn entry point for the read-only planetGen API.

Points Apache's `WSGIScriptAlias` (or a `gunicorn wsgi:application`
invocation) at this file's `application` object. See `docs/api.md` for
the full deployment story alongside the existing `examples/apache/` vhost
that serves the `html/` CGI browser this file now lives alongside. Also
runnable directly (`python src/html/wsgi.py`) to start Flask's own dev
server locally.

This file lives alongside `api/` under `html/`, so a plain `python
src/html/wsgi.py` invocation already makes `api` importable for free
(Python inserts the running script's own directory as `sys.path[0]`).
mod_wsgi's `WSGIScriptAlias`, though, doesn't execute this file as a
normal `__main__` script -- it loads it as a WSGI script module via its
own machinery, which does *not* reliably add this file's directory to
`sys.path` first (confirmed in production: `WSGIDaemonProcess`/
`WSGIScriptAlias` here raised `ModuleNotFoundError: No module named
'api'` on this very `from api.app import create_app` line). Both `html/`
(for `api`) and `src/` (for `queryDb`/`stellarObjects`, which
`api/routes.py`/`api/config.py` import -- src layout, same as every CGI
script under this directory reaches for them, see e.g. `nav.py`) are
therefore added explicitly below, rather than leaning on either
interpreter's own implicit sys.path setup.
"""

import os
import sys

_HTML_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HTML_DIR))
sys.path.insert(0, _HTML_DIR)

from api.app import create_app

application = create_app()

if __name__ == "__main__":
    application.run(debug=True)
