# wsgi.py

"""
mod_wsgi/gunicorn entry point for the read-only planetGen API.

Points Apache's `WSGIScriptAlias` (or a `gunicorn wsgi:application`
invocation) at this file's `application` object. See `api/README.md` for
the full deployment story alongside the existing `apache/` vhost that
serves the `html/` CGI browser. Also runnable directly (`python wsgi.py`)
to start Flask's own dev server locally.
"""

import os
import sys

# api/ lives at src/api (src layout) -- add src/ to the import path so this
# keeps working without requiring `pip install .` first, matching how
# html/'s CGI scripts fall back to a no-install layout.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))

from api.app import create_app

application = create_app()

if __name__ == "__main__":
    application.run(debug=True)
