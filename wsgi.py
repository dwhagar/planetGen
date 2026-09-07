# wsgi.py

"""
mod_wsgi/gunicorn entry point for the read-only planetGen API.

Points Apache's `WSGIScriptAlias` (or a `gunicorn wsgi:application`
invocation) at this file's `application` object. See `api/README.md` for
the full deployment story alongside the existing `apache/` vhost that
serves the `html/` CGI browser.
"""

from api.app import create_app

application = create_app()
