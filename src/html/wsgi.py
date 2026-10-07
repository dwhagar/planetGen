# html/wsgi.py

"""
WSGI entry point for the planetGen Flask app: the JSON API under `/api`
and the HTML pages (`planetgen.web`, served at `/`, `/sectors`, ...).

Apache's `WSGIScriptAlias` points at this file's `application` object;
gunicorn (`gunicorn --pythonpath .../src/html wsgi:application`) and
waitress (`waitress-serve wsgi:application` run from `src/html`) load it
the same way. Behind a separate reverse proxy (nginx, Caddy, IIS), set
`config.json`'s `proxy_fix` so the app sees the client's address and
scheme (`planetgen/web/app.py`'s `_apply_proxy_fix`). See `docs/deployment/` for
every platform. Also runnable directly (`python src/html/wsgi.py`) to
start Flask's own dev server locally.

The planetgen package lives in `src/`, one folder up. mod_wsgi's
`WSGIScriptAlias` loads this file as a WSGI script module and doesn't
reliably put any folder on `sys.path` first, so `src/` is added
explicitly below. (Step 14 of the package move replaces this with an
editable install.)
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from planetgen.web.app import create_app
from planetgen.util import log

# The API's log lines are named "api" in the debug log; mod_wsgi owns
# stdout, so nothing goes to the console.
log.set_component("api")
log.configure(log.NORMAL, console=False)


def _restore_mount_prefix(wsgi_app):
    """
    Undoes what a `WSGIScriptAlias /api ...` mount does to each request
    before it reaches Flask.

    The example vhost now mounts the app at `/` (`WSGIScriptAlias /
    .../wsgi.py`, so the same app serves the HTML pages too). Under that
    mount `SCRIPT_NAME` is empty and this is a no-op. It stays for a
    server still using the older `/api` mount, which keeps working for
    the API (the HTML pages need the `/` mount to be reachable at all).
    The history below explains the `/api` case.

    mod_wsgi mounts this app the same way Apache's `ScriptAlias` mounts a
    CGI script: a request for `/api/health` arrives here with
    `SCRIPT_NAME="/api"` and `PATH_INFO="/health"`, the `/api` already
    stripped off. But every route in `api/routes.py`/`api/auth.py` is
    itself registered under a blueprint `url_prefix` of `/api`(`/auth`) --
    chosen so the *documented* URLs (`docs/api.md`) and a direct
    `python wsgi.py` dev-server run (where Flask's own server puts the
    whole path in `PATH_INFO` and leaves `SCRIPT_NAME` empty) both match
    without a mount in front. Under the real mod_wsgi mount, Werkzeug's
    routing matches a rule's full path against `PATH_INFO` alone --
    `SCRIPT_NAME` isn't reattached first -- so `/api/health`'s rule never
    matches the already-stripped `/health`, and every `/api/*` request
    404s from Flask's own handler (confirmed in production: `curl
    https://.../api/health` -> 404 `{"error": "not found"}`, the JSON
    shape and headers from `api/app.py`'s own error handler, not Apache's
    default 404 page).

    Reattaching `SCRIPT_NAME` onto `PATH_INFO` and clearing `SCRIPT_NAME`
    restores the full path mod_wsgi took apart, so the blueprint's
    `/api`-prefixed rules match again. A no-op wherever `SCRIPT_NAME` is
    already empty (the dev server, `python -m pytest`'s `test_client()`),
    so this only changes behavior under the real Apache mount.
    """

    def middleware(environ, start_response):
        environ["PATH_INFO"] = environ.get("SCRIPT_NAME", "") + environ.get("PATH_INFO", "")
        environ["SCRIPT_NAME"] = ""
        return wsgi_app(environ, start_response)

    return middleware


application = create_app()
application.wsgi_app = _restore_mount_prefix(application.wsgi_app)

if __name__ == "__main__":
    application.run(debug=True)
