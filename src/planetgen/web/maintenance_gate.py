# planetgen/web/maintenance_gate.py

"""
The request hook for maintenance mode (UX.112; `admin/maintenance.py` holds
the switch). While it is on, anything but a signed-in admin or an admin's
API key gets a 503: the "being updated" page for a page request, JSON for
`/api`. The sign-in pages, the sign-in calls and the static files stay
reachable, so an admin can still log in and the maintenance page can load
its own styling.

It is the app's first `before_request` hook, so nothing behind it (the
database, the ID resolver, the CSRF check) runs for a refused visitor. The
pages' own in-process API calls skip it: the page request that made them
has already passed.
"""

from flask import jsonify, make_response, redirect, render_template, request, url_for

from planetgen.admin import maintenance
from planetgen.api.limiter import is_in_process_call

RETRY_AFTER_SECONDS = 600
"""int: What `Retry-After` says, and how often the maintenance page asks again."""

MESSAGE = ("The site is being updated and is unavailable for now: the galaxy's sector maps are being rebuilt.")
"""str: What a refused visitor is told."""

OPEN_PATHS = frozenset({
    "/maintenance", "/login", "/login/code", "/logout",
    "/api/auth/login", "/api/auth/login/totp", "/api/auth/logout", "/api/auth/me",
})
"""frozenset: The paths a visitor can still reach while maintenance is on."""


def _is_api():
    return request.path == "/api" or request.path.startswith("/api/")


def is_admin():
    """Whether the request carries a valid admin session or API key. A control database that doesn't
    answer counts as not an admin: the safe mistake is a closed site."""
    from planetgen.api.authz import _current_admin

    try:
        return _current_admin() is not None
    except Exception:  # noqa: BLE001
        return False


def page_response(state, status=503):
    """The "being updated" page (a complete page of its own that reads no database): 503 for a visitor,
    200 for an admin looking at what visitors see."""
    response = make_response(render_template(
        "maintenance.html", state=state, message=state.get("reason") or MESSAGE, refresh=RETRY_AFTER_SECONDS, admin=status == 200), status)
    if status == 503:
        response.headers["Retry-After"] = str(RETRY_AFTER_SECONDS)
    response.headers["Cache-Control"] = "no-store"
    return response


def gate():
    """The `before_request` hook: `None` lets the request through, a response stops it."""
    state = maintenance.read()
    if not state["enabled"] or is_in_process_call():
        return None
    path = request.path
    if path in OPEN_PATHS or path.startswith("/static/") or is_admin():
        return None
    if _is_api():
        response = jsonify({"error": state["reason"] or MESSAGE, "maintenance": True, "eta": state["eta"] or None})
        response.status_code = 503
        response.headers["Retry-After"] = str(RETRY_AFTER_SECONDS)
        response.headers["Cache-Control"] = "no-store"
        return response
    return redirect(url_for("web.maintenance_page"), code=302)


def install(app):
    """Registers `gate` on `app`; call it before any other `before_request` hook."""
    app.before_request(gate)
