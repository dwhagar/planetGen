# html/web/csrf.py

"""
CSRF protection for the Flask-served pages' POST forms (signed double
submit, no extra dependency).

How it works: the first page that renders a form gives the browser a
random nonce in the `pg_csrf` cookie (HttpOnly, SameSite=Strict, Secure
unless `SESSION_COOKIE_SECURE` is off). The form carries
`HMAC-SHA256(SECRET_KEY, nonce)` in a hidden `csrf_token` field. On any
POST/PUT/PATCH/DELETE to a page (any path outside `/api`), `protect()` recomputes the HMAC
from the cookie and compares it with the field in constant time; a
missing or wrong token gets a 400 page and the view never runs.

An attacker's page can make the browser send the cookie (SameSite=Strict
already stops that for cross-site requests), but it can neither read the
cookie nor compute the HMAC without the app's secret, so it cannot fill
in the field.

Use in a template:

    <form method="post" action="{{ url_for('web.login') }}">
      {{ csrf_field() }}
      ...
    </form>

Nothing else is needed: the check runs app-wide for every unsafe request
outside `/api` (registered by `web.init_app`), so a new form route is
covered without doing anything. The JSON API under `/api` is not affected (it has its own CSRF defense: JSON
bodies plus a SameSite=Strict session cookie, see `api/auth.py`).
"""

import hashlib
import hmac
import secrets

from flask import abort, current_app, g, request
from markupsafe import Markup, escape

COOKIE_NAME = "pg_csrf"
FIELD_NAME = "csrf_token"
UNSAFE_METHODS = frozenset({"POST", "PUT", "PATCH", "DELETE"})


def _sign(nonce):
    key = current_app.config["SECRET_KEY"]
    if isinstance(key, str):
        key = key.encode("utf-8")
    return hmac.new(key, nonce.encode("utf-8"), hashlib.sha256).hexdigest()


def _nonce():
    """This request's nonce: the browser's cookie, or a new one that
    `set_cookie` sends back with the response."""
    nonce = g.get("csrf_nonce")
    if nonce is None:
        nonce = request.cookies.get(COOKIE_NAME) or ""
        if len(nonce) < 32:
            nonce = secrets.token_urlsafe(32)
            g.csrf_new_nonce = True
        g.csrf_nonce = nonce
    return nonce


def csrf_token():
    """The token value for this request's forms."""
    return _sign(_nonce())


def csrf_field():
    """A hidden `<input>` carrying `csrf_token()` (a template global)."""
    return Markup(f'<input type="hidden" name="{FIELD_NAME}" value="{escape(csrf_token())}">')


def valid(submitted, nonce):
    """Whether `submitted` is the right token for `nonce`."""
    if not submitted or not nonce:
        return False
    return hmac.compare_digest(str(submitted), _sign(nonce))


def protect():
    """App-wide `before_request` hook: rejects an unsafe request to any
    page (anything outside `/api`) whose form token doesn't match its
    cookie."""
    if request.method not in UNSAFE_METHODS:
        return None
    if request.path == "/api" or request.path.startswith("/api/"):
        return None
    submitted = request.form.get(FIELD_NAME) or request.headers.get("X-CSRF-Token")
    if not valid(submitted, request.cookies.get(COOKIE_NAME)):
        abort(400, description="This form has expired or was not sent from this site. "
                               "Go back, reload the page and try again.")
    return None


def set_cookie(response):
    """`after_request` hook: sends a newly made nonce to the browser."""
    if g.get("csrf_new_nonce"):
        response.set_cookie(
            COOKIE_NAME, g.csrf_nonce,
            httponly=True, samesite="Strict", path="/",
            secure=current_app.config.get("SESSION_COOKIE_SECURE", True),
        )
    return response
