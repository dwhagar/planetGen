# html/web/csrf.py

"""
CSRF protection for the Flask-served pages' POST forms (signed double
submit, no extra dependency).

How it works: the first page that renders a form gives the browser a
random nonce in the `pg_csrf` cookie (HttpOnly, SameSite=Strict, Secure
unless `SESSION_COOKIE_SECURE` is off). The form carries
`HMAC-SHA256(SECRET_KEY, nonce + SHA-256(admin session cookie))` in a
hidden `csrf_token` field. On any POST/PUT/PATCH/DELETE to a page (any
path outside `/api`), `protect()` recomputes the HMAC from the two
cookies and compares it with the field in constant time; a missing or
wrong token gets a 400 page and the view never runs.

The token is bound to the login session (`api.authz.SESSION_COOKIE_NAME`,
empty when nobody is logged in, as on the login form): a token minted
for one session, or before logging in, fails for any other session.
Tokens are computed per render, so the pages a browser sees after
logging in (or after `/account` re-issues its session) already carry
tokens for the new session; a form left open in another tab from before
that change fails once and works after a reload.

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

from flask import abort, current_app, g, redirect, request, url_for
from markupsafe import Markup, escape

from planetgen.api.authz import SESSION_COOKIE_NAME
from planetgen.admin import activity_log

COOKIE_NAME = "pg_csrf"
FIELD_NAME = "csrf_token"
UNSAFE_METHODS = frozenset({"POST", "PUT", "PATCH", "DELETE"})
SIGN_IN_ENDPOINTS = frozenset({"web.login", "web.login_code"})
"""frozenset: The sign-in forms' endpoints (SEC.31): see `protect`."""


def _session():
    """This request's admin session cookie (`""` when not logged in)."""
    return request.cookies.get(SESSION_COOKIE_NAME) or ""


def _sign(nonce, session=""):
    """The token for `nonce` bound to the admin session cookie value
    `session` (`""` for none). Only a hash of the session goes into the
    HMAC, so its length and characters don't matter."""
    key = current_app.config["SECRET_KEY"]
    if isinstance(key, str):
        key = key.encode("utf-8")
    session_hash = hashlib.sha256(str(session or "").encode("utf-8", errors="surrogatepass")).hexdigest()
    message = f"{nonce}\n{session_hash}".encode("utf-8", errors="surrogatepass")
    return hmac.new(key, message, hashlib.sha256).hexdigest()


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
    """The token value for this request's forms (bound to this request's
    login session)."""
    return _sign(_nonce(), _session())


def csrf_field():
    """A hidden `<input>` carrying `csrf_token()` (a template global)."""
    return Markup(f'<input type="hidden" name="{FIELD_NAME}" value="{escape(csrf_token())}">')


def valid(submitted, nonce, session=""):
    """Whether `submitted` is the right token for `nonce` and the admin
    session cookie value `session`."""
    if not submitted or not nonce:
        return False
    # Bytes, not str: compare_digest refuses a str with non-ASCII
    # characters (TypeError), and a submitted token can hold anything.
    return hmac.compare_digest(str(submitted).encode("utf-8", errors="surrogatepass"),
                               _sign(nonce, session).encode("ascii"))


def protect():
    """App-wide `before_request` hook: rejects an unsafe request to any
    page (anything outside `/api`) whose form token doesn't match its
    cookie and login session."""
    if request.method not in UNSAFE_METHODS:
        return None
    if request.path == "/api" or request.path.startswith("/api/"):
        return None
    submitted = request.form.get(FIELD_NAME) or request.headers.get("X-CSRF-Token")
    if not valid(submitted, request.cookies.get(COOKIE_NAME), _session()):
        activity_log.event("AUTHZ", "csrf.failed", path=request.path)
        if (request.endpoint in SIGN_IN_ENDPOINTS and _session()
                and valid(submitted, request.cookies.get(COOKIE_NAME), "")):
            # SEC.31: a sign-in form sent again after the browser already
            # signed in (a second click or Enter, or a password manager
            # submitting too) carries the token minted before sign-in,
            # which no longer matches the new session (and only that
            # token: no token, or any other, is still a 400). Nothing is done
            # with it; the login page forwards a signed-in admin to
            # `next`, or shows the form again if the session is stale.
            return redirect(url_for("web.login", next=request.values.get("next") or None), code=303)
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
