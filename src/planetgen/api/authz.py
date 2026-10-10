# planetgen/api/authz.py

"""
Request-level admin authentication/authorization: resolving the calling
admin from either the session cookie (browser) or an `Authorization:
Bearer <api-key>` header (programmatic callers), the `require_admin`
route decorator built on that, and the audit-log helper every write
route calls after succeeding. See `planetgen/admin/auth.py` for the
actual credential/session/key logic this wraps.
"""

from functools import wraps

from flask import g, request

from planetgen.admin import activity_log, auth as adminAuth
from planetgen.util import log

from . import ids
from .common import ApiError, get_control_db

SESSION_COOKIE_NAME = "pg_admin_session"

_BEARER_PREFIX = "Bearer "


def api_key_admin(raw_key):
    """
    `adminAuth.validate_api_key` for `raw_key`, looked up once per request:
    the rate limiter asks before `require_admin` does, and both get the
    same answer (API.9).
    """
    cached = getattr(g, "_api_key_lookup", None)
    if cached is not None and cached[0] == raw_key:
        return cached[1]
    admin = adminAuth.validate_api_key(get_control_db(), raw_key)
    g._api_key_lookup = (raw_key, admin)
    return admin


def request_api_key_id():
    """The id of the valid API key this request carries, else `None`
    (no `Authorization: Bearer` header, or a key that doesn't work). Used
    by the rate limiter, which runs before `require_admin`."""
    header = request.headers.get("Authorization", "")
    if not header.startswith(_BEARER_PREFIX):
        return None
    admin = api_key_admin(header[len(_BEARER_PREFIX):].strip())
    return admin["api_key_id"] if admin else None


def _current_admin():
    """
    Resolves the calling admin, if any: an `Authorization: Bearer <key>`
    header (an API key) takes priority over the session cookie when both
    are somehow present, since a caller sending a header did so
    deliberately.

    Returns:
        dict or None: The `admin_users` row, or `None` if neither
            credential is present/valid.
    """
    conn = get_control_db()
    auth_header = request.headers.get("Authorization", "")
    if auth_header.startswith(_BEARER_PREFIX):
        return api_key_admin(auth_header[len(_BEARER_PREFIX):].strip())
    return adminAuth.validate_session(conn, request.cookies.get(SESSION_COOKIE_NAME))


def require_admin(fresh=False, session_only=False, scope="admin"):
    """
    Route decorator: requires a valid session cookie or API key, storing
    the resolved admin on `g.admin_user` for the view (and for `audit`
    below) to use.

    Args:
        fresh (bool): If `True`, also requires the admin's
            `must_change_credentials` flag to be clear -- every write/
            admin route except `/api/auth/me`, `/api/auth/logout`, and
            `/api/auth/change-credentials` itself sets this, so the
            seeded first `admin` login (random password printed once by
            `planetgen.cli.migrate`) can authenticate but
            can't do anything else until credentials are actually
            changed (see `adminAuth`'s module docstring and
            `control_schema.sql`'s `admin_users` comment).
        session_only (bool): If `True`, an API key isn't enough: the
            route manages the account itself (new API keys, credentials,
            two-factor sign-in, logging out), so it needs a signed-in
            browser session (TEST.44). Otherwise a leaked key could mint
            a fresh key and outlive its own revocation.
        scope (str): The scope an API key must hold (API.9; `adminAuth.SCOPES`,
            where `admin` implies all and `generate` and `upload` imply
            `read`). A signed-in admin has every scope. Every route that
            changes anything needs `admin` unless it says otherwise.

    Raises (at request time, via `ApiError`):
        401: No valid session/API key.
        403: The key lacks `scope`; `fresh=True` and the admin still has the default,
            never-rotated credentials; or `session_only=True` and the
            caller used an API key.
    """
    def decorator(view):
        @wraps(view)
        def wrapped(*args, **kwargs):
            admin = _current_admin()
            bearer = request.headers.get("Authorization", "").startswith(_BEARER_PREFIX)
            via = "API key" if bearer else "session cookie"
            if admin is None:
                log.debug(f"Access to {request.path} denied: no valid {via} (401)")
                _log_refusal(bearer)
                raise ApiError("authentication required", status_code=401)
            if session_only and bearer:
                log.debug(f"Access to {request.path} denied for {admin['username']!r}: API keys can't manage "
                          f"the account (403)")
                activity_log.event("AUTHZ", "apikey.refused", user=admin["username"], path=request.path)
                raise ApiError("an API key can't do this; sign in with your password instead", status_code=403)
            if bearer and not adminAuth.scope_allows(admin["api_scopes"], scope):
                log.debug(f"Access to {request.path} denied for {admin['username']!r}: the key lacks the "
                          f"{scope!r} scope (403)")
                activity_log.event("AUTHZ", "apikey.scope", user=admin["username"], path=request.path,
                                   key_id=admin["api_key_id"], key_prefix=admin["api_key_prefix"], required=scope)
                raise ApiError(f"this key lacks the '{scope}' scope", status_code=403,
                               extra={"required_scope": scope})
            if fresh and admin["must_change_credentials"]:
                log.debug(f"Access to {request.path} denied for {admin['username']!r}: default credentials "
                          f"not changed yet (403)")
                activity_log.event("AUTHZ", "credentials.unchanged", user=admin["username"], path=request.path)
                raise ApiError(
                    "default credentials must be changed (POST /api/auth/change-credentials) "
                    "before this action is allowed",
                    status_code=403,
                )
            g.admin_user = admin
            g.api_key_id = admin["api_key_id"] if bearer else None
            log.debug(f"Access to {request.path} granted to admin {admin['username']!r} via {via}")
            ids.resolve_kwargs(kwargs)  # after the checks above: who may ask comes before what exists
            return view(*args, **kwargs)
        # Read by the route sweep test (tests/test_api_auth_sweep.py) so
        # a new route can't quietly skip the check.
        wrapped.admin_required = {"fresh": fresh, "session_only": session_only, "scope": scope}
        return wrapped
    return decorator


def _log_refusal(bearer):
    """
    The activity-log line for a 401 (SEC.28): `apikey.invalid` for an
    unknown or revoked API key, `session.invalid` for a session cookie
    that names no live session (expired, logged out, or forged), and
    `login.required` when an admin route got neither. `/api/auth/me` with
    no credential at all is how a page asks "is anyone logged in?", so it
    isn't logged.
    """
    if bearer:
        activity_log.event("AUTHZ", "apikey.invalid", path=request.path)
    elif request.cookies.get(SESSION_COOKIE_NAME):
        activity_log.event("AUTHZ", "session.invalid", path=request.path)
    elif request.path != "/api/auth/me":
        activity_log.event("AUTHZ", "login.required", path=request.path)


def audit(action, target=None, detail=None):
    """
    Writes one `admin_audit_log` row for the current request's already-
    authenticated admin (`g.admin_user`, set by `require_admin` above --
    every write route calls this only after `require_admin`/
    `require_fresh_credentials` has run and only after the action itself
    has actually succeeded, never speculatively).

    Args:
        action (str): Short machine-readable label, e.g. `"sector.create"`.
        target (str, optional): What was acted on, e.g. `"sector:42"`.
        detail (str, optional): Free-form extra context.

    Also writes the change to the activity log (`DB <action>`, SEC.28).
    """
    admin = g.admin_user
    adminAuth.record_audit(get_control_db(), admin["id"], admin["username"], action, target=target, detail=detail)
    activity_log.event("DB", action, user=admin["username"], target=target, detail=detail)
