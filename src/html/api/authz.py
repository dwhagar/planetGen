# html/api/authz.py

"""
Request-level admin authentication/authorization: resolving the calling
admin from either the session cookie (browser) or an `Authorization:
Bearer <api-key>` header (programmatic callers), the `require_admin`
route decorator built on that, and the audit-log helper every write
route calls after succeeding. See `stellarObjects/adminAuth.py` for the
actual credential/session/key logic this wraps.
"""

from functools import wraps

from flask import g, request

from stellarObjects import adminAuth

from .common import ApiError, get_control_db

SESSION_COOKIE_NAME = "pg_admin_session"

_BEARER_PREFIX = "Bearer "


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
        raw_key = auth_header[len(_BEARER_PREFIX):].strip()
        return adminAuth.validate_api_key(conn, raw_key)
    return adminAuth.validate_session(conn, request.cookies.get(SESSION_COOKIE_NAME))


def require_admin(fresh=False):
    """
    Route decorator: requires a valid session cookie or API key, storing
    the resolved admin on `g.admin_user` for the view (and for `audit`
    below) to use.

    Args:
        fresh (bool): If `True`, also requires the admin's
            `must_change_credentials` flag to be clear -- every write/
            admin route except `/api/auth/me`, `/api/auth/logout`, and
            `/api/auth/change-credentials` itself sets this, so the
            seeded default admin/password login can authenticate but
            can't do anything else until credentials are actually
            changed (see `adminAuth`'s module docstring and
            `control_schema.sql`'s `admin_users` comment).

    Raises (at request time, via `ApiError`):
        401: No valid session/API key.
        403: `fresh=True` and the admin still has the default,
            never-rotated credentials.
    """
    def decorator(view):
        @wraps(view)
        def wrapped(*args, **kwargs):
            admin = _current_admin()
            if admin is None:
                raise ApiError("authentication required", status_code=401)
            if fresh and admin["must_change_credentials"]:
                raise ApiError(
                    "default credentials must be changed (POST /api/auth/change-credentials) "
                    "before this action is allowed",
                    status_code=403,
                )
            g.admin_user = admin
            return view(*args, **kwargs)
        return wrapped
    return decorator


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
    """
    admin = g.admin_user
    adminAuth.record_audit(get_control_db(), admin["id"], admin["username"], action, target=target, detail=detail)
