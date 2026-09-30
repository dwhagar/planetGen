# html/api/auth.py

"""
Admin authentication endpoints: login/logout, forced credential rotation
off the seeded default, "who am I", and API key management. See
`stellarObjects/adminAuth.py` for the underlying credential/session/key
logic and `authz.py` for the `require_admin` decorator every other
write/admin route in this API uses.

The session cookie (`authz.SESSION_COOKIE_NAME`) is `HttpOnly` (never
readable from JS -- irrelevant to XSS exfiltration), `Secure` (never sent
over plain HTTP -- see `config.SESSION_COOKIE_SECURE`, on by default
since this deployment already terminates TLS via Let's Encrypt/Apache),
and `SameSite=Strict` (never sent on a cross-site request at all, which
-- combined with every state-changing route here requiring a JSON body,
something a cross-site HTML form cannot send without triggering a CORS
preflight this API doesn't allow -- is this API's CSRF defense; no
separate CSRF token scheme is layered on top of that).
"""

from flask import Blueprint, current_app, g, jsonify, request

from stellarObjects import adminAuth

from .authz import SESSION_COOKIE_NAME, require_admin
from .common import ApiError, get_control_db, require_json_body
from .limiter import limiter
from .loginbackoff import backoff

bp = Blueprint("auth", __name__, url_prefix="/api/auth")

MAX_API_KEY_LABEL_LENGTH = 128
"""int: `admin_api_keys.label` is VARCHAR(128)."""

LOGIN_RATE_LIMIT = "10 per minute"
"""str: Applied to `POST /api/auth/login` on top of the app-wide default
(`config.Config.RATELIMIT_DEFAULT`) -- the one endpoint in this API an
attacker has any reason to hammer (password guessing against the one
login form), so it gets its own tight per-IP limit
regardless of how the global default is configured. Failed logins are
also counted per username (`loginbackoff.py`), which locks a username
for a growing time after `loginbackoff.FREE_FAILURES` failures whatever
addresses they came from."""


def _admin_public_dict(admin):
    """
    The subset of an `admin_users` row safe to return to the client --
    never `password_hash`, `id` only where a caller needs to reference it
    (e.g. audit context), which none of the response shapes below do.
    """
    return {
        "username": admin["username"],
        "must_change_credentials": bool(admin["must_change_credentials"]),
    }


def _set_session_cookie(resp, raw_token):
    resp.set_cookie(
        SESSION_COOKIE_NAME,
        raw_token,
        max_age=adminAuth.SESSION_TTL_HOURS * 3600,
        httponly=True,
        secure=current_app.config.get("SESSION_COOKIE_SECURE", True),
        samesite="Strict",
        # Not "/api": the admin pages (/admin, /account, ...,
        # web/admin_pages.py) live at the site root and relay this cookie
        # to the API themselves. A cookie scoped to /api is never
        # attached by the browser to a request for /admin, so a path
        # scoped that narrowly locks every admin page out immediately
        # after a successful login.
        path="/",
    )


def _clear_session_cookie(resp):
    # Must match _set_session_cookie's own path exactly -- a browser scopes
    # a cookie by (name, domain, path) together, so deleting with a
    # different path wouldn't clear the one actually set on login; it would
    # just set a second, separate empty cookie under the old path.
    resp.delete_cookie(SESSION_COOKIE_NAME, path="/")


def _too_many_failures(wait):
    """The 429 for a username locked by `loginbackoff`: the wait in the
    body (`retry_after`) and the standard `Retry-After` header."""
    resp = jsonify({
        "error": f"too many failed logins for this username; try again in {wait} second{'' if wait == 1 else 's'}",
        "retry_after": wait,
    })
    resp.status_code = 429
    resp.headers["Retry-After"] = str(wait)
    return resp


@bp.route("/login", methods=["POST"])
@limiter.limit(LOGIN_RATE_LIMIT)
def login():
    """
    `POST /api/auth/login` `{"username": str, "password": str}` -> sets
    the session cookie and returns `{"username", "must_change_credentials"}`.
    A wrong username or password both get the same generic 401 (see
    `adminAuth.authenticate`) -- this endpoint never reveals whether a
    given username exists. A username locked by too many failures
    (`loginbackoff.py`) gets a 429 with `retry_after` and `Retry-After`,
    known or not.
    """
    body = require_json_body()
    username = body.get("username")
    password = body.get("password")
    if not isinstance(username, str) or not isinstance(password, str) or not username.strip() or not password:
        raise ApiError("username and password are required")
    username = username.strip()

    # Checked before the password, so a guess made during a lock never
    # learns anything, right or wrong.
    use_backoff = current_app.config.get("LOGIN_BACKOFF_ENABLED", True)
    wait = backoff.retry_after(username) if use_backoff else 0
    if wait:
        return _too_many_failures(wait)

    conn = get_control_db()
    try:
        admin = adminAuth.authenticate(conn, username, password)
    except adminAuth.AuthError as exc:
        if use_backoff:
            backoff.record_failure(username)
        raise ApiError(str(exc), status_code=401)
    if use_backoff:
        backoff.record_success(username)

    raw_token = adminAuth.create_session(conn, admin["id"])
    resp = jsonify(_admin_public_dict(admin))
    _set_session_cookie(resp, raw_token)
    return resp


@bp.route("/logout", methods=["POST"])
@require_admin()
def logout():
    """`POST /api/auth/logout` -- ends the current session, clears the cookie."""
    adminAuth.end_session(get_control_db(), request.cookies.get(SESSION_COOKIE_NAME))
    resp = jsonify({"status": "ok"})
    _clear_session_cookie(resp)
    return resp


@bp.route("/me")
@require_admin()
def me():
    """`GET /api/auth/me` -- the calling admin's identity, for the web UI
    (and any API-key caller) to check who's logged in and whether
    credentials still need changing."""
    return jsonify(_admin_public_dict(g.admin_user))


@bp.route("/change-credentials", methods=["POST"])
@require_admin()
def change_credentials():
    """
    `POST /api/auth/change-credentials`
    `{"current_password": str, "new_username": str, "new_password": str}`
    -- requires re-entering the current password regardless of whether
    `must_change_credentials` is set (see `adminAuth.change_credentials`).
    On success, clears `must_change_credentials`, ends every session this
    admin had (`adminAuth.change_credentials` deletes them, so any other
    browser is logged out; API keys stay) and issues this caller a fresh
    session cookie, so the browser that made the change stays logged in.
    """
    body = require_json_body()
    current_password = body.get("current_password") or ""
    new_username = body.get("new_username")
    new_password = body.get("new_password") or ""

    conn = get_control_db()
    try:
        adminAuth.change_credentials(conn, g.admin_user["id"], current_password, new_username, new_password)
    except adminAuth.AuthError as exc:
        raise ApiError(str(exc), status_code=400)

    updated = conn.execute("SELECT * FROM admin_users WHERE id = ?", (g.admin_user["id"],)).fetchone()
    raw_token = adminAuth.create_session(conn, updated["id"])
    resp = jsonify(_admin_public_dict(updated))
    _set_session_cookie(resp, raw_token)
    return resp


@bp.route("/api-keys", methods=["GET"])
@require_admin(fresh=True)
def list_api_keys():
    """`GET /api/auth/api-keys` -- the calling admin's own API keys
    (label/timestamps only, never the key/its hash -- see
    `adminAuth.list_api_keys`)."""
    rows = adminAuth.list_api_keys(get_control_db(), g.admin_user["id"])
    return jsonify({"items": [
        {
            "id": row["id"],
            "label": row["label"],
            # UTC (the connection's zone), with an explicit offset.
            "created_at": row["created_at"].isoformat() + "Z" if row["created_at"] else None,
            "last_used_at": row["last_used_at"].isoformat() + "Z" if row["last_used_at"] else None,
            "revoked_at": row["revoked_at"].isoformat() + "Z" if row["revoked_at"] else None,
        }
        for row in rows
    ]})


@bp.route("/api-keys", methods=["POST"])
@require_admin(fresh=True)
def create_api_key():
    """
    `POST /api/auth/api-keys` `{"label": str}` -> `{"id", "label", "key"}`.
    `key` is the raw API key, returned exactly this once -- only its hash
    is ever persisted (see `adminAuth.create_api_key`), so a caller that
    loses it has no way to recover it and must revoke and create another.
    """
    body = require_json_body()
    label = body.get("label")
    if not isinstance(label, str) or not label.strip():
        raise ApiError("'label' is required")
    label = label.strip()
    if len(label) > MAX_API_KEY_LABEL_LENGTH:
        raise ApiError(f"'label' must be at most {MAX_API_KEY_LABEL_LENGTH} characters")

    key_id, raw_key = adminAuth.create_api_key(get_control_db(), g.admin_user["id"], label)
    return jsonify({"id": key_id, "label": label, "key": raw_key}), 201


@bp.route("/api-keys/<int:key_id>", methods=["DELETE"])
@require_admin(fresh=True)
def revoke_api_key(key_id):
    """`DELETE /api/auth/api-keys/<id>` -- revokes one of the calling
    admin's own keys (never another admin's, see `adminAuth.revoke_api_key`)."""
    revoked = adminAuth.revoke_api_key(get_control_db(), g.admin_user["id"], key_id)
    if not revoked:
        raise ApiError(f"no active API key {key_id} for this admin", status_code=404)
    return jsonify({"status": "ok"})
