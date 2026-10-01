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

import hashlib

from flask import Blueprint, current_app, g, jsonify, request
from itsdangerous import BadSignature, SignatureExpired, URLSafeTimedSerializer

from stellarObjects import activitylog, adminAuth, log

from .authz import SESSION_COOKIE_NAME, audit, require_admin
from .common import ApiError, get_control_db, require_json_body
from .limiter import limiter
from .loginguard import LoginGuard

bp = Blueprint("auth", __name__, url_prefix="/api/auth")

MAX_API_KEY_LABEL_LENGTH = 128
"""int: `admin_api_keys.label` is VARCHAR(128)."""

LOGIN_RATE_LIMIT = "10 per minute"
"""str: Applied to `POST /api/auth/login` on top of the app-wide default
(`config.Config.RATELIMIT_DEFAULT`) -- the one endpoint in this API an
attacker has any reason to hammer (password guessing against the one
login form), so it gets its own tight per-IP limit
regardless of how the global default is configured. Failed logins are
also counted per address and per username in the control database
(`loginguard.py`, `stellarObjects/loginThrottle.py`): 3 from one address
lock it for 5 minutes doubling to a day, and 10 for one username lock it
for 1 s doubling to 15 minutes."""


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


DEVICE_COOKIE_NAME = "pg_admin_device"
"""str: The trusted-device cookie (SEC.22): set by every successful
login, it lets that browser past its username's lock (`loginguard.py`)."""


def trusted_device_for(username):
    """
    Whether this request carries a valid device cookie for `username`
    (case-insensitively). A database error counts as no.
    """
    raw = request.cookies.get(DEVICE_COOKIE_NAME)
    if not raw:
        return False
    try:
        owner = adminAuth.device_username(get_control_db(), raw)
    except Exception:  # noqa: BLE001 -- e.g. admin_devices not created yet
        return False
    return owner is not None and owner.casefold() == username.strip().casefold()


def _issue_device_cookie(resp, conn, admin_user_id):
    """Replaces the browser's device cookie with a fresh one for this
    admin. Skipped (logged) if `admin_devices` can't be written yet."""
    try:
        adminAuth.revoke_device(conn, request.cookies.get(DEVICE_COOKIE_NAME))
        raw = adminAuth.create_device(conn, admin_user_id)
    except Exception as exc:  # noqa: BLE001
        try:
            conn.rollback()
        except Exception:  # noqa: BLE001
            pass
        log.error(f"Could not record a trusted device (run update.sh?): {exc}")
        return
    resp.set_cookie(
        DEVICE_COOKIE_NAME,
        raw,
        max_age=adminAuth.DEVICE_TTL_DAYS * 86400,
        httponly=True,
        secure=current_app.config.get("SESSION_COOKIE_SECURE", True),
        samesite="Strict",
        path="/",
    )


def _clear_session_cookie(resp):
    # Must match _set_session_cookie's own path exactly -- a browser scopes
    # a cookie by (name, domain, path) together, so deleting with a
    # different path wouldn't clear the one actually set on login; it would
    # just set a second, separate empty cookie under the old path.
    resp.delete_cookie(SESSION_COOKIE_NAME, path="/")


@bp.route("/login", methods=["POST"])
@limiter.limit(LOGIN_RATE_LIMIT)
def login():
    """
    `POST /api/auth/login` `{"username": str, "password": str}` -> sets
    the session cookie and returns `{"username", "must_change_credentials"}`.
    A wrong username or password both get the same generic 401 (see
    `adminAuth.authenticate`) -- this endpoint never reveals whether a
    given username exists. An address or username locked by too many
    failures (`loginguard.py`) gets a 429 with `retry_after`, `scope`
    (`ip` or `user`) and `Retry-After`, known username or not.
    """
    body = require_json_body()
    username = body.get("username")
    password = body.get("password")
    if not isinstance(username, str) or not isinstance(password, str) or not username.strip() or not password:
        raise ApiError("username and password are required")
    username = username.strip()

    # Checked before the password, so a guess made during a lock never
    # learns anything, right or wrong.
    guard = LoginGuard(username, trusted_device=trusted_device_for(username))
    refused = guard.refusal()
    if refused is not None:
        return refused

    conn = get_control_db()
    try:
        admin = adminAuth.authenticate(conn, username, password)
    except adminAuth.AuthError as exc:
        guard.failed("login.failed")
        raise ApiError(str(exc), status_code=401)

    if _totp_on(conn, admin["id"]):
        # Right password, second step still to come (SEC.26): no session
        # yet, and the failure counts stay until the code is right too.
        activitylog.event("AUTH", "login.password_ok", user=admin["username"])
        return jsonify({"totp_required": True, "pending": _pending_serializer().dumps(
            {"id": admin["id"], "h": _hash_fingerprint(admin["password_hash"])})})
    guard.succeeded()
    return _finish_login(conn, admin)


def _finish_login(conn, admin):
    activitylog.event("AUTH", "login.ok", user=admin["username"])
    raw_token = adminAuth.create_session(conn, admin["id"])
    resp = jsonify(_admin_public_dict(admin))
    _set_session_cookie(resp, raw_token)
    _issue_device_cookie(resp, conn, admin["id"])
    return resp


PENDING_LOGIN_SECONDS = 300
"""int: How long the second step of a two-factor login may wait after
the password (SEC.26)."""


def _pending_serializer():
    return URLSafeTimedSerializer(current_app.config["SECRET_KEY"], salt="planetgen-totp-login")


def _hash_fingerprint(password_hash):
    """Ties a pending login to the password it was made with, so a
    password change voids it."""
    return hashlib.sha256(password_hash.encode("utf-8")).hexdigest()[:16]


def _totp_on(conn, admin_user_id):
    """Whether this admin has two-factor sign-in on. While `admin_totp`
    doesn't exist yet (before `update.sh`), nobody has."""
    try:
        return adminAuth.totp_enabled(conn, admin_user_id)
    except Exception as exc:  # noqa: BLE001
        try:
            conn.rollback()
        except Exception:  # noqa: BLE001
            pass
        log.error(f"Could not read admin_totp (run update.sh?): {exc}")
        return False


@bp.route("/login/totp", methods=["POST"])
@limiter.limit(LOGIN_RATE_LIMIT)
def login_totp():
    """
    `POST /api/auth/login/totp` `{"pending": str, "code": str}` -- the
    second step of a two-factor login (SEC.26): `pending` from the first
    step's response (good for `PENDING_LOGIN_SECONDS`), `code` a current
    authenticator code or an unused recovery code. Sets the session cookie
    like `login`. A wrong code is a 401 and counts as a failed login
    (`totp.failed`) for the address and username; a stale or forged
    `pending` is a 401 asking to sign in again.
    """
    body = require_json_body()
    pending = body.get("pending")
    code = body.get("code")
    if not isinstance(pending, str) or not isinstance(code, str) or not code.strip():
        raise ApiError("pending and code are required")
    try:
        data = _pending_serializer().loads(pending, max_age=PENDING_LOGIN_SECONDS)
    except SignatureExpired:
        raise ApiError("the sign-in took too long; enter your password again", status_code=401)
    except BadSignature:
        raise ApiError("sign in with your password first", status_code=401)
    conn = get_control_db()
    admin = conn.execute("SELECT * FROM admin_users WHERE id = ?", (data.get("id"),)).fetchone()
    if admin is None or _hash_fingerprint(admin["password_hash"]) != data.get("h"):
        raise ApiError("sign in with your password first", status_code=401)

    guard = LoginGuard(admin["username"], trusted_device=trusted_device_for(admin["username"]))
    refused = guard.refusal()
    if refused is not None:
        return refused
    used = adminAuth.check_second_factor(conn, admin["id"], code)
    if used is None:
        guard.failed("totp.failed")
        raise ApiError("that code isn't right", status_code=401)
    guard.succeeded()
    if used == "recovery":
        activitylog.event("AUTH", "totp.recovery_used", user=admin["username"],
                          left=adminAuth.totp_status(conn, admin["id"])["recovery_codes_left"])
    return _finish_login(conn, admin)


@bp.route("/logout", methods=["POST"])
@require_admin(session_only=True)
def logout():
    """`POST /api/auth/logout` -- ends the current session, clears the cookie."""
    adminAuth.end_session(get_control_db(), request.cookies.get(SESSION_COOKIE_NAME))
    activitylog.event("AUTH", "logout", user=g.admin_user["username"])
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
@limiter.limit(LOGIN_RATE_LIMIT)
@require_admin(session_only=True)
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

    The current password is guarded like a login (SEC.23): the same
    per-address limit (`LOGIN_RATE_LIMIT`), refused with 429 while the
    address or the admin's username is locked, and a wrong one counts as
    a failed login (`password.failed`) for both.
    """
    body = require_json_body()
    guard = LoginGuard(g.admin_user["username"], admin_user_id=g.admin_user["id"],
                       trusted_device=trusted_device_for(g.admin_user["username"]))
    refused = guard.refusal()
    if refused is not None:
        return refused
    current_password = body.get("current_password") or ""
    new_username = body.get("new_username")
    new_password = body.get("new_password") or ""

    conn = get_control_db()
    try:
        adminAuth.change_credentials(conn, g.admin_user["id"], current_password, new_username, new_password)
    except adminAuth.WrongPasswordError as exc:
        guard.failed("password.failed")
        raise ApiError(str(exc), status_code=400)
    except adminAuth.AuthError as exc:
        raise ApiError(str(exc), status_code=400)

    guard.succeeded()
    updated = conn.execute("SELECT * FROM admin_users WHERE id = ?", (g.admin_user["id"],)).fetchone()
    activitylog.event("AUTH", "credentials.changed", user=updated["username"],
                      old_user=g.admin_user["username"] if g.admin_user["username"] != updated["username"] else None)
    raw_token = adminAuth.create_session(conn, updated["id"])
    resp = jsonify(_admin_public_dict(updated))
    _set_session_cookie(resp, raw_token)
    # The change deleted every device of this admin; this browser gets a
    # new one.
    _issue_device_cookie(resp, conn, updated["id"])
    return resp


def _check_current_password(body):
    """Re-checks `current_password` for a two-factor change, guarded like
    change-credentials (SEC.23)."""
    guard = LoginGuard(g.admin_user["username"], admin_user_id=g.admin_user["id"],
                       trusted_device=trusted_device_for(g.admin_user["username"]))
    refused = guard.refusal()
    if refused is not None:
        return refused
    password = body.get("current_password")
    row = get_control_db().execute("SELECT password_hash FROM admin_users WHERE id = ?",
                                   (g.admin_user["id"],)).fetchone()
    if not isinstance(password, str) or row is None or not adminAuth.verify_password(password, row["password_hash"]):
        guard.failed("password.failed")
        raise ApiError("current password is incorrect", status_code=400)
    guard.succeeded()
    return None


@bp.route("/totp", methods=["GET"])
@require_admin(fresh=True)
def totp_status():
    """`GET /api/auth/totp` -- `{"enabled", "recovery_codes_left"}` for the
    calling admin (SEC.26)."""
    return jsonify(adminAuth.totp_status(get_control_db(), g.admin_user["id"]))


@bp.route("/totp/setup", methods=["POST"])
@limiter.limit(LOGIN_RATE_LIMIT)
@require_admin(fresh=True, session_only=True)
def totp_setup():
    """
    `POST /api/auth/totp/setup` `{"current_password"}` -- starts setting up
    an authenticator app: returns `{"secret", "uri", "qr_svg"}` (the QR
    code encodes `uri`). Nothing changes at sign-in until `confirm`.
    """
    from stellarObjects import totp
    body = require_json_body()
    refused = _check_current_password(body)
    if refused is not None:
        return refused
    try:
        secret = adminAuth.begin_totp_setup(get_control_db(), g.admin_user["id"])
    except adminAuth.AuthError as exc:
        raise ApiError(str(exc), status_code=400)
    uri = totp.provisioning_uri(secret, g.admin_user["username"])
    return jsonify({"secret": secret, "uri": uri, "qr_svg": totp.qr_svg(uri)})


@bp.route("/totp/confirm", methods=["POST"])
@limiter.limit(LOGIN_RATE_LIMIT)
@require_admin(fresh=True, session_only=True)
def totp_confirm():
    """
    `POST /api/auth/totp/confirm` `{"code"}` -- turns two-factor sign-in on
    once a code from the app matches; returns `{"recovery_codes": [...]}`,
    shown this once.
    """
    body = require_json_body()
    try:
        codes = adminAuth.confirm_totp_setup(get_control_db(), g.admin_user["id"], body.get("code"))
    except adminAuth.AuthError as exc:
        raise ApiError(str(exc), status_code=400)
    activitylog.event("AUTH", "totp.enabled", user=g.admin_user["username"])
    audit("totp.enable", target=f"admin:{g.admin_user['id']}")
    return jsonify({"recovery_codes": codes})


@bp.route("/totp/disable", methods=["POST"])
@limiter.limit(LOGIN_RATE_LIMIT)
@require_admin(fresh=True, session_only=True)
def totp_disable():
    """
    `POST /api/auth/totp/disable` `{"current_password", "code"}` -- turns
    two-factor sign-in off; needs the password and a current code (or a
    recovery code). Lost both? `src/loginLockouts.py --reset-two-factor`.
    Forgets every trusted device of this admin and gives the caller a
    new one.
    """
    body = require_json_body()
    refused = _check_current_password(body)
    if refused is not None:
        return refused
    conn = get_control_db()
    if adminAuth.totp_enabled(conn, g.admin_user["id"]):
        if adminAuth.check_second_factor(conn, g.admin_user["id"], body.get("code")) is None:
            guard = LoginGuard(g.admin_user["username"], admin_user_id=g.admin_user["id"],
                               trusted_device=trusted_device_for(g.admin_user["username"]))
            guard.failed("totp.failed")
            raise ApiError("that code isn't right", status_code=400)
    adminAuth.disable_totp(conn, g.admin_user["id"])
    # Devices trusted while the second factor was on aren't trusted any
    # more (TEST.46); this browser gets a fresh one, like after a
    # credentials change.
    adminAuth.revoke_devices(conn, g.admin_user["id"])
    activitylog.event("AUTH", "totp.disabled", user=g.admin_user["username"])
    audit("totp.disable", target=f"admin:{g.admin_user['id']}")
    resp = jsonify({"enabled": False})
    _issue_device_cookie(resp, conn, g.admin_user["id"])
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
@require_admin(fresh=True, session_only=True)
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
    activitylog.event("AUTH", "apikey.create", user=g.admin_user["username"], key_id=key_id, label=label)
    return jsonify({"id": key_id, "label": label, "key": raw_key}), 201


@bp.route("/api-keys/<int:key_id>", methods=["DELETE"])
@require_admin(fresh=True)
def revoke_api_key(key_id):
    """`DELETE /api/auth/api-keys/<id>` -- revokes one of the calling
    admin's own keys (never another admin's, see `adminAuth.revoke_api_key`)."""
    revoked = adminAuth.revoke_api_key(get_control_db(), g.admin_user["id"], key_id)
    if not revoked:
        raise ApiError(f"no active API key {key_id} for this admin", status_code=404)
    activitylog.event("AUTH", "apikey.revoke", user=g.admin_user["username"], key_id=key_id)
    return jsonify({"status": "ok"})
