# stellarObjects/adminAuth.py

"""
Admin authentication/authorization: password hashing, web-session and
API-key lifecycle, credential rotation, and the audit log -- all against
the control schema (`control_schema.sql`), never the per-galaxy content
schema. See that file's header comment for why the two are separate.

Session tokens and API keys are `secrets.token_urlsafe` values, shown to
the caller exactly once (at login / at key creation) and never persisted
-- only their SHA-256 hash (`_hash_token`) is stored, so a database read
alone can't be used to impersonate either. Passwords are hashed with
`werkzeug.security.generate_password_hash` (salted, whichever algorithm
the installed `werkzeug` version defaults to), imported lazily inside the
two functions that need it rather than at module import time, so this
module (and everything in `stellarObjects` that imports it transitively)
stays importable without the `api` extra installed -- matching this
package's existing boundary where `flask`/`werkzeug` are optional, needed
only by `html/api/`.

Every function here takes an already-open `Connection` (to the control
schema -- see `stellarObjects._db.get_control_connection`) rather than
opening its own, the same division of responsibility `_db.py`'s own
`insert_*` functions already use -- callers (Flask routes) own connection
lifetime/pooling.
"""

import hashlib
import secrets

from . import _db

SESSION_TTL_HOURS = 12
"""int: Fixed lifetime for a web-UI session, computed server-side
(`DATE_ADD(CURRENT_TIMESTAMP, ...)` in `create_session`, not a Python
clock read -- same "don't trust the calling process' clock" principle
`_db.get_orbit_update_elapsed_years` already uses `TIMESTAMPDIFF` for).
No sliding renewal (see `control_schema.sql`'s `admin_sessions` comment)
-- short enough that a stolen cookie has bounded value, long enough not
to constantly re-prompt the handful of admins this is built for."""

MIN_PASSWORD_LENGTH = 12
"""int: Minimum length `validate_password_policy` accepts for a new
password -- applies to every credential change, not just the forced one
off the seeded default."""

MAX_USERNAME_LENGTH = 64
"""int: `admin_users.username` is VARCHAR(64)."""

DEFAULT_ADMIN_USERNAME = "admin"
"""str: The username `bootstrap_control_schema` seeds when `admin_users`
is empty."""

INITIAL_PASSWORD_BYTES = 16
"""int: Entropy of the random first password `bootstrap_control_schema`
seeds (`secrets.token_urlsafe(16)`: 22 characters, comfortably over
`MIN_PASSWORD_LENGTH`). There is no published default password: the
seeded one is printed once by `migrateDb.py` (install/update) and never
stored anywhere but as its hash. `admin_users.must_change_credentials`
still starts `TRUE` for this row, and every write/admin endpoint refuses
to work until it's changed (see `html/api/authz.py`'s `require_admin`)."""


class AuthError(Exception):
    """
    Raised for any credential/session/key/policy failure. Callers (Flask
    routes) turn this into a 400/401/403 as appropriate; the message is
    always safe to return to the client as-is -- `authenticate` in
    particular uses the exact same message for "no such username" and
    "wrong password" so a caller can never use this API to enumerate
    valid admin usernames.
    """


class WrongPasswordError(AuthError):
    """The `AuthError` for a wrong current password on a credential
    change, so the caller can count it as a failed sign-in (SEC.20)."""


def hash_password(password):
    """Returns a salted hash of `password` (`werkzeug.security.
    generate_password_hash`) -- never the plaintext, never a reversible
    encoding."""
    from werkzeug.security import generate_password_hash
    return generate_password_hash(password)


def verify_password(password, password_hash):
    """Returns whether `password` matches a hash `hash_password` produced
    (`werkzeug.security.check_password_hash`, constant-time)."""
    from werkzeug.security import check_password_hash
    return check_password_hash(password_hash, password)


def _hash_token(raw_token):
    """SHA-256 hex digest of a session token or API key -- what's actually
    stored (`admin_sessions.token_hash`/`admin_api_keys.key_hash`), never
    the raw value itself."""
    return hashlib.sha256(raw_token.encode("utf-8")).hexdigest()


def _new_token():
    return secrets.token_urlsafe(32)


def validate_password_policy(password, username=None):
    """
    Raises `AuthError` if `password` fails the minimum policy: at least
    `MIN_PASSWORD_LENGTH` characters and not (case-insensitively) equal
    to `username` -- the weakest, most-guessable choices an admin could
    otherwise pick.
    """
    if not isinstance(password, str) or len(password) < MIN_PASSWORD_LENGTH:
        raise AuthError(f"password must be at least {MIN_PASSWORD_LENGTH} characters")
    if username and password.lower() == username.lower():
        raise AuthError("password must not match the username")


def bootstrap_control_schema(config=None):
    """
    Ensures the control schema's *database* (MySQL schema) itself exists,
    then ensures its tables (DDL -- needs a `CREATE`-capable account; see
    `_db.get_control_connection`'s `ensure_schema` note) and, if
    `admin_users` is empty, seeds an `admin` row with a random first
    password (`INITIAL_PASSWORD_BYTES`) and `must_change_credentials`
    set. Idempotent
    -- safe to call on every deploy (`migrateDb.py` calls this right
    after its own content-schema migration, using the same full-access
    account pointed at the control schema instead).

    Unlike a content schema (this project's own convention assumes that
    database already exists, created once by a human/DBA, before any
    tool connects to it -- pymysql can't even open a connection
    *selecting* a database that doesn't exist yet), the control schema's
    name is a fixed, predictable, deployment-wide constant rather than an
    arbitrary per-galaxy choice, so creating it automatically here (via a
    connection that selects no specific database first, the same pattern
    `_db.list_databases` uses) is a reasonable convenience rather than
    something that needs its own manual step.

    Args:
        config (MySQLConfig, optional): Connection parameters for the
            control schema (typically `_db.control_mysql_config(...)`).
            Defaults to `_db.control_mysql_config()`.

    Returns:
        tuple[str, str] or None: `(username, password)` of the admin
            this call just seeded -- the only time the plaintext exists
            anywhere, so the caller must show it to the operator (see
            `migrateDb.py`) -- or `None` when `admin_users` already had a
            row and nothing was seeded. An older install's admin still on
            the published `admin`/`password` login gets a random password
            the same way (`_rotate_published_default_password`).
    """
    config = config or _db.control_mysql_config()
    unselected = _db.get_connection(
        _db.MySQLConfig(host=config.host, port=config.port, user=config.user, password=config.password, database=""),
        ensure_schema=False,
    )
    try:
        unselected.execute(f"CREATE DATABASE IF NOT EXISTS `{config.database}`")
        unselected.commit()
    finally:
        unselected.close()

    conn = _db.get_control_connection(config, ensure_schema=True)
    try:
        row = conn.execute("SELECT COUNT(*) AS n FROM admin_users").fetchone()
        if row["n"] != 0:
            return _rotate_published_default_password(conn)
        password = secrets.token_urlsafe(INITIAL_PASSWORD_BYTES)
        conn.execute(
            "INSERT INTO admin_users (username, password_hash, must_change_credentials) VALUES (?, ?, 1)",
            (DEFAULT_ADMIN_USERNAME, hash_password(password)),
        )
        conn.commit()
        return DEFAULT_ADMIN_USERNAME, password
    finally:
        conn.close()


_PUBLISHED_DEFAULT_PASSWORD = "password"
"""str: The first password older installs seeded, published in the
repository. Only ever checked against, never set."""


def _rotate_published_default_password(conn):
    """
    An install from before random first passwords may still have its
    seeded admin on the published `admin`/`password` login, waiting for
    its forced change -- claimable by whoever logs in first. Gives any such
    row (still `must_change_credentials`, still verifying against the
    published password) a random password instead, ending its sessions.

    Returns:
        tuple[str, str] or None: `(username, new password)` for the
            operator, as `bootstrap_control_schema` returns a new seed, or
            `None` when no row needed it.
    """
    rows = conn.execute(
        "SELECT id, username, password_hash FROM admin_users WHERE must_change_credentials = 1"
    ).fetchall()
    for row in rows:
        if not verify_password(_PUBLISHED_DEFAULT_PASSWORD, row["password_hash"]):
            continue
        password = secrets.token_urlsafe(INITIAL_PASSWORD_BYTES)
        conn.execute("UPDATE admin_users SET password_hash = ? WHERE id = ?", (hash_password(password), row["id"]))
        conn.execute("DELETE FROM admin_sessions WHERE admin_user_id = ?", (row["id"],))
        conn.commit()
        return row["username"], password
    return None


_dummy_password_hash = None
"""str: A real `hash_password` hash of a random value, checked against
when the submitted username doesn't exist so that path costs the same
hash work as a wrong password (see `authenticate`). Computed at import
when `werkzeug` is installed, else on first use."""


def _get_dummy_password_hash():
    global _dummy_password_hash
    if _dummy_password_hash is None:
        _dummy_password_hash = hash_password(secrets.token_urlsafe(16))
    return _dummy_password_hash


try:
    _get_dummy_password_hash()
except ImportError:  # no `api` extra: nothing here authenticates anyway
    pass


def authenticate(conn, username, password):
    """
    Verifies a username/password pair against `admin_users`, updating
    `last_login_at` on success.

    Args:
        conn (Connection): Open control-schema connection.
        username (str): The submitted username.
        password (str): The submitted password.

    Returns:
        dict: The `admin_users` row on success.

    Raises:
        AuthError: On an unknown username or a wrong password -- same
            message either way (see the class docstring), and the same
            time: an unknown username is still checked against a dummy
            hash, so the response time doesn't reveal which usernames
            exist.
    """
    row = conn.execute("SELECT * FROM admin_users WHERE username = ?", (username,)).fetchone()
    if row is None:
        verify_password(password, _get_dummy_password_hash())
        raise AuthError("invalid username or password")
    if not verify_password(password, row["password_hash"]):
        raise AuthError("invalid username or password")
    conn.execute("UPDATE admin_users SET last_login_at = CURRENT_TIMESTAMP WHERE id = ?", (row["id"],))
    conn.commit()
    # Re-fetched rather than patched onto the pre-update `row` in memory --
    # `last_login_at` is a server-evaluated CURRENT_TIMESTAMP, so only a
    # fresh read reflects the real value the UPDATE just wrote.
    return conn.execute("SELECT * FROM admin_users WHERE id = ?", (row["id"],)).fetchone()


def create_session(conn, admin_user_id):
    """
    Creates a new session row for `admin_user_id`.

    Returns:
        str: The raw session token -- persist only as the cookie value;
            the database only ever sees/stores its hash.
    """
    raw_token = _new_token()
    conn.execute(
        "INSERT INTO admin_sessions (admin_user_id, token_hash, expires_at) "
        "VALUES (?, ?, DATE_ADD(CURRENT_TIMESTAMP, INTERVAL ? HOUR))",
        (admin_user_id, _hash_token(raw_token), SESSION_TTL_HOURS),
    )
    conn.commit()
    return raw_token


def validate_session(conn, raw_token):
    """
    Looks up the admin owning a still-valid (unexpired) session token.

    Returns:
        dict or None: The `admin_users` row, or `None` for a missing/
            expired/unknown token -- never raises, since "not logged in"
            is an ordinary outcome here, not an error. Refreshes
            `last_seen_at` on success (informational only; does not
            extend `expires_at` -- see `SESSION_TTL_HOURS`).
    """
    if not raw_token:
        return None
    token_hash = _hash_token(raw_token)
    row = conn.execute(
        """
        SELECT au.* FROM admin_sessions s
        JOIN admin_users au ON au.id = s.admin_user_id
        WHERE s.token_hash = ? AND s.expires_at > CURRENT_TIMESTAMP
        """,
        (token_hash,),
    ).fetchone()
    if row is None:
        return None
    conn.execute("UPDATE admin_sessions SET last_seen_at = CURRENT_TIMESTAMP WHERE token_hash = ?", (token_hash,))
    conn.commit()
    return row


def end_session(conn, raw_token):
    """Deletes a session row (logout) -- a no-op if it's already gone/expired."""
    if not raw_token:
        return
    conn.execute("DELETE FROM admin_sessions WHERE token_hash = ?", (_hash_token(raw_token),))
    conn.commit()


def create_api_key(conn, admin_user_id, label):
    """
    Creates a new API key for `admin_user_id`.

    Returns:
        tuple[int, str]: `(key_id, raw_key)` -- `raw_key` (`pg_<random>`)
            is shown to the caller exactly once; only its hash is
            persisted. Returning `key_id` from the same insert (rather
            than a caller re-querying "this admin's newest key"
            afterward) avoids a race against a concurrent key creation by
            the same admin picking up the wrong id.
    """
    raw_key = f"pg_{_new_token()}"
    cur = conn.execute(
        "INSERT INTO admin_api_keys (admin_user_id, label, key_hash) VALUES (?, ?, ?)",
        (admin_user_id, label, _hash_token(raw_key)),
    )
    conn.commit()
    return cur.lastrowid, raw_key


def validate_api_key(conn, raw_key):
    """
    Looks up the admin owning a still-active (non-revoked) API key.

    Returns:
        dict or None: The `admin_users` row, or `None` for a missing/
            unknown/revoked key. Refreshes `last_used_at` on success.
    """
    if not raw_key:
        return None
    key_hash = _hash_token(raw_key)
    row = conn.execute(
        """
        SELECT au.* FROM admin_api_keys k
        JOIN admin_users au ON au.id = k.admin_user_id
        WHERE k.key_hash = ? AND k.revoked_at IS NULL
        """,
        (key_hash,),
    ).fetchone()
    if row is None:
        return None
    conn.execute("UPDATE admin_api_keys SET last_used_at = CURRENT_TIMESTAMP WHERE key_hash = ?", (key_hash,))
    conn.commit()
    return row


def list_api_keys(conn, admin_user_id):
    """
    Returns every API key belonging to `admin_user_id`, newest first --
    never the key itself or its hash, only display metadata (`GET
    /api/auth/api-keys`'s response shape).
    """
    return conn.execute(
        "SELECT id, label, created_at, last_used_at, revoked_at FROM admin_api_keys "
        "WHERE admin_user_id = ? ORDER BY created_at DESC",
        (admin_user_id,),
    ).fetchall()


def revoke_api_key(conn, admin_user_id, key_id):
    """
    Revokes one of `admin_user_id`'s own API keys -- never another
    admin's; there's no cross-admin management surface at all (per the
    brief: a few admins, no general user-accounts system, each manages
    only their own keys).

    Returns:
        bool: `True` if a key was revoked, `False` if `key_id` doesn't
            exist, doesn't belong to `admin_user_id`, or was already
            revoked.
    """
    cur = conn.execute(
        "UPDATE admin_api_keys SET revoked_at = CURRENT_TIMESTAMP "
        "WHERE id = ? AND admin_user_id = ? AND revoked_at IS NULL",
        (key_id, admin_user_id),
    )
    conn.commit()
    return cur.rowcount > 0


def change_credentials(conn, admin_user_id, current_password, new_username, new_password):
    """
    Changes an admin's username and password together, clearing
    `must_change_credentials`. Always requires re-proving the current
    password (regardless of whether `must_change_credentials` is set) --
    a session cookie/API key alone is never enough to rotate credentials,
    so a hijacked-but-not-fully-compromised session can't lock the real
    admin out.

    Ends every one of this admin's sessions (`admin_sessions`) in the same
    transaction, so a stolen or forgotten browser session stops working
    the moment the password changes; the caller issues a fresh session
    for the browser that made the change (`POST
    /api/auth/change-credentials` does). API keys are kept: they're
    separate, individually revocable credentials.

    Args:
        conn (Connection): Open control-schema connection.
        admin_user_id (int): The admin making the change (from the
            already-authenticated session/API key -- never caller-supplied).
        current_password (str): Must match the admin's existing password.
        new_username (str): The new username (must be non-empty and not
            already used by another admin).
        new_password (str): Must pass `validate_password_policy`.

    Raises:
        AuthError: On a wrong current password, an empty/duplicate
            `new_username`, or a `new_password` failing policy.
    """
    row = conn.execute("SELECT * FROM admin_users WHERE id = ?", (admin_user_id,)).fetchone()
    if row is None or not isinstance(current_password, str) \
            or not verify_password(current_password, row["password_hash"]):
        raise WrongPasswordError("current password is incorrect")

    if new_username is not None and not isinstance(new_username, str):
        raise AuthError("username must be a string")
    new_username = (new_username or "").strip()
    if not new_username:
        raise AuthError("username must not be empty")
    if len(new_username) > MAX_USERNAME_LENGTH:
        raise AuthError(f"username must be at most {MAX_USERNAME_LENGTH} characters")
    validate_password_policy(new_password, username=new_username)

    existing = conn.execute(
        "SELECT id FROM admin_users WHERE username = ? AND id != ?", (new_username, admin_user_id),
    ).fetchone()
    if existing is not None:
        raise AuthError("username is already in use")

    conn.execute(
        "UPDATE admin_users SET username = ?, password_hash = ?, must_change_credentials = 0 WHERE id = ?",
        (new_username, hash_password(new_password), admin_user_id),
    )
    conn.execute("DELETE FROM admin_sessions WHERE admin_user_id = ?", (admin_user_id,))
    conn.commit()


def record_audit(conn, admin_user_id, admin_username, action, target=None, detail=None):
    """
    Writes one `admin_audit_log` row. Called by every write/admin
    endpoint after it succeeds (see `html/api/routes.py`) -- `admin_username`
    is captured here (denormalized alongside `admin_user_id`) rather than
    joined at read time, so the trail stays readable even if that admin is
    later renamed (see `control_schema.sql`'s comment on this table).

    Args:
        conn (Connection): Open control-schema connection.
        admin_user_id (int): The acting admin's id.
        admin_username (str): The acting admin's username, as of this
            action (not re-read from `admin_users` later).
        action (str): A short machine-readable label, e.g. `"sector.create"`.
        target (str, optional): What was acted on, e.g. `"sector:42"`.
        detail (str, optional): Free-form extra context (e.g. the fields
            changed).
    """
    conn.execute(
        "INSERT INTO admin_audit_log (admin_user_id, admin_username, action, target, detail) "
        "VALUES (?, ?, ?, ?, ?)",
        (admin_user_id, admin_username, action, target, detail),
    )
    conn.commit()


LOGIN_FAILURE_ACTIONS = ("login.failed", "login.locked", "password.failed")
"""tuple: The `admin_audit_log` actions a refused sign-in writes (SEC.20):
a wrong username or password, a locked login, and a wrong current
password on a credential change."""

LOGIN_FAILURE_RETENTION_DAYS = 90
"""int: Failure rows older than this are deleted as new ones are written,
so a flood of failed logins can't grow `admin_audit_log` without bound.
Every other audit row is kept."""


def record_login_failure(conn, action, username, ip=None, detail=None, admin_user_id=None):
    """
    Writes one `admin_audit_log` row for a refused sign-in (one of
    `LOGIN_FAILURE_ACTIONS`) and prunes such rows past
    `LOGIN_FAILURE_RETENTION_DAYS`. Never the password.

    Args:
        conn (Connection): Open control-schema connection.
        action (str): One of `LOGIN_FAILURE_ACTIONS`.
        username (str): The username as typed (cut to
            `MAX_USERNAME_LENGTH`); it may not exist.
        ip (str, optional): The client address, kept as the row's target
            (`ip:<address>`).
        detail (str, optional): Free-form extra context.
        admin_user_id (int, optional): The admin, when known (a wrong
            current password comes from a logged-in admin).
    """
    if action not in LOGIN_FAILURE_ACTIONS:
        raise ValueError(f"not a login failure action: {action!r}")
    name = (username or "")[:MAX_USERNAME_LENGTH] or "-"
    conn.execute(
        "INSERT INTO admin_audit_log (admin_user_id, admin_username, action, target, detail) VALUES (?, ?, ?, ?, ?)",
        (admin_user_id, name, action, f"ip:{ip}" if ip else None, detail),
    )
    placeholders = ", ".join("?" for _ in LOGIN_FAILURE_ACTIONS)
    conn.execute(
        f"DELETE FROM admin_audit_log WHERE action IN ({placeholders}) "
        f"AND created_at < DATE_SUB(CURRENT_TIMESTAMP, INTERVAL ? DAY)",
        (*LOGIN_FAILURE_ACTIONS, LOGIN_FAILURE_RETENTION_DAYS),
    )
    conn.commit()


def recent_login_failures(conn, limit=20):
    """
    The newest refused sign-ins (`LOGIN_FAILURE_ACTIONS` rows), newest
    first, for the admin stats page.

    Returns:
        list[dict]: `{"action", "username", "ip", "created_at"}` rows.
    """
    placeholders = ", ".join("?" for _ in LOGIN_FAILURE_ACTIONS)
    rows = conn.execute(
        f"SELECT action, admin_username, target, created_at FROM admin_audit_log "
        f"WHERE action IN ({placeholders}) ORDER BY created_at DESC, id DESC LIMIT ?",
        (*LOGIN_FAILURE_ACTIONS, int(limit)),
    ).fetchall()
    return [{
        "action": row["action"],
        "username": row["admin_username"],
        "ip": row["target"][3:] if row["target"] and row["target"].startswith("ip:") else None,
        "created_at": row["created_at"],
    } for row in rows]
