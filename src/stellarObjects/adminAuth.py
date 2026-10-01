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
`werkzeug.security.generate_password_hash` (salted, PBKDF2-SHA256 with
600,000 rounds: `PASSWORD_HASH_METHOD`), imported lazily inside the
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

import gzip
import hashlib
import os
import re
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


PASSWORD_HASH_METHOD = "pbkdf2:sha256:600000"
"""str: How new password hashes are made (SEC.25): PBKDF2-HMAC-SHA256 with
600,000 rounds, OWASP's recommendation. werkzeug's default until
2026-10-01 was scrypt with N = 2^15 (about 100 ms here, under OWASP's
N = 2^17 minimum); N = 2^17 needs 128 MiB per check, which the web
server's five threads at once could turn into a memory spike on a host
that has already been OOM-killed once, so PBKDF2 at about the same cost
in time and next to nothing in memory was chosen instead. A successful
login re-hashes a password stored with any other method
(`needs_rehash`)."""


def hash_password(password):
    """Returns a salted hash of `password` (`werkzeug.security.
    generate_password_hash` with `PASSWORD_HASH_METHOD`) -- never the
    plaintext, never a reversible encoding."""
    from werkzeug.security import generate_password_hash
    return generate_password_hash(password, method=PASSWORD_HASH_METHOD)


def needs_rehash(password_hash):
    """Whether a stored hash was made with settings other than
    `PASSWORD_HASH_METHOD` (an older install's scrypt, for example)."""
    return not (password_hash or "").startswith(PASSWORD_HASH_METHOD + "$")


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


COMMON_PASSWORDS_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "common_passwords.txt.gz")
"""str: The bundled blocklist (SEC.24): every password of 12 or more
characters in the "xato-net 10 million passwords" top 1,000,000 and the
UK NCSC's 100,000 most common (Have I Been Pwned) lists, as published in
SecLists (https://github.com/danielmiessler/SecLists, MIT licence),
case-folded, one per line, gzipped. Shorter ones are refused by length
anyway. Checked offline: nothing is sent anywhere."""

SITE_WORDS = ("planetgen", "password", "admin")
"""tuple: Words a password may not be built around (with the username):
what's left once they are taken out must be at least
`MIN_LEFT_AFTER_SITE_WORDS` characters."""

MIN_LEFT_AFTER_SITE_WORDS = 8

_common_passwords = None


def _load_common_passwords():
    global _common_passwords
    if _common_passwords is None:
        try:
            with gzip.open(COMMON_PASSWORDS_PATH, "rt", encoding="utf-8") as f:
                _common_passwords = frozenset(line.rstrip("\n") for line in f if line.strip())
        except OSError:
            _common_passwords = frozenset()
    return _common_passwords


def is_common_password(password):
    """Whether `password` (case-insensitively) is on the bundled blocklist."""
    return password.casefold() in _load_common_passwords()


def validate_password_policy(password, username=None):
    """
    Raises `AuthError` if `password` fails the policy: at least
    `MIN_PASSWORD_LENGTH` characters, not (case-insensitively) the
    username, not on the bundled list of common and breached passwords
    (SEC.24, NIST SP 800-63B-4), and not just the username or a site word
    (`SITE_WORDS`) with a few characters added. No composition rules
    (NIST advises against them).
    """
    if not isinstance(password, str) or len(password) < MIN_PASSWORD_LENGTH:
        raise AuthError(f"password must be at least {MIN_PASSWORD_LENGTH} characters")
    if username and password.lower() == username.lower():
        raise AuthError("password must not match the username")
    if is_common_password(password):
        raise AuthError("that password is on a list of common or breached passwords; choose another")
    words = [w for w in (*SITE_WORDS, (username or "").strip().casefold()) if len(w) >= 3]
    left = password.casefold()
    for word in sorted(words, key=len, reverse=True):
        left = left.replace(word, "")
    if len(re.sub(r"\s", "", left)) < MIN_LEFT_AFTER_SITE_WORDS:
        raise AuthError("password is too close to the username or the site's name; choose another")


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
    if needs_rehash(row["password_hash"]):
        # Stored with older settings (SEC.25): re-hash now, while the
        # plaintext is at hand.
        conn.execute("UPDATE admin_users SET password_hash = ? WHERE id = ?", (hash_password(password), row["id"]))
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
    revoke_devices(conn, admin_user_id, commit=False)
    conn.commit()


DEVICE_TTL_DAYS = 90
"""int: How long a trusted-device cookie (SEC.22) lasts. Fixed at
creation, like a session; each successful login from that browser issues
a fresh one."""


def create_device(conn, admin_user_id):
    """
    Records a trusted device for `admin_user_id` (SEC.22), deleting that
    admin's expired ones.

    Returns:
        str: The raw device token, for the cookie only; the database keeps
            its hash.
    """
    raw_token = _new_token()
    conn.execute("DELETE FROM admin_devices WHERE admin_user_id = ? AND expires_at <= CURRENT_TIMESTAMP",
                 (admin_user_id,))
    conn.execute(
        "INSERT INTO admin_devices (admin_user_id, token_hash, expires_at) "
        "VALUES (?, ?, DATE_ADD(CURRENT_TIMESTAMP, INTERVAL ? DAY))",
        (admin_user_id, _hash_token(raw_token), DEVICE_TTL_DAYS),
    )
    conn.commit()
    return raw_token


def device_username(conn, raw_token):
    """
    The username of the admin a still-valid device token belongs to, or
    `None` (missing, unknown or expired token -- never raises).
    """
    if not raw_token or not isinstance(raw_token, str):
        return None
    row = conn.execute(
        "SELECT u.username FROM admin_devices d JOIN admin_users u ON u.id = d.admin_user_id "
        "WHERE d.token_hash = ? AND d.expires_at > CURRENT_TIMESTAMP",
        (_hash_token(raw_token),),
    ).fetchone()
    return row["username"] if row else None


def revoke_device(conn, raw_token):
    """Deletes one device token's row (a replaced cookie)."""
    if raw_token and isinstance(raw_token, str):
        conn.execute("DELETE FROM admin_devices WHERE token_hash = ?", (_hash_token(raw_token),))
        conn.commit()


def revoke_devices(conn, admin_user_id, commit=True):
    """Deletes every trusted device of `admin_user_id` (on a credentials
    change, or from the command line). Returns how many."""
    cur = conn.execute("DELETE FROM admin_devices WHERE admin_user_id = ?", (admin_user_id,))
    if commit:
        conn.commit()
    return cur.rowcount


RECOVERY_CODE_COUNT = 10
"""int: Recovery codes made when two-factor sign-in is turned on."""

_RECOVERY_ALPHABET = "abcdefghjkmnpqrstuvwxyz23456789"


def _new_recovery_code():
    """`xxxxx-xxxxx` from an alphabet without look-alikes (about 49 bits)."""
    chars = "".join(secrets.choice(_RECOVERY_ALPHABET) for _ in range(10))
    return f"{chars[:5]}-{chars[5:]}"


def _normalize_recovery_code(code):
    return code.strip().lower().replace(" ", "").replace("-", "") if isinstance(code, str) else ""


def totp_status(conn, admin_user_id):
    """`{"enabled": bool, "recovery_codes_left": int}` for one admin's
    two-factor sign-in (SEC.26)."""
    row = conn.execute("SELECT enabled_at FROM admin_totp WHERE admin_user_id = ?", (admin_user_id,)).fetchone()
    left = conn.execute(
        "SELECT COUNT(*) AS n FROM admin_recovery_codes WHERE admin_user_id = ? AND used_at IS NULL",
        (admin_user_id,)).fetchone()["n"]
    return {"enabled": bool(row and row["enabled_at"]), "recovery_codes_left": int(left)}


def totp_enabled(conn, admin_user_id):
    row = conn.execute("SELECT enabled_at FROM admin_totp WHERE admin_user_id = ?", (admin_user_id,)).fetchone()
    return bool(row and row["enabled_at"])


def begin_totp_setup(conn, admin_user_id):
    """
    Starts (or restarts) setting up an authenticator app: stores a new
    secret, not yet enabled, and returns it. Refused while two-factor
    sign-in is already on (turn it off first).
    """
    from . import totp
    if totp_enabled(conn, admin_user_id):
        raise AuthError("two-factor sign-in is already on")
    secret = totp.new_secret()
    conn.execute("DELETE FROM admin_totp WHERE admin_user_id = ?", (admin_user_id,))
    conn.execute("INSERT INTO admin_totp (admin_user_id, secret) VALUES (?, ?)", (admin_user_id, secret))
    conn.commit()
    return secret


def pending_totp_secret(conn, admin_user_id):
    """The secret of a setup not yet confirmed, or `None`."""
    row = conn.execute("SELECT secret, enabled_at FROM admin_totp WHERE admin_user_id = ?",
                       (admin_user_id,)).fetchone()
    return row["secret"] if row and not row["enabled_at"] else None


def confirm_totp_setup(conn, admin_user_id, code):
    """
    Turns two-factor sign-in on once `code` matches the pending secret.

    Returns:
        list[str]: The new recovery codes, to show once.

    Raises:
        AuthError: No setup started, or the code doesn't match.
    """
    from . import totp
    secret = pending_totp_secret(conn, admin_user_id)
    if secret is None:
        raise AuthError("start setting up two-factor sign-in first")
    step = totp.verify(secret, code)
    if step is None:
        raise AuthError("that code doesn't match; check the app's clock and try the newest code")
    conn.execute("UPDATE admin_totp SET enabled_at = CURRENT_TIMESTAMP, last_step = ? WHERE admin_user_id = ?",
                 (step, admin_user_id))
    codes = _replace_recovery_codes(conn, admin_user_id)
    conn.commit()
    return codes


def _replace_recovery_codes(conn, admin_user_id):
    conn.execute("DELETE FROM admin_recovery_codes WHERE admin_user_id = ?", (admin_user_id,))
    codes = [_new_recovery_code() for _ in range(RECOVERY_CODE_COUNT)]
    for code in codes:
        conn.execute("INSERT INTO admin_recovery_codes (admin_user_id, code_hash) VALUES (?, ?)",
                     (admin_user_id, _hash_token(_normalize_recovery_code(code))))
    return codes


def check_second_factor(conn, admin_user_id, code):
    """
    Whether `code` is a current authenticator code (not used before) or
    an unused recovery code for this admin; using either uses it up.

    Returns:
        str or None: `"totp"` or `"recovery"`, or `None` for a wrong code.
    """
    from . import totp
    row = conn.execute("SELECT secret, last_step FROM admin_totp WHERE admin_user_id = ? AND enabled_at IS NOT NULL "
                       "FOR UPDATE", (admin_user_id,)).fetchone()
    if row is None:
        conn.rollback()
        return None
    step = totp.verify(row["secret"], code, last_step=row["last_step"])
    if step is not None:
        conn.execute("UPDATE admin_totp SET last_step = ? WHERE admin_user_id = ?", (step, admin_user_id))
        conn.commit()
        return "totp"
    normalized = _normalize_recovery_code(code)
    if normalized:
        cur = conn.execute(
            "UPDATE admin_recovery_codes SET used_at = CURRENT_TIMESTAMP "
            "WHERE admin_user_id = ? AND code_hash = ? AND used_at IS NULL",
            (admin_user_id, _hash_token(normalized)))
        if cur.rowcount:
            conn.commit()
            return "recovery"
    conn.rollback()
    return None


def disable_totp(conn, admin_user_id):
    """Turns two-factor sign-in off for one admin (and drops its recovery
    codes). Returns whether it was on or being set up."""
    cur = conn.execute("DELETE FROM admin_totp WHERE admin_user_id = ?", (admin_user_id,))
    conn.execute("DELETE FROM admin_recovery_codes WHERE admin_user_id = ?", (admin_user_id,))
    conn.commit()
    return bool(cur.rowcount)


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


LOGIN_FAILURE_ACTIONS = ("login.failed", "login.locked", "password.failed", "totp.failed")
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
