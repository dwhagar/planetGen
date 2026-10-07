# tests/test_password_hardening.py

"""
The password checks added after the lockouts: wrong current passwords on
change-credentials count like failed logins (SEC.23), the trusted-device
cookie (SEC.22), the common-password blocklist and site-word rule
(SEC.24), and the hash settings with re-hash on login (SEC.25).
"""

import pytest

from api import loginguard
from api.auth import DEVICE_COOKIE_NAME
from stellarObjects import _db
from planetgen.admin import auth as adminAuth

# --- SEC.24: the policy -------------------------------------------------------


def test_blocklist_is_bundled_and_case_insensitive():
    assert len(adminAuth._load_common_passwords()) > 40000
    assert adminAuth.is_common_password("qwertyuiopasdfgh")
    assert adminAuth.is_common_password("QwertyuiopASDFGH")
    assert not adminAuth.is_common_password("violet-orbit-ledger-91")


@pytest.mark.parametrize("password", ["qwertyuiopasdfgh", "QWERTYUIOPASDFGH"])
def test_policy_refuses_common_passwords(password):
    with pytest.raises(adminAuth.AuthError, match="common or breached"):
        adminAuth.validate_password_policy(password, username="boss")


@pytest.mark.parametrize("password, username", [
    ("planetgenadmin1", "boss"),
    ("Password-Admin!", "boss"),
    ("captainkirk-123", "captainkirk"),
    ("PlanetGen 2026 !", "x"),
])
def test_policy_refuses_site_words_with_a_few_extras(password, username):
    with pytest.raises(adminAuth.AuthError, match="too close"):
        adminAuth.validate_password_policy(password, username=username)


@pytest.mark.parametrize("password", ["violet-orbit-ledger-91", "planetgen-is-my-hobby-now", "a-strong-test-password"])
def test_policy_accepts_long_uncommon_passwords(password):
    adminAuth.validate_password_policy(password, username="boss")


def test_missing_blocklist_file_does_not_break_the_policy(monkeypatch):
    monkeypatch.setattr(adminAuth, "COMMON_PASSWORDS_PATH", "/nonexistent/list.gz")
    monkeypatch.setattr(adminAuth, "_common_passwords", None)
    adminAuth.validate_password_policy("qwertyuiopasdfgh", username="boss")
    monkeypatch.setattr(adminAuth, "_common_passwords", None)


# --- SEC.25: hashing ------------------------------------------------------------


@pytest.mark.real_password_hashing
def test_new_hashes_use_pbkdf2_600000():
    hashed = adminAuth.hash_password("violet-orbit-ledger-91")
    assert hashed.startswith("pbkdf2:sha256:600000$")
    assert not adminAuth.needs_rehash(hashed)
    assert adminAuth.needs_rehash("scrypt:32768:8:1$salt$abc")
    assert adminAuth.needs_rehash("")


@pytest.fixture
def control_conn(mysql_config):
    conn = _db.get_control_connection(mysql_config, ensure_schema=True)
    try:
        yield conn
    finally:
        conn.close()


@pytest.mark.real_password_hashing
def test_login_rehashes_an_older_hash(control_conn):
    from werkzeug.security import generate_password_hash
    old = generate_password_hash("violet-orbit-ledger-91", method="scrypt")
    cur = control_conn.execute(
        "INSERT INTO admin_users (username, password_hash, must_change_credentials) VALUES (?, ?, 0)",
        ("oldhash", old))
    control_conn.commit()
    with pytest.raises(adminAuth.AuthError):
        adminAuth.authenticate(control_conn, "oldhash", "wrong-password-here")
    assert control_conn.execute("SELECT password_hash FROM admin_users WHERE id = ?",
                                (cur.lastrowid,)).fetchone()["password_hash"] == old
    adminAuth.authenticate(control_conn, "oldhash", "violet-orbit-ledger-91")
    new = control_conn.execute("SELECT password_hash FROM admin_users WHERE id = ?",
                               (cur.lastrowid,)).fetchone()["password_hash"]
    assert new.startswith("pbkdf2:sha256:600000$")
    assert adminAuth.authenticate(control_conn, "oldhash", "violet-orbit-ledger-91")["username"] == "oldhash"


@pytest.mark.real_password_hashing
def test_dummy_hash_uses_the_same_method():
    assert adminAuth._get_dummy_password_hash().startswith(adminAuth.PASSWORD_HASH_METHOD + "$")


# --- SEC.22: trusted devices (database) ----------------------------------------


def test_device_tokens_belong_to_one_admin_and_are_revoked_on_change(control_conn):
    cur = control_conn.execute(
        "INSERT INTO admin_users (username, password_hash, must_change_credentials) VALUES (?, ?, 0)",
        ("devowner", adminAuth.hash_password("violet-orbit-ledger-91")))
    control_conn.commit()
    admin_id = cur.lastrowid
    raw = adminAuth.create_device(control_conn, admin_id)
    assert adminAuth.device_username(control_conn, raw) == "devowner"
    assert adminAuth.device_username(control_conn, raw + "x") is None
    assert adminAuth.device_username(control_conn, None) is None
    row = control_conn.execute("SELECT token_hash FROM admin_devices WHERE admin_user_id = ?", (admin_id,)).fetchone()
    assert raw not in row["token_hash"]

    control_conn.execute("UPDATE admin_devices SET expires_at = DATE_SUB(CURRENT_TIMESTAMP, INTERVAL 1 DAY)")
    control_conn.commit()
    assert adminAuth.device_username(control_conn, raw) is None

    raw = adminAuth.create_device(control_conn, admin_id)
    # The expired one was deleted when the new one was made.
    assert control_conn.execute("SELECT COUNT(*) AS n FROM admin_devices").fetchone()["n"] == 1
    adminAuth.change_credentials(control_conn, admin_id, "violet-orbit-ledger-91", "devowner", "amber-comet-harbor-42")
    assert adminAuth.device_username(control_conn, raw) is None


# --- The routes, against a real control database -------------------------------


@pytest.fixture
def real_app(mysql_config):
    from api.app import create_app
    from api.config import Config

    class RealConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "test-secret"

    _username, password = adminAuth.bootstrap_control_schema(mysql_config)
    application = create_app(RealConfig)
    application.testing = True
    application.first_password = password
    return application


_counter = iter(range(1, 10 ** 6))


def _address():
    n = next(_counter)
    return {"REMOTE_ADDR": f"10.77.{n // 256 % 256}.{n % 256}"}


def _login(client, username, password, environ=None):
    return client.post("/api/auth/login", json={"username": username, "password": password},
                       environ_base=environ or _address())


def test_login_sets_a_device_cookie_that_skips_the_username_lock(real_app):
    browser = real_app.test_client()
    response = _login(browser, "admin", real_app.first_password)
    assert response.status_code == 200
    cookies = [h for h in response.headers.getlist("Set-Cookie") if h.startswith(DEVICE_COOKIE_NAME + "=")]
    assert len(cookies) == 1
    assert "HttpOnly" in cookies[0] and "SameSite=Strict" in cookies[0] and "Path=/" in cookies[0]

    # Someone else fails on purpose until the username is locked.
    attacker = real_app.test_client()
    for _ in range(11):
        _login(attacker, "admin", "wrong", _address())
    assert _login(attacker, "admin", real_app.first_password).status_code == 429
    # The browser with the device cookie still gets in, and its own
    # mistakes aren't counted against the username.
    assert _login(browser, "admin", "wrong").status_code == 401
    assert _login(browser, "admin", real_app.first_password).status_code == 200
    # But not from an address that is itself locked.
    locked = {"REMOTE_ADDR": "93.184.216.35"}
    for _ in range(3):
        _login(attacker, "nobody", "wrong", locked)
    assert _login(browser, "admin", real_app.first_password, locked).status_code == 429


def test_a_device_cookie_for_another_username_does_not_help(real_app):
    browser = real_app.test_client()
    assert _login(browser, "admin", real_app.first_password).status_code == 200
    for _ in range(11):
        _login(browser, "someone-else", "wrong")
    assert _login(browser, "someone-else", "wrong").status_code == 429


def test_wrong_current_password_counts_and_locks(real_app):
    """SEC.23: three wrong current passwords from one address lock it,
    and the right one is then refused unchecked."""
    client = real_app.test_client()
    where = {"REMOTE_ADDR": "93.184.216.36"}
    # Logged in from elsewhere; the change attempts come from `where`.
    assert _login(client, "admin", real_app.first_password).status_code == 200
    for _ in range(3):
        response = client.post("/api/auth/change-credentials", environ_base=where, json={
            "current_password": "wrong-password-xx", "new_username": "admin",
            "new_password": "violet-orbit-ledger-91"})
        assert response.status_code == 400
    response = client.post("/api/auth/change-credentials", environ_base=where, json={
        "current_password": real_app.first_password, "new_username": "admin",
        "new_password": "violet-orbit-ledger-91"})
    assert response.status_code == 429
    assert response.get_json()["scope"] == "ip"
    assert len(loginguard.memory_store) == 0


def test_change_credentials_refuses_a_common_password_and_reissues_the_device(real_app):
    client = real_app.test_client()
    assert _login(client, "admin", real_app.first_password).status_code == 200
    old_device = client.get_cookie(DEVICE_COOKIE_NAME).value
    response = client.post("/api/auth/change-credentials", json={
        "current_password": real_app.first_password, "new_username": "admin",
        "new_password": "qwertyuiopasdfgh"})
    assert response.status_code == 400
    assert "common or breached" in response.get_json()["error"]
    response = client.post("/api/auth/change-credentials", json={
        "current_password": real_app.first_password, "new_username": "admin",
        "new_password": "violet-orbit-ledger-91"})
    assert response.status_code == 200
    new_device = client.get_cookie(DEVICE_COOKIE_NAME).value
    assert new_device != old_device
    conn = _db.get_control_connection(real_app.config["CONTROL_MYSQL_CONFIG"])
    try:
        assert adminAuth.device_username(conn, old_device) is None
        assert adminAuth.device_username(conn, new_device) == "admin"
    finally:
        conn.close()


def test_command_line_forgets_devices(mysql_config, capsys, monkeypatch):
    import loginLockouts
    adminAuth.bootstrap_control_schema(mysql_config)
    conn = _db.get_control_connection(mysql_config)
    try:
        admin_id = conn.execute("SELECT id FROM admin_users WHERE username = 'admin'").fetchone()["id"]
        raw = adminAuth.create_device(conn, admin_id)
    finally:
        conn.close()
    args = ["--mysql-host", mysql_config.host, "--mysql-port", str(mysql_config.port),
            "--mysql-user", mysql_config.user, "--mysql-password", mysql_config.password,
            "--mysql-database", mysql_config.database]
    monkeypatch.setattr(_db, "configured_control_database", lambda: mysql_config.database)
    assert loginLockouts.main(["--forget-devices", "admin"] + args) == 0
    assert "Revoked 1 trusted device of admin." in capsys.readouterr().out
    conn = _db.get_control_connection(mysql_config)
    try:
        assert adminAuth.device_username(conn, raw) is None
    finally:
        conn.close()
    with pytest.raises(SystemExit):
        loginLockouts.main(["--forget-devices", "nobody"] + args)
