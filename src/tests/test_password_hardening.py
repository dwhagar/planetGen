# tests/test_password_hardening.py

"""
The password checks added after the lockouts: wrong current passwords on
change-credentials count like failed logins (SEC.23), the trusted-device
cookie (SEC.22), the common-password blocklist and site-word rule
(SEC.24), and the hash settings with re-hash on login (SEC.25, then Argon2id
and the length bound: SEC.32).
"""

import pytest

from planetgen.api import loginguard
from planetgen.api.auth import DEVICE_COOKIE_NAME
from planetgen.db import store
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


# --- SEC.25, SEC.32: hashing ----------------------------------------------------


@pytest.mark.real_password_hashing
def test_new_hashes_are_argon2id_at_the_owasp_minimum():
    hashed = adminAuth.hash_password("violet-orbit-ledger-91")
    assert hashed.startswith("$argon2id$")
    assert "m=19456,t=2,p=1" in hashed
    assert adminAuth.hash_parameters() == (19456, 2, 1)
    assert not adminAuth.needs_rehash(hashed)
    assert adminAuth.verify_password("violet-orbit-ledger-91", hashed)
    assert not adminAuth.verify_password("violet-orbit-ledger-92", hashed)


@pytest.mark.real_password_hashing
def test_verify_sends_a_hash_by_its_prefix():
    from argon2 import PasswordHasher, Type
    from werkzeug.security import generate_password_hash
    password = "violet-orbit-ledger-91"
    assert adminAuth.verify_password(password, generate_password_hash(password, method="pbkdf2:sha256:1000"))
    assert adminAuth.verify_password(password, generate_password_hash(password, method="scrypt"))
    assert not adminAuth.verify_password("another-long-password", generate_password_hash(password, method="scrypt"))
    # Only Argon2id is ours: Argon2i and anything that is no hash at all never match.
    other_variant = PasswordHasher(type=Type.I, memory_cost=8, time_cost=1, parallelism=1).hash(password)
    for stored in (other_variant, "", None, "plain", "$argon2id$garbage", "pbkdf2:sha256:1000$only-two"):
        assert not adminAuth.verify_password(password, stored)


@pytest.mark.real_password_hashing
def test_needs_rehash_for_old_schemes_and_weaker_argon2_parameters(monkeypatch):
    from argon2 import PasswordHasher
    from werkzeug.security import generate_password_hash
    assert adminAuth.needs_rehash("")
    assert adminAuth.needs_rehash(None)
    assert adminAuth.needs_rehash("scrypt:32768:8:1$salt$abc")
    assert adminAuth.needs_rehash(generate_password_hash("x", method="pbkdf2:sha256:600000"))
    assert adminAuth.needs_rehash("$argon2id$not-a-hash")
    weaker_memory = PasswordHasher(memory_cost=8192, time_cost=2, parallelism=1).hash("x")
    weaker_passes = PasswordHasher(memory_cost=19456, time_cost=1, parallelism=1).hash("x")
    assert adminAuth.needs_rehash(weaker_memory)
    assert adminAuth.needs_rehash(weaker_passes)
    current = adminAuth.hash_password("x")
    assert not adminAuth.needs_rehash(current)
    # The settings name the cost: raising it makes today's hashes weaker.
    monkeypatch.setattr(adminAuth, "hash_parameters", lambda: (32768, 3, 1))
    assert adminAuth.needs_rehash(current)
    assert not adminAuth.needs_rehash(adminAuth.hash_password("x"))


@pytest.mark.real_password_hashing
def test_the_cost_comes_from_the_password_hash_settings(monkeypatch):
    from planetgen.util import settings
    monkeypatch.setenv("PLANETGEN_PASSWORD_HASH_MEMORY_KIB", "32768")
    monkeypatch.setenv("PLANETGEN_PASSWORD_HASH_TIME_COST", "3")
    monkeypatch.setenv("PLANETGEN_PASSWORD_HASH_PARALLELISM", "2")
    settings.reset_cache()
    try:
        assert adminAuth.hash_parameters() == (32768, 3, 2)
        assert "m=32768,t=3,p=2" in adminAuth.hash_password("violet-orbit-ledger-91")
    finally:
        monkeypatch.undo()
        settings.reset_cache()


def test_hash_parameter_settings_have_the_owasp_defaults_and_floors():
    from planetgen.util import settings
    section = settings.build({}, {}).password_hash
    assert (section.memory_kib, section.time_cost, section.parallelism) == (19456, 2, 1)
    for name, bad in (("memory_kib", 4), ("time_cost", 0), ("parallelism", 0)):
        with pytest.raises(settings.SettingsError):
            settings.build({"password_hash": {name: bad}}, {})


def test_a_password_is_at_most_1024_characters():
    assert adminAuth.MAX_PASSWORD_LENGTH == 1024
    longest = ("violet-orbit-ledger-91-" * 50)[:1024]
    adminAuth.validate_password_policy(longest, username="boss")
    with pytest.raises(adminAuth.AuthError, match="at most 1024"):
        adminAuth.validate_password_policy(longest + "x", username="boss")
    assert adminAuth.verify_password(longest, adminAuth.hash_password(longest))
    with pytest.raises(ValueError, match="at most 1024"):
        adminAuth.hash_password(longest + "x")


def test_an_oversize_password_is_refused_without_hashing(monkeypatch):
    hashed = adminAuth.hash_password("violet-orbit-ledger-91")

    def explode():
        raise AssertionError("an oversize password must not reach the hasher")
    monkeypatch.setattr(adminAuth, "_hasher", explode)
    assert not adminAuth.verify_password("x" * 1025, hashed)
    assert not adminAuth.verify_password(b"bytes", hashed)


@pytest.fixture
def control_conn(mysql_config):
    conn = store.get_control_connection(mysql_config, ensure_schema=True)
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
    assert new.startswith("$argon2id$")
    assert adminAuth.authenticate(control_conn, "oldhash", "violet-orbit-ledger-91")["username"] == "oldhash"


@pytest.mark.real_password_hashing
def test_login_rehashes_a_weaker_argon2_hash(control_conn):
    from argon2 import PasswordHasher
    weak = PasswordHasher(memory_cost=8192, time_cost=1, parallelism=1).hash("violet-orbit-ledger-91")
    cur = control_conn.execute(
        "INSERT INTO admin_users (username, password_hash, must_change_credentials) VALUES (?, ?, 0)",
        ("weakhash", weak))
    control_conn.commit()
    adminAuth.authenticate(control_conn, "weakhash", "violet-orbit-ledger-91")
    new = control_conn.execute("SELECT password_hash FROM admin_users WHERE id = ?",
                               (cur.lastrowid,)).fetchone()["password_hash"]
    assert "m=19456,t=2,p=1" in new
    assert not adminAuth.needs_rehash(new)


def test_login_with_an_oversize_password_fails_like_a_wrong_one(control_conn):
    cur = control_conn.execute(
        "INSERT INTO admin_users (username, password_hash, must_change_credentials) VALUES (?, ?, 0)",
        ("longpass", adminAuth.hash_password("violet-orbit-ledger-91")))
    control_conn.commit()
    with pytest.raises(adminAuth.AuthError, match="invalid username or password"):
        adminAuth.authenticate(control_conn, "longpass", "violet-orbit-ledger-91" + "x" * 1100)
    with pytest.raises(adminAuth.AuthError, match="invalid username or password"):
        adminAuth.authenticate(control_conn, "nobody", "x" * 1100)
    assert cur.lastrowid


@pytest.mark.real_password_hashing
def test_dummy_hash_uses_the_same_cost():
    assert "m=19456,t=2,p=1" in adminAuth._get_dummy_password_hash()


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
def real_app(mysql_config, redis_server, key_prefix):
    from planetgen.web.app import create_app
    from planetgen.api.config import Config

    class RealConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "test-secret"
        RATELIMIT_STORAGE_URI = redis_server  # the lockouts count in Redis (SEC.30)
        RATELIMIT_KEY_PREFIX = key_prefix

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


def _session_headers(response):
    return [h for h in response.headers.getlist("Set-Cookie") if "pg_admin_session" in h]


def test_secure_session_cookie_carries_the_host_prefix(real_app):
    """SEC.32: the prefix needs Secure, Path=/ and no Domain, and its removal at logout is Secure too."""
    real_app.config["SESSION_COOKIE_SECURE"] = True
    https = {"base_url": "https://localhost"}
    client = real_app.test_client()
    response = client.post("/api/auth/login", json={"username": "admin", "password": real_app.first_password},
                           environ_base=_address(), **https)
    assert response.status_code == 200
    (header,) = _session_headers(response)
    assert header.startswith("__Host-pg_admin_session=")
    assert "; Secure" in header and "; Path=/" in header and "Domain" not in header and "HttpOnly" in header
    assert client.get("/api/auth/me", **https).status_code == 200
    out = client.post("/api/auth/logout", **https)
    (cleared,) = _session_headers(out)
    assert cleared.startswith("__Host-pg_admin_session=;") and "; Secure" in cleared and "Max-Age=0" in cleared
    assert client.get("/api/auth/me", **https).status_code == 401


def test_insecure_development_session_cookie_keeps_the_plain_name(real_app):
    client = real_app.test_client()
    response = _login(client, "admin", real_app.first_password)
    (header,) = _session_headers(response)
    assert header.startswith("pg_admin_session=") and "Secure" not in header
    assert client.get("/api/auth/me").status_code == 200


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
    conn = store.get_control_connection(real_app.config["CONTROL_MYSQL_CONFIG"])
    try:
        assert adminAuth.device_username(conn, old_device) is None
        assert adminAuth.device_username(conn, new_device) == "admin"
    finally:
        conn.close()


def test_command_line_forgets_devices(mysql_config, capsys, monkeypatch):
    from planetgen.cli import forget_devices
    adminAuth.bootstrap_control_schema(mysql_config)
    conn = store.get_control_connection(mysql_config)
    try:
        admin_id = conn.execute("SELECT id FROM admin_users WHERE username = 'admin'").fetchone()["id"]
        raw = adminAuth.create_device(conn, admin_id)
    finally:
        conn.close()
    args = ["--mysql-host", mysql_config.host, "--mysql-port", str(mysql_config.port),
            "--mysql-user", mysql_config.user, "--mysql-password", mysql_config.password,
            "--mysql-database", mysql_config.database]
    monkeypatch.setattr(store, "configured_control_database", lambda: mysql_config.database)
    assert forget_devices.main(["admin"] + args) == 0
    assert "Revoked 1 trusted device of admin." in capsys.readouterr().out
    conn = store.get_control_connection(mysql_config)
    try:
        assert adminAuth.device_username(conn, raw) is None
    finally:
        conn.close()
    with pytest.raises(SystemExit):
        forget_devices.main(["nobody"] + args)
