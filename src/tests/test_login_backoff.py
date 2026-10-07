# tests/test_login_backoff.py

"""
Login lockouts (`planetgen/admin/throttle.py`, `api/loginguard.py`):
per address (SEC.1: 3 failures, 5 minutes doubling to a day) and per
username (SEC.21, formerly SEC.17's in-memory backoff: 10 free failures,
1 s doubling to 15 minutes). First the rules with a fake clock, then the
login route and page against a fake `adminAuth.authenticate` (no
database: the in-memory fallback), then the real `login_throttle` table.
"""

import time

import pytest

from api import auth, loginguard
from planetgen.admin import auth as adminAuth, throttle
from planetgen.admin.throttle import IP_POLICY, SCOPE_IP, SCOPE_USER, USER_POLICY

FREE_FAILURES = USER_POLICY.free_failures


@pytest.fixture
def store():
    return throttle.MemoryStore()


def _fail(store, scope, subject, now):
    return throttle.record_failure(store, scope, subject, now=now)


# --- Per username -------------------------------------------------------------

def test_free_failures_then_doubling_locks(store):
    now = 1000.0
    for _ in range(FREE_FAILURES):
        assert _fail(store, SCOPE_USER, "admin", now) == 0
        assert throttle.check(store, SCOPE_USER, "admin", now) == 0
    assert _fail(store, SCOPE_USER, "admin", now) == 1
    assert throttle.check(store, SCOPE_USER, "admin", now) == 1
    now += 1
    assert throttle.check(store, SCOPE_USER, "admin", now) == 0
    assert _fail(store, SCOPE_USER, "admin", now) == 2
    now += 2
    assert _fail(store, SCOPE_USER, "admin", now) == 4


def test_username_lock_is_capped(store):
    for _ in range(FREE_FAILURES + 60):
        seconds = _fail(store, SCOPE_USER, "admin", 1000.0)
    assert seconds == USER_POLICY.max_lock_seconds
    assert throttle.check(store, SCOPE_USER, "admin", 1000.0) == int(USER_POLICY.max_lock_seconds)


def test_usernames_are_counted_case_and_space_insensitively():
    assert throttle.normalize_username(" Admin ") == throttle.normalize_username("ADMIN") == "admin"


def test_success_clears_a_usernames_count(store):
    for _ in range(FREE_FAILURES):
        _fail(store, SCOPE_USER, "admin", 1000.0)
    throttle.record_success(store, SCOPE_USER, "admin", now=1000.0)
    assert store.get(SCOPE_USER, "admin") is None
    assert _fail(store, SCOPE_USER, "admin", 1000.0) == 0


def test_an_old_failure_is_forgotten(store):
    for _ in range(FREE_FAILURES):
        _fail(store, SCOPE_USER, "admin", 1000.0)
    assert _fail(store, SCOPE_USER, "admin", 1000.0 + USER_POLICY.forget_after_seconds) == 0


def test_memory_store_is_bounded_and_keeps_locked_ones(store, monkeypatch):
    monkeypatch.setattr(throttle, "MAX_MEMORY_ENTRIES", 20)
    now = time.time()
    for _ in range(FREE_FAILURES + 1):
        _fail(store, SCOPE_USER, "target", now)
    for i in range(100):
        _fail(store, SCOPE_USER, f"junk-{i}", now)
    assert len(store) <= 20
    assert throttle.check(store, SCOPE_USER, "target", now) >= 1


# --- Per address --------------------------------------------------------------

def test_three_failures_lock_an_address_for_five_minutes_doubling_to_a_day(store):
    now = 1000.0
    locks = []
    for _ in range(10):
        for attempt in range(IP_POLICY.free_failures):
            seconds = _fail(store, SCOPE_IP, "203.0.113.5", now)
            if attempt < IP_POLICY.free_failures - 1:
                assert seconds == 0
        locks.append(seconds)
        assert throttle.check(store, SCOPE_IP, "203.0.113.5", now) == int(seconds)
        now += seconds  # the lock runs out; the count starts over
    assert locks == [300, 600, 1200, 2400, 4800, 9600, 19200, 38400, 76800, 86400]
    # A day-long lock is itself a day without a lockout, so the level has
    # halved by the time it runs out.
    for _ in range(3):
        seconds = _fail(store, SCOPE_IP, "203.0.113.5", now)
    assert seconds == 9600


def test_success_keeps_an_addresss_doubling_level(store):
    now = 1000.0
    for _ in range(3):
        _fail(store, SCOPE_IP, "203.0.113.5", now)
    now += 300
    _fail(store, SCOPE_IP, "203.0.113.5", now)
    throttle.record_success(store, SCOPE_IP, "203.0.113.5", now=now)
    assert store.get(SCOPE_IP, "203.0.113.5")["failures"] == 0
    _fail(store, SCOPE_IP, "203.0.113.5", now)
    _fail(store, SCOPE_IP, "203.0.113.5", now)
    assert _fail(store, SCOPE_IP, "203.0.113.5", now) == 600  # not back to 300


def test_doubling_level_halves_each_day_without_a_lockout(store):
    now = 1000.0
    for _level in range(4):  # 300, 600, 1200, 2400: level 4 now
        for _ in range(3):
            seconds = _fail(store, SCOPE_IP, "198.51.100.1", now)
        now += seconds
    now += 2 * 86400  # two halvings: level 1
    for _ in range(3):
        seconds = _fail(store, SCOPE_IP, "198.51.100.1", now)
    assert seconds == 600


def test_address_subjects():
    assert throttle.ip_subject("203.0.113.5") == "203.0.113.5"
    assert throttle.ip_subject("2001:db8:1:2:3:4:5:6") == "2001:db8:1:2::/64"
    assert throttle.ip_subject("2001:db8:1:2::ffff") == "2001:db8:1:2::/64"
    assert throttle.ip_subject("::ffff:203.0.113.5") == "203.0.113.5"
    for exempt in ("127.0.0.1", "127.8.9.10", "::1", "nonsense", "", None):
        assert throttle.ip_subject(exempt) is None
    networks, bad = throttle.parse_allowlist(["198.51.100.0/24", "2001:db8::1", "junk"])
    assert bad == ["junk"]
    assert throttle.ip_subject("198.51.100.77", networks) is None
    assert throttle.ip_subject("2001:db8::1", networks) is None
    assert throttle.ip_subject("198.51.101.1", networks) == "198.51.101.1"
    assert len(throttle.parse_allowlist("10.0.0.1, 10.0.0.2")[0]) == 2


def test_private_addresses():
    assert throttle.is_private_address("10.1.2.3")
    assert throttle.is_private_address("192.168.0.4")
    assert not throttle.is_private_address("127.0.0.1")
    assert not throttle.is_private_address("8.8.8.8")
    assert not throttle.is_private_address("bogus")


def test_wait_text():
    assert loginguard.wait_text(1) == "1 second"
    assert loginguard.wait_text(45) == "45 seconds"
    assert loginguard.wait_text(300) == "5 minutes"
    assert loginguard.wait_text(86400) == "24 hours"


# --- The login route and page (no database: the in-memory fallback) ---------

@pytest.fixture
def app(monkeypatch):
    from api.app import create_app
    from api.config import Config

    class _Config(Config):
        TESTING = True
        SECRET_KEY = "test"

    monkeypatch.setattr(auth, "get_control_db", lambda: object())
    monkeypatch.setattr(loginguard, "get_control_db", lambda: object())
    good = {"id": 1, "username": "admin", "must_change_credentials": 0}

    def authenticate(_conn, username, password):
        if username == "admin" and password == "right":
            return good
        raise auth.adminAuth.AuthError("invalid username or password")

    monkeypatch.setattr(auth.adminAuth, "authenticate", authenticate)
    monkeypatch.setattr(auth.adminAuth, "create_session", lambda _conn, _id: "token")
    application = create_app(_Config)
    application.testing = True
    return application


_addresses = iter(range(1, 10 ** 6))


def _fresh_address():
    """A new client address per request, so neither the per-IP rate limit
    nor the per-address lockout answers first: these tests are about the
    per-username count."""
    n = next(_addresses)
    return {"REMOTE_ADDR": f"10.{n // 65536 % 256}.{n // 256 % 256}.{n % 256}"}


def _login(client, username, password, environ=None):
    return client.post("/api/auth/login", json={"username": username, "password": password},
                       environ_base=environ or _fresh_address())


def test_route_locks_a_username_from_any_address(app):
    client = app.test_client()
    for _ in range(FREE_FAILURES):
        assert _login(client, "admin", "wrong").status_code == 401
    assert _login(client, "admin", "wrong").status_code == 401
    # Locked now: even the right password is refused, unchecked.
    response = _login(client, "admin", "right")
    assert response.status_code == 429
    # Flask-Limiter's headers raise Retry-After to its own window's reset
    # when that is later; never below the lock.
    assert int(response.headers["Retry-After"]) >= 1
    assert response.get_json()["retry_after"] == 1
    assert response.get_json()["scope"] == "user"
    # Other usernames are unaffected.
    assert _login(client, "someone", "wrong").status_code == 401


def test_route_locks_unknown_usernames_the_same_way(app):
    client = app.test_client()
    for _ in range(FREE_FAILURES + 1):
        assert _login(client, "no-such-admin", "wrong").status_code == 401
    assert _login(client, "no-such-admin", "wrong").status_code == 429


def test_route_success_resets_the_count(app):
    client = app.test_client()
    for _ in range(FREE_FAILURES):
        _login(client, "admin", "wrong")
    assert _login(client, "admin", "right").status_code == 200
    assert _login(client, "admin", "wrong").status_code == 401
    assert _login(client, "admin", "wrong").status_code == 401


def test_route_locks_an_address_after_three_failures(app):
    client = app.test_client()
    where = {"REMOTE_ADDR": "203.0.113.50"}
    for name in ("a", "b", "c"):
        assert _login(client, name, "wrong", where).status_code == 401
    response = _login(client, "admin", "right", where)
    assert response.status_code == 429
    assert response.get_json()["scope"] == "ip"
    assert response.get_json()["retry_after"] == 300
    assert "try again in 5 minutes" in response.get_json()["error"]
    # Another address is unaffected, and loopback is never locked.
    assert _login(client, "admin", "right").status_code == 200
    for _ in range(5):
        _login(client, "x", "wrong", {"REMOTE_ADDR": "127.0.0.1"})
    assert _login(client, "admin", "right", {"REMOTE_ADDR": "127.0.0.1"}).status_code == 200


def test_allowlisted_addresses_are_never_locked(app):
    app.config["LOGIN_ALLOWLIST"] = ["203.0.113.0/24"]
    client = app.test_client()
    where = {"REMOTE_ADDR": "203.0.113.60"}
    for name in ("a", "b", "c", "d"):
        assert _login(client, name, "wrong", where).status_code == 401
    assert _login(client, "admin", "right", where).status_code == 200


def test_login_page_names_the_wait(app):
    from web import csrf
    client = app.test_client()
    client.get("/login")
    nonce = client.get_cookie(csrf.COOKIE_NAME).value
    with app.test_request_context():
        token = csrf._sign(nonce, "")
    form = {"username": "admin", "password": "wrong", "csrf_token": token}
    for _ in range(FREE_FAILURES + 1):
        assert client.post("/login", data=form, environ_base=_fresh_address()).status_code == 401
    response = client.post("/login", data=form, environ_base=_fresh_address())
    assert response.status_code == 429
    assert "Too many failed logins for this username. Try again in 1 second." in response.get_data(as_text=True)


# --- The login_throttle table -------------------------------------------------

@pytest.fixture
def control(mysql_config):
    adminAuth.bootstrap_control_schema(mysql_config)
    conn = adminAuth._db.get_control_connection(mysql_config)
    yield conn
    conn.close()


def test_db_store_applies_the_same_rules(control):
    store = throttle.DbStore(control)
    now = time.time()
    assert _fail(store, SCOPE_IP, "203.0.113.5", now) == 0
    assert _fail(store, SCOPE_IP, "203.0.113.5", now) == 0
    assert _fail(store, SCOPE_IP, "203.0.113.5", now) == 300
    assert throttle.check(store, SCOPE_IP, "203.0.113.5", now + 10) == 290
    for _ in range(FREE_FAILURES + 1):
        _fail(store, SCOPE_USER, "admin", now)
    assert sorted((r["scope"], r["subject"]) for r in throttle.locked_subjects(store, now=now)) == [
        (SCOPE_IP, "203.0.113.5"), (SCOPE_USER, "admin")]
    throttle.record_success(store, SCOPE_USER, "admin", now=now)
    assert store.get(SCOPE_USER, "admin") is None
    assert store.lift(SCOPE_IP, "203.0.113.5") == 1
    assert throttle.check(store, SCOPE_IP, "203.0.113.5", now) == 0


def test_control_schema_has_login_throttle(control):
    row = control.execute("SELECT MAX(version) AS v FROM control_schema_migrations").fetchone()
    assert row["v"] == adminAuth._db.CONTROL_SCHEMA_VERSION >= 2


def test_older_control_schema_gets_the_table(mysql_config):
    adminAuth.bootstrap_control_schema(mysql_config)
    conn = adminAuth._db.get_control_connection(mysql_config)
    try:
        conn.execute("DROP TABLE login_throttle")
        conn.execute("DELETE FROM control_schema_migrations")
        conn.execute("INSERT INTO control_schema_migrations (version) VALUES (1)")
        conn.commit()
    finally:
        conn.close()
    adminAuth.bootstrap_control_schema(mysql_config)
    conn = adminAuth._db.get_control_connection(mysql_config)
    try:
        assert throttle.DbStore(conn).get(SCOPE_IP, "203.0.113.5") is None
        versions = [r["version"] for r in conn.execute("SELECT version FROM control_schema_migrations").fetchall()]
        assert sorted(versions) == [1, adminAuth._db.CONTROL_SCHEMA_VERSION]
    finally:
        conn.close()


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
    adminAuth._db.get_connection(mysql_config).close()
    application = create_app(RealConfig)
    application.testing = True
    application.first_password = password
    return application


def _admin_client(real_app):
    admin = real_app.test_client()
    assert _login(admin, "admin", real_app.first_password, {"REMOTE_ADDR": "127.0.0.1"}).status_code == 200
    assert admin.post("/api/auth/change-credentials", json={
        "current_password": real_app.first_password, "new_username": "admin",
        "new_password": "a-new-long-password-1"}).status_code == 200
    return admin


def test_lockouts_shared_through_the_database_and_lifted_by_an_admin(real_app):
    client = real_app.test_client()
    attacker = {"REMOTE_ADDR": "93.184.216.34"}
    for _ in range(3):
        assert _login(client, "admin", "wrong", attacker).status_code == 401
    assert _login(client, "admin", real_app.first_password, attacker).status_code == 429
    # The in-memory fallback wasn't used: the count is in the table.
    assert len(loginguard.memory_store) == 0

    admin = _admin_client(real_app)
    body = admin.get("/api/admin/lockouts").get_json()
    assert [(i["scope"], i["subject"]) for i in body["items"]] == [("ip", "93.184.216.34")]
    assert body["items"][0]["retry_after"] > 290 and body["proxy_warning"] is False
    assert admin.post("/api/admin/lockouts/lift", json={"scope": "ip", "subject": "93.184.216.34"}).get_json() == {
        "lifted": 1}
    assert admin.get("/api/admin/lockouts").get_json()["items"] == []
    assert admin.post("/api/admin/lockouts/lift", json={"scope": "nope", "subject": "x"}).status_code == 400
    assert admin.post("/api/admin/lockouts/lift", json={"all": True}).status_code == 200
    assert _login(client, "admin", "a-new-long-password-1", attacker).status_code == 200


def test_private_lockout_warns_about_a_proxy(real_app):
    client = real_app.test_client()
    for _ in range(3):
        _login(client, "x", "wrong", {"REMOTE_ADDR": "10.0.0.2"})
    admin = _admin_client(real_app)
    assert admin.get("/api/admin/lockouts").get_json()["proxy_warning"] is True
    real_app.config["PROXY_FIX"] = {"x_for": 1, "x_proto": 1, "x_host": 0}
    assert admin.get("/api/admin/lockouts").get_json()["proxy_warning"] is False


def test_command_line_lists_and_lifts(mysql_config, capsys, monkeypatch):
    import loginLockouts
    adminAuth.bootstrap_control_schema(mysql_config)
    conn = adminAuth._db.get_control_connection(mysql_config)
    try:
        store = throttle.DbStore(conn)
        for _ in range(3):
            throttle.record_failure(store, SCOPE_IP, "2001:db8:1:2::/64")
    finally:
        conn.close()
    args = ["--mysql-host", mysql_config.host, "--mysql-port", str(mysql_config.port),
            "--mysql-user", mysql_config.user, "--mysql-password", mysql_config.password,
            "--mysql-database", mysql_config.database]
    monkeypatch.setattr(adminAuth._db, "configured_control_database", lambda: mysql_config.database)
    assert loginLockouts.main(args) == 0
    assert "2001:db8:1:2::/64" in capsys.readouterr().out
    assert loginLockouts.main(["--ip", "2001:db8:1:2::99"] + args) == 0
    assert "Lifted 1 lockout (ip:2001:db8:1:2::/64)" in capsys.readouterr().out
    assert loginLockouts.main(args) == 0
    assert "Nothing is locked." in capsys.readouterr().out
