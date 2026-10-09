# tests/test_two_factor.py

"""
Two-factor sign-in for admins (SEC.26): the TOTP codes themselves, the
setup and recovery codes in the control database, the two-step API
login, the web pages, and the command-line reset.
"""

import base64
import re
import time

import pyotp
import pytest

from planetgen.api import loginguard
from planetgen.db import store
from planetgen.admin import auth as adminAuth

# --- Codes (RFC 6238) -----------------------------------------------------------

RFC_SECRET = base64.b32encode(b"12345678901234567890").decode()


@pytest.mark.parametrize("when, code", [(59, "287082"), (1111111109, "081804"), (1234567890, "005924"),
                                        (2000000000, "279037")])
def test_rfc_6238_vectors(when, code):
    assert pyotp.TOTP(RFC_SECRET).at(when) == code


def _code(step):
    return pyotp.TOTP(RFC_SECRET).at(step * 30)


def test_verify_allows_one_step_of_drift_and_no_replay():
    verify = adminAuth.verify_totp
    now = 1_000_000_000
    step = now // 30
    assert verify(RFC_SECRET, _code(step), 0, now) == step
    assert verify(RFC_SECRET, _code(step - 1), 0, now) == step - 1
    assert verify(RFC_SECRET, _code(step + 1), 0, now) == step + 1
    assert verify(RFC_SECRET, _code(step - 2), 0, now) is None
    # A step already used (or older) never works again.
    assert verify(RFC_SECRET, _code(step), step, now) is None
    assert verify(RFC_SECRET, "abcdef", 0, now) is None
    assert verify(RFC_SECRET, None, 0, now) is None
    # Spaces and dashes in what was typed are ignored.
    code = _code(step)
    assert verify(RFC_SECRET, f"{code[:3]} {code[3:]}", 0, now) == step
    assert verify(RFC_SECRET, f"{code[:3]}-{code[3:]}", 0, now) == step


@pytest.mark.parametrize("offset", range(-3, 4))
def test_verify_is_exact_around_a_step_boundary(offset):
    """TEST.72: the window is the same just before, on and just after a step boundary."""
    boundary = 1_000_000_020  # a multiple of 30
    for now in (boundary - 1, boundary, boundary + 1):
        step = now // 30
        expected = step + offset if abs(offset) <= 1 else None
        assert adminAuth.verify_totp(RFC_SECRET, _code(step + offset), 0, now) == expected


def test_new_secret_and_enrolment():
    secret = adminAuth.new_totp_secret()
    assert re.fullmatch(r"[A-Z2-7]{32}", secret)
    uri, svg = adminAuth.totp_enrolment(secret, "boss admin")
    assert uri.startswith("otpauth://totp/planetGen:boss%20admin?secret=" + secret)
    assert "issuer=planetGen" in uri
    assert svg.startswith("<svg") and "<path" in svg and "<script" not in svg and "style=" not in svg


# --- The control database ----------------------------------------------------------


@pytest.fixture
def control_conn(mysql_config):
    conn = store.get_control_connection(mysql_config, ensure_schema=True)
    try:
        yield conn
    finally:
        conn.close()


@pytest.fixture
def admin_id(control_conn):
    cur = control_conn.execute(
        "INSERT INTO admin_users (username, password_hash, must_change_credentials) VALUES (?, ?, 0)",
        ("twofa", adminAuth.hash_password("violet-orbit-ledger-91")))
    control_conn.commit()
    return cur.lastrowid


class _Clock:
    """`time` as the auth code sees it: `time()` is held at the middle of the
    TOTP step the test started in, everything else is the real module."""

    def __init__(self, real, now):
        self._real, self.now = real, now

    def time(self):
        return self.now

    def __getattr__(self, name):
        return getattr(self._real, name)


_clock = _Clock(time, time.time())


@pytest.fixture(autouse=True)
def _one_totp_step(monkeypatch):
    """Holds the clock the TOTP checks read inside one 30 s step, so a test
    run slowly (a loaded machine) cannot cross a step boundary between
    making a code and checking it, which turned "used once" codes into
    fresh ones (and the other way round)."""
    _clock.now = (int(time.time() // 30) * 30) + 15
    monkeypatch.setattr(adminAuth, "time", _clock)


def _now_code(secret, offset=0):
    return pyotp.TOTP(secret).at((int(_clock.now // 30) + offset) * 30)


def test_setup_confirm_and_check(control_conn, admin_id):
    assert adminAuth.totp_status(control_conn, admin_id) == {"enabled": False, "recovery_codes_left": 0}
    with pytest.raises(adminAuth.AuthError):
        adminAuth.confirm_totp_setup(control_conn, admin_id, "123456")
    secret = adminAuth.begin_totp_setup(control_conn, admin_id)
    assert not adminAuth.totp_enabled(control_conn, admin_id)
    with pytest.raises(adminAuth.AuthError, match="doesn't match"):
        adminAuth.confirm_totp_setup(control_conn, admin_id, "000000" if _now_code(secret) != "000000" else "111111")
    codes = adminAuth.confirm_totp_setup(control_conn, admin_id, _now_code(secret, -1))
    assert len(codes) == 10 and len(set(codes)) == 10
    assert all(re.fullmatch(r"[a-z2-9]{5}-[a-z2-9]{5}", c) for c in codes)
    assert adminAuth.totp_status(control_conn, admin_id) == {"enabled": True, "recovery_codes_left": 10}
    with pytest.raises(adminAuth.AuthError, match="already on"):
        adminAuth.begin_totp_setup(control_conn, admin_id)
    stored = [r["code_hash"] for r in control_conn.execute("SELECT code_hash FROM admin_recovery_codes").fetchall()]
    assert not any(c.replace("-", "") in stored for c in codes)

    # The confirming step is used up; the next one works once.
    assert adminAuth.check_second_factor(control_conn, admin_id, _now_code(secret, -1)) is None
    assert adminAuth.check_second_factor(control_conn, admin_id, _now_code(secret)) == "totp"
    assert adminAuth.check_second_factor(control_conn, admin_id, _now_code(secret)) is None
    # A recovery code works once, typed any which way.
    assert adminAuth.check_second_factor(control_conn, admin_id, " " + codes[0].upper().replace("-", "") + " ") \
        == "recovery"
    assert adminAuth.check_second_factor(control_conn, admin_id, codes[0]) is None
    assert adminAuth.totp_status(control_conn, admin_id)["recovery_codes_left"] == 9
    assert adminAuth.check_second_factor(control_conn, admin_id, "") is None

    assert adminAuth.disable_totp(control_conn, admin_id) is True
    assert adminAuth.totp_status(control_conn, admin_id) == {"enabled": False, "recovery_codes_left": 0}
    assert adminAuth.check_second_factor(control_conn, admin_id, codes[1]) is None


# --- The API ---------------------------------------------------------------------------


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
    client = application.test_client()
    assert _login(client, "admin", password).status_code == 200
    assert client.post("/api/auth/change-credentials", json={
        "current_password": password, "new_username": "admin",
        "new_password": "violet-orbit-ledger-91"}).status_code == 200
    application.admin_client = client
    return application


_counter = iter(range(1, 10 ** 6))


def _address():
    n = next(_counter)
    return {"REMOTE_ADDR": f"10.88.{n // 256 % 256}.{n % 256}"}


def _login(client, username, password, environ=None):
    return client.post("/api/auth/login", json={"username": username, "password": password},
                       environ_base=environ or _address())


def _turn_on(client):
    setup = client.post("/api/auth/totp/setup", json={"current_password": "violet-orbit-ledger-91"})
    assert setup.status_code == 200
    body = setup.get_json()
    assert body["uri"].startswith("otpauth://totp/") and body["qr_svg"].startswith("<svg")
    confirmed = client.post("/api/auth/totp/confirm", json={"code": _now_code(body["secret"], -1)})
    assert confirmed.status_code == 200
    return body["secret"], confirmed.get_json()["recovery_codes"]


def test_two_step_login(real_app):
    admin = real_app.admin_client
    assert admin.get("/api/auth/totp").get_json() == {"enabled": False, "recovery_codes_left": 0}
    assert admin.post("/api/auth/totp/setup", json={"current_password": "wrong-password"}).status_code == 400
    secret, codes = _turn_on(admin)
    assert admin.get("/api/auth/totp").get_json() == {"enabled": True, "recovery_codes_left": 10}

    browser = real_app.test_client()
    first = _login(browser, "admin", "violet-orbit-ledger-91")
    assert first.status_code == 200
    body = first.get_json()
    assert body["totp_required"] is True and "username" not in body
    assert not [h for h in first.headers.getlist("Set-Cookie") if h.startswith("pg_admin_session=")]
    pending = body["pending"]

    wrong = browser.post("/api/auth/login/totp", json={"pending": pending, "code": "000000"
                                                       if _now_code(secret) != "000000" else "111111"},
                         environ_base=_address())
    assert wrong.status_code == 401
    assert browser.post("/api/auth/login/totp", json={"pending": pending + "x", "code": "123456"}).status_code == 401
    ok = browser.post("/api/auth/login/totp", json={"pending": pending, "code": _now_code(secret)},
                      environ_base=_address())
    assert ok.status_code == 200 and ok.get_json()["username"] == "admin"
    assert browser.get("/api/auth/me").status_code == 200

    # A recovery code also works, once.
    other = real_app.test_client()
    pending = _login(other, "admin", "violet-orbit-ledger-91").get_json()["pending"]
    assert other.post("/api/auth/login/totp", json={"pending": pending, "code": codes[3]},
                      environ_base=_address()).status_code == 200
    pending = _login(other, "admin", "violet-orbit-ledger-91").get_json()["pending"]
    assert other.post("/api/auth/login/totp", json={"pending": pending, "code": codes[3]},
                      environ_base=_address()).status_code == 401


def test_wrong_codes_lock_the_address(real_app):
    secret, _codes = _turn_on(real_app.admin_client)
    attacker = real_app.test_client()
    where = {"REMOTE_ADDR": "93.184.216.40"}
    pending = _login(attacker, "admin", "violet-orbit-ledger-91", where).get_json()["pending"]
    bad = "000000" if _now_code(secret) != "000000" else "111111"
    for _ in range(3):
        assert attacker.post("/api/auth/login/totp", json={"pending": pending, "code": bad},
                             environ_base=where).status_code == 401
    response = attacker.post("/api/auth/login/totp", json={"pending": pending, "code": _now_code(secret)},
                             environ_base=where)
    assert response.status_code == 429 and response.get_json()["scope"] == "ip"
    assert len(loginguard.memory_store) == 0


def test_a_password_change_voids_pending_logins(real_app):
    secret, _codes = _turn_on(real_app.admin_client)
    browser = real_app.test_client()
    pending = _login(browser, "admin", "violet-orbit-ledger-91").get_json()["pending"]
    assert real_app.admin_client.post("/api/auth/change-credentials", json={
        "current_password": "violet-orbit-ledger-91", "new_username": "admin",
        "new_password": "amber-comet-harbor-42"}).status_code == 200
    response = browser.post("/api/auth/login/totp", json={"pending": pending, "code": _now_code(secret, 1)},
                            environ_base=_address())
    assert response.status_code == 401 and "password" in response.get_json()["error"]


def test_turning_it_off_needs_password_and_code(real_app):
    admin = real_app.admin_client
    secret, codes = _turn_on(admin)
    assert admin.post("/api/auth/totp/disable", json={"current_password": "violet-orbit-ledger-91",
                                                      "code": "not-a-code"}).status_code == 400
    assert admin.post("/api/auth/totp/disable", json={"current_password": "nope-nope-nope",
                                                      "code": codes[0]}).status_code == 400
    assert admin.post("/api/auth/totp/disable", json={"current_password": "violet-orbit-ledger-91",
                                                      "code": codes[0]}).get_json() == {"enabled": False}
    first = _login(real_app.test_client(), "admin", "violet-orbit-ledger-91")
    assert first.get_json()["username"] == "admin"


def test_api_keys_skip_the_second_factor(real_app):
    admin = real_app.admin_client
    _turn_on(admin)
    key = admin.post("/api/auth/api-keys", json={"label": "script"}).get_json()["key"]
    caller = real_app.test_client()
    assert caller.get("/api/auth/me", headers={"Authorization": f"Bearer {key}"}).status_code == 200


def test_command_line_reset(mysql_config, real_app, capsys, monkeypatch):
    from planetgen.cli import lockouts
    _turn_on(real_app.admin_client)
    args = ["--mysql-host", mysql_config.host, "--mysql-port", str(mysql_config.port),
            "--mysql-user", mysql_config.user, "--mysql-password", mysql_config.password,
            "--mysql-database", mysql_config.database]
    monkeypatch.setattr(store, "configured_control_database", lambda: mysql_config.database)
    assert lockouts.main(["--reset-two-factor", "admin"] + args) == 0
    assert "Two-factor sign-in for admin is off." in capsys.readouterr().out
    assert _login(real_app.test_client(), "admin", "violet-orbit-ledger-91").get_json()["username"] == "admin"


# --- The web pages -------------------------------------------------------------------


def _csrf(app, client):
    from planetgen.api.authz import SESSION_COOKIE_NAME
    from planetgen.web import csrf
    nonce = "n" * 43
    client.set_cookie(csrf.COOKIE_NAME, nonce)
    session = client.get_cookie(SESSION_COOKIE_NAME)
    with app.app_context():
        return {csrf.FIELD_NAME: csrf._sign(nonce, session.value if session else "")}


def test_web_account_turns_it_on_and_login_asks_for_the_code(real_app):
    from html import unescape
    admin = real_app.admin_client
    page = admin.get("/account").get_data(as_text=True)
    assert 'id="two-factor"' in page and "Set up" in page

    page = admin.post("/account/two-factor", data={"action": "setup", "current_password": "wrong-one",
                                                   **_csrf(real_app, admin)})
    assert page.status_code == 400 and "Current password is incorrect." in page.get_data(as_text=True)
    page = admin.post("/account/two-factor", data={"action": "setup", "current_password": "violet-orbit-ledger-91",
                                                   **_csrf(real_app, admin)}).get_data(as_text=True)
    assert '<div class="qr-code"><svg' in page
    secret = re.search(r'name="secret" value="([A-Z2-7]+)"', page).group(1)
    page = admin.post("/account/two-factor", data={"action": "confirm", "code": "abc", "secret": secret,
                                                   **_csrf(real_app, admin)}).get_data(as_text=True)
    assert "doesn" in page and '<div class="qr-code"><svg' in page
    page = admin.post("/account/two-factor", data={"action": "confirm", "code": _now_code(secret, -1),
                                                   "secret": secret, **_csrf(real_app, admin)}).get_data(as_text=True)
    codes = re.findall(r"<li><code>([a-z2-9]{5}-[a-z2-9]{5})</code></li>", page)
    assert len(codes) == 10 and "Two-factor sign-in is on." in page
    assert "Recovery codes left: 10" in admin.get("/account").get_data(as_text=True)

    browser = real_app.test_client()
    page = browser.post("/login", data={"username": "admin", "password": "violet-orbit-ledger-91",
                                        **_csrf(real_app, browser)}, environ_base=_address())
    assert page.status_code == 200
    html = page.get_data(as_text=True)
    assert 'autocomplete="one-time-code"' in html
    pending = unescape(re.search(r'name="pending" value="([^"]+)"', html).group(1))
    assert not browser.get_cookie("pg_admin_session")
    wrong = browser.post("/login/code", data={"pending": pending, "code": "nope", "username": "admin",
                                              **_csrf(real_app, browser)}, environ_base=_address())
    assert wrong.status_code == 401 and "That code isn" in wrong.get_data(as_text=True)
    ok = browser.post("/login/code", data={"pending": pending, "code": codes[0], "username": "admin",
                                           "next": "/admin", **_csrf(real_app, browser)}, environ_base=_address())
    assert ok.status_code == 303 and ok.headers["Location"].endswith("/admin")
    assert browser.get_cookie("pg_admin_session")

    page = admin.post("/account/two-factor", data={"action": "disable", "current_password": "violet-orbit-ledger-91",
                                                   "code": codes[1], **_csrf(real_app, admin)}).get_data(as_text=True)
    assert "Two-factor sign-in is off." in page
