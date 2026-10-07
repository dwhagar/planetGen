# tests/test_api_auth_sweep.py

"""
Who may call what, across the whole API (TEST.43 to TEST.46):

- TEST.43: a sweep generated from `app.url_map`. Every API write route
  (POST/PATCH/PUT/DELETE, apart from the two sign-in steps) is behind
  `authz.require_admin`; every route behind it gives 401 to an anonymous
  caller, a garbage `Bearer` key and a revoked key; every one that needs
  fresh credentials gives 403 to an admin still on the seeded ones (and
  every admin page sends that admin to `/account`). A new route is swept
  without editing this file.
- TEST.44: what an API key may do. It can read (`/me`, its keys, the
  two-factor status), make content changes and revoke a key, but not
  manage the account: making keys, changing credentials, setting up or
  turning off two-factor sign-in and logging out need a browser session
  (403 with a key), so a leaked key can't mint another and outlive its
  own revocation.
- TEST.45: more than one admin. Admin B can't revoke admin A's key;
  lifting another admin's lockout is in the audit log under the admin
  who lifted it; two admins editing the same system both succeed, last
  write wins, each audited under its own name, and an edit to a system
  the other just deleted is a 404.
- TEST.46: trusted-device and two-factor edge cases: expired, tampered
  and other-admin device cookies don't skip a username lock; turning
  two-factor off forgets every trusted device (the browser that did it
  gets a new one); a code used at the API's second step can't be used
  again at `/login/code`; a pending login past its time is refused.

Skipped, not failed, without a reachable MySQL test server.
"""

import re

import pytest

from api import auth as auth_routes
from api.app import create_app
from api.config import Config
from api.limiter import limiter
from stellarObjects import _db, adminAuth
from planetgen.generation.config import SystemConfig
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.system import StarSystem

from tests.test_fuzz_web_routes import csrf_pair, session_of
from tests.test_two_factor import _now_code

PASSWORD_A = "violet-orbit-ledger-91"
PASSWORD_B = "amber-comet-harbor-42"
PASSWORD_NEW = "quiet-meadow-lantern-58"
SECRET = "test-secret"

PUBLIC_WRITES = {"/api/auth/login", "/api/auth/login/totp"}
"""The only API write routes anyone may call without signing in."""

_WRITE_METHODS = {"POST", "PATCH", "PUT", "DELETE"}


class _Admins:
    """The app and three admins: `a` and `b` past the credential change,
    `fresh_pending` still on seeded credentials."""


@pytest.fixture
def admins(mysql_config):
    class TestConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = SECRET
        # The sweep sends hundreds of requests from one address. The
        # limiter is one object for the whole process and keeps this
        # setting for every later app, so the fixture turns it back on.
        RATELIMIT_ENABLED = False

    _username, first_password = adminAuth.bootstrap_control_schema(mysql_config)
    _db.get_connection(mysql_config).close()  # the content schema
    app = create_app(TestConfig)
    app.testing = True

    state = _Admins()
    state.app = app
    state.config = mysql_config
    state.a = app.test_client()
    assert _login(state.a, "admin", first_password).status_code == 200
    assert state.a.post("/api/auth/change-credentials", json={
        "current_password": first_password, "new_username": "alice", "new_password": PASSWORD_A,
    }).status_code == 200

    conn = _db.get_control_connection(mysql_config)
    try:
        conn.execute("INSERT INTO admin_users (username, password_hash, must_change_credentials) VALUES (?, ?, 0)",
                     ("bob", adminAuth.hash_password(PASSWORD_B)))
        conn.execute("INSERT INTO admin_users (username, password_hash, must_change_credentials) VALUES (?, ?, 1)",
                     ("newcomer", adminAuth.hash_password(PASSWORD_B)))
        conn.commit()
    finally:
        conn.close()
    state.b = app.test_client()
    assert _login(state.b, "bob", PASSWORD_B).status_code == 200
    state.fresh_pending = app.test_client()
    response = _login(state.fresh_pending, "newcomer", PASSWORD_B)
    assert response.get_json()["must_change_credentials"] is True
    try:
        yield state
    finally:
        limiter.enabled = True


_counter = iter(range(1, 10 ** 6))


def _address():
    n = next(_counter)
    return {"REMOTE_ADDR": f"10.66.{n // 256 % 256}.{n % 256}"}


def _login(client, username, password, environ=None):
    return client.post("/api/auth/login", json={"username": username, "password": password},
                       environ_base=environ or _address())


def _bearer(key):
    return {"Authorization": f"Bearer {key}"}


def _new_key(client, label="script"):
    response = client.post("/api/auth/api-keys", json={"label": label})
    assert response.status_code == 201
    return response.get_json()


def _control(state):
    return _db.get_control_connection(state.config)


# --- Building the sweep from the URL map --------------------------------------

_ARGUMENT = re.compile(r"<(?:(\w+):)?(\w+)>")


def _concrete(rule):
    """`rule` with every `<converter:name>` replaced by a value it
    accepts (ids that don't exist; a real phenomenon type)."""
    def value(match):
        converter, name = match.group(1), match.group(2)
        if converter == "int":
            return "999999"
        if name == "phenomenon_type":
            return "nebula"
        return "x"
    return _ARGUMENT.sub(value, rule)


def _api_routes(app):
    """`(method, path, admin_required)` for every API route and method."""
    routes = []
    for rule in app.url_map.iter_rules():
        if not rule.rule.startswith("/api/"):
            continue
        view = app.view_functions[rule.endpoint]
        for method in sorted(rule.methods - {"HEAD", "OPTIONS"}):
            routes.append((method, _concrete(rule.rule), getattr(view, "admin_required", None)))
    return routes


def _sweep(app, required):
    return [(method, path, needs) for method, path, needs in _api_routes(app)
            if (needs is not None) == required]


# --- TEST.43: the sweep ---------------------------------------------------------

def test_every_api_write_route_requires_an_admin(admins):
    unguarded = [(method, path) for method, path, needs in _api_routes(admins.app)
                 if method in _WRITE_METHODS and needs is None and path not in PUBLIC_WRITES]
    assert unguarded == []
    # The sweep found the routes it's meant to cover.
    guarded = _sweep(admins.app, True)
    assert len(guarded) > 30
    assert ("POST", "/api/systems", {"fresh": True, "session_only": False}) in guarded


def test_public_routes_are_reads_or_sign_in(admins):
    for method, path, _needs in _sweep(admins.app, False):
        assert method == "GET" or path in PUBLIC_WRITES, (method, path)


@pytest.mark.parametrize("label,headers", [
    ("anonymous", {}),
    ("garbage key", _bearer("pg_not-a-real-key")),
    ("empty key", _bearer("")),
    ("session token as key", None),  # filled in below
])
def test_admin_routes_refuse_callers_without_valid_credentials(admins, label, headers):
    if headers is None:
        headers = _bearer(session_of(admins.a))
    client = admins.app.test_client(use_cookies=False)
    for method, path, _needs in _sweep(admins.app, True):
        response = client.open(path, method=method, headers=headers, json={})
        assert response.status_code == 401, (label, method, path, response.status_code)
        assert response.get_json() == {"error": "authentication required"}


def test_admin_routes_refuse_a_revoked_key(admins):
    created = _new_key(admins.a)
    client = admins.app.test_client(use_cookies=False)
    assert client.get("/api/auth/me", headers=_bearer(created["key"])).status_code == 200
    assert admins.a.delete(f"/api/auth/api-keys/{created['id']}").status_code == 200
    for method, path, _needs in _sweep(admins.app, True):
        response = client.open(path, method=method, headers=_bearer(created["key"]), json={})
        assert response.status_code == 401, (method, path, response.status_code)


def test_fresh_routes_refuse_an_admin_on_seeded_credentials(admins):
    fresh_routes = [(m, p) for m, p, needs in _sweep(admins.app, True) if needs["fresh"]]
    assert fresh_routes
    for method, path in fresh_routes:
        response = admins.fresh_pending.open(path, method=method, json={})
        assert response.status_code == 403, (method, path, response.status_code)
        assert "default credentials must be changed" in response.get_json()["error"]
    # The few that aren't fresh let that admin through (identity, the
    # credential change itself, logging out, a read of class options).
    assert admins.fresh_pending.get("/api/auth/me").status_code == 200


def _is_script_endpoint(app, rule, response):
    """The JSON endpoints under /admin a page's script polls (marked
    `json_only`) answer a plain 403 instead of a redirect."""
    if getattr(app.view_functions[rule.endpoint], "json_only", False):
        assert response.status_code == 403, rule.rule
        return True
    return False


def test_admin_pages_send_an_admin_on_seeded_credentials_to_account(admins):
    pages = [rule for rule in admins.app.url_map.iter_rules()
             if rule.rule.startswith("/admin") and "GET" in rule.methods]
    assert len(pages) >= 5
    for rule in pages:
        response = admins.fresh_pending.get(_concrete(rule.rule))
        if _is_script_endpoint(admins.app, rule, response):
            continue
        assert response.status_code == 302, (rule.rule, response.status_code)
        assert response.headers["Location"].startswith("/account"), rule.rule


def test_admin_pages_send_anonymous_visitors_to_login(admins):
    client = admins.app.test_client()
    for rule in admins.app.url_map.iter_rules():
        if rule.rule.startswith("/admin") and "GET" in rule.methods:
            response = client.get(_concrete(rule.rule))
            if _is_script_endpoint(admins.app, rule, response):
                continue
            assert response.status_code == 302, rule.rule
            assert response.headers["Location"].startswith("/login"), rule.rule


# --- TEST.44: what an API key may do ----------------------------------------------

def test_an_api_key_can_read_and_change_content(admins):
    key = _new_key(admins.a)["key"]
    caller = admins.app.test_client(use_cookies=False)
    assert caller.get("/api/auth/me", headers=_bearer(key)).get_json()["username"] == "alice"
    assert caller.get("/api/auth/api-keys", headers=_bearer(key)).status_code == 200
    assert caller.get("/api/auth/totp", headers=_bearer(key)).status_code == 200
    assert caller.get("/api/admin/stats", headers=_bearer(key)).status_code == 200
    created = caller.post("/api/sectors", headers=_bearer(key), json={"name": "Keyed Sector", "edge_ly": 10})
    assert created.status_code == 201


@pytest.mark.parametrize("path,body", [
    ("/api/auth/api-keys", {"label": "minted by a key"}),
    ("/api/auth/change-credentials", {"current_password": PASSWORD_A, "new_username": "mallory",
                                      "new_password": PASSWORD_NEW}),
    ("/api/auth/totp/setup", {"current_password": PASSWORD_A}),
    ("/api/auth/totp/confirm", {"code": "123456"}),
    ("/api/auth/totp/disable", {"current_password": PASSWORD_A, "code": "123456"}),
    ("/api/auth/logout", {}),
])
def test_an_api_key_cannot_manage_the_account(admins, path, body):
    key = _new_key(admins.a)["key"]
    caller = admins.app.test_client(use_cookies=False)
    response = caller.post(path, headers=_bearer(key), json=body)
    assert response.status_code == 403
    assert "API key" in response.get_json()["error"]
    # Nothing changed: still one key, same name and password, no
    # two-factor setup started, the key still works.
    assert len(admins.a.get("/api/auth/api-keys").get_json()["items"]) == 1
    assert caller.get("/api/auth/me", headers=_bearer(key)).get_json()["username"] == "alice"
    assert _login(admins.app.test_client(), "alice", PASSWORD_A).status_code == 200
    conn = _control(admins)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM admin_totp").fetchone()["n"] == 0
    finally:
        conn.close()


def test_the_same_things_work_from_a_browser_session(admins):
    assert admins.a.post("/api/auth/totp/setup", json={"current_password": PASSWORD_A}).status_code == 200
    assert admins.a.post("/api/auth/totp/disable", json={"current_password": PASSWORD_A}).status_code == 200
    assert admins.a.post("/api/auth/change-credentials", json={
        "current_password": PASSWORD_A, "new_username": "alice", "new_password": PASSWORD_NEW,
    }).status_code == 200
    assert admins.a.post("/api/auth/logout").status_code == 200


def test_an_api_key_can_revoke_itself(admins):
    created = _new_key(admins.a)
    caller = admins.app.test_client(use_cookies=False)
    assert caller.delete(f"/api/auth/api-keys/{created['id']}", headers=_bearer(created["key"])).status_code == 200
    assert caller.get("/api/auth/me", headers=_bearer(created["key"])).status_code == 401


def test_a_key_with_a_session_cookie_is_still_a_key(admins):
    # The header wins over the cookie (authz._current_admin), so a
    # browser that also sends a key can't use the cookie to get round it.
    key = _new_key(admins.a)["key"]
    response = admins.a.post("/api/auth/api-keys", headers=_bearer(key), json={"label": "both"})
    assert response.status_code == 403


# --- TEST.45: more than one admin --------------------------------------------------

def test_admin_b_cannot_revoke_admin_a_key(admins):
    a_key = _new_key(admins.a)
    response = admins.b.delete(f"/api/auth/api-keys/{a_key['id']}")
    assert response.status_code == 404
    assert admins.app.test_client().get("/api/auth/me", headers=_bearer(a_key["key"])).status_code == 200
    # Nor see it.
    assert admins.b.get("/api/auth/api-keys").get_json()["items"] == []
    # Nor can A's key revoke B's.
    b_key = _new_key(admins.b)
    caller = admins.app.test_client(use_cookies=False)
    assert caller.delete(f"/api/auth/api-keys/{b_key['id']}", headers=_bearer(a_key["key"])).status_code == 404
    assert caller.get("/api/auth/me", headers=_bearer(b_key["key"])).get_json()["username"] == "bob"


def test_lifting_another_admins_lockout_is_audited_under_the_lifter(admins):
    attacker = admins.app.test_client()
    for _ in range(11):
        _login(attacker, "alice", "wrong-password-here")
    assert _login(admins.app.test_client(), "alice", PASSWORD_A).status_code == 429

    lifted = admins.b.post("/api/admin/lockouts/lift", json={"scope": "user", "subject": "Alice"})
    assert lifted.status_code == 200 and lifted.get_json()["lifted"] == 1
    assert _login(admins.app.test_client(), "alice", PASSWORD_A).status_code == 200

    conn = _control(admins)
    try:
        rows = conn.execute("SELECT * FROM admin_audit_log WHERE action = 'lockout.lift'").fetchall()
        bob_id = conn.execute("SELECT id FROM admin_users WHERE username = 'bob'").fetchone()["id"]
    finally:
        conn.close()
    assert len(rows) == 1
    assert rows[0]["admin_username"] == "bob" and rows[0]["admin_user_id"] == bob_id
    assert rows[0]["target"] == "user:alice"


def _seed_system(config):
    sector = SpaceSector("Shared Sector", edge_ly=10.0)
    system_config = SystemConfig()
    system_config.STAR_TYPE = "G2V"
    system_config.PLANETS = False
    system_config.BINARY_SYSTEM = False
    sector.add_system(StarSystem(system_config=system_config), position=(1.0, 1.0, 1.0),
                      system_config=system_config)
    sector_id = _db.save_sector(sector, config=config)
    conn = _db.get_connection(config)
    try:
        return conn.execute("SELECT id FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()["id"]
    finally:
        conn.close()


def test_two_admins_editing_the_same_system(admins):
    system_id = _seed_system(admins.config)
    assert admins.a.patch(f"/api/systems/{system_id}", json={"name": "Alpha Haven"}).status_code == 200
    assert admins.b.patch(f"/api/systems/{system_id}", json={"name": "Beta Haven"}).status_code == 200
    assert admins.a.get(f"/api/systems/{system_id}").get_json()["name"] == "Beta Haven"

    conn = _control(admins)
    try:
        rows = conn.execute("SELECT admin_username, target FROM admin_audit_log WHERE action LIKE ? ORDER BY id",
                            ("system.%",)).fetchall()
    finally:
        conn.close()
    assert [row["admin_username"] for row in rows] == ["alice", "bob"]
    assert {row["target"] for row in rows} == {f"system:{system_id}"}

    # B deletes it; A's edit, made from a page that still showed it, is a
    # clean 404.
    assert admins.b.delete(f"/api/systems/{system_id}").status_code == 200
    late = admins.a.patch(f"/api/systems/{system_id}", json={"name": "Gamma Haven"})
    assert late.status_code == 404
    assert "error" in late.get_json()


# --- TEST.46: trusted devices and two-factor sign-in ---------------------------------

_DEVICE = auth_routes.DEVICE_COOKIE_NAME


def _device_cookie(response):
    for header in response.headers.getlist("Set-Cookie"):
        if header.startswith(_DEVICE + "=") and not header.startswith(_DEVICE + "=;"):
            return header.split(";", 1)[0].split("=", 1)[1]
    return None


def _lock_username(app, username):
    attacker = app.test_client()
    for _ in range(11):
        _login(attacker, username, "wrong-password-here")
    assert _login(app.test_client(), username, PASSWORD_A).status_code == 429


def _login_with_device(app, username, password, device):
    browser = app.test_client()
    browser.set_cookie(_DEVICE, device, domain="localhost")
    return _login(browser, username, password)


def test_a_good_device_cookie_skips_the_username_lock(admins):
    device = _device_cookie(_login(admins.app.test_client(), "alice", PASSWORD_A))
    assert device
    _lock_username(admins.app, "alice")
    assert _login_with_device(admins.app, "alice", PASSWORD_A, device).status_code == 200


def test_an_expired_device_cookie_does_not(admins):
    device = _device_cookie(_login(admins.app.test_client(), "alice", PASSWORD_A))
    conn = _control(admins)
    try:
        conn.execute("UPDATE admin_devices SET expires_at = DATE_SUB(CURRENT_TIMESTAMP, INTERVAL 1 MINUTE)")
        conn.commit()
    finally:
        conn.close()
    _lock_username(admins.app, "alice")
    assert _login_with_device(admins.app, "alice", PASSWORD_A, device).status_code == 429


@pytest.mark.parametrize("tamper", [
    lambda raw: raw[:-1] + ("A" if raw[-1] != "A" else "B"),
    lambda raw: raw + "x",
    lambda raw: raw[1:],
    lambda raw: raw.upper(),
    lambda raw: "",
    lambda raw: "' OR '1'='1",
])
def test_a_tampered_device_cookie_does_not(admins, tamper):
    device = _device_cookie(_login(admins.app.test_client(), "alice", PASSWORD_A))
    _lock_username(admins.app, "alice")
    response = _login_with_device(admins.app, "alice", PASSWORD_A, tamper(device))
    assert response.status_code == 429


def test_another_admins_device_cookie_does_not(admins):
    bobs_device = _device_cookie(_login(admins.app.test_client(), "bob", PASSWORD_B))
    _lock_username(admins.app, "alice")
    assert _login_with_device(admins.app, "alice", PASSWORD_A, bobs_device).status_code == 429


def _turn_on_two_factor(client, password=PASSWORD_A):
    setup = client.post("/api/auth/totp/setup", json={"current_password": password})
    assert setup.status_code == 200
    secret = setup.get_json()["secret"]
    confirmed = client.post("/api/auth/totp/confirm", json={"code": _now_code(secret, -1)})
    assert confirmed.status_code == 200
    return secret, confirmed.get_json()["recovery_codes"]


def test_turning_two_factor_off_forgets_trusted_devices(admins):
    _secret, codes = _turn_on_two_factor(admins.a)
    # Another browser signed in with the second factor and is trusted.
    laptop = admins.app.test_client()
    pending = _login(laptop, "alice", PASSWORD_A).get_json()["pending"]
    signed_in = laptop.post("/api/auth/login/totp", json={"pending": pending, "code": codes[0]},
                            environ_base=_address())
    assert signed_in.status_code == 200
    laptop_device = _device_cookie(signed_in)
    conn = _control(admins)
    try:
        assert adminAuth.device_username(conn, laptop_device) == "alice"
    finally:
        conn.close()

    turned_off = admins.a.post("/api/auth/totp/disable", json={"current_password": PASSWORD_A, "code": codes[1]})
    assert turned_off.status_code == 200
    new_device = _device_cookie(turned_off)
    assert new_device and new_device != laptop_device

    conn = _control(admins)
    try:
        assert adminAuth.device_username(conn, laptop_device) is None
        assert adminAuth.device_username(conn, new_device) == "alice"
        audited = conn.execute("SELECT COUNT(*) AS n FROM admin_audit_log WHERE action = 'totp.disable'").fetchone()
    finally:
        conn.close()
    assert audited["n"] == 1
    _lock_username(admins.app, "alice")
    assert _login_with_device(admins.app, "alice", PASSWORD_A, laptop_device).status_code == 429
    assert _login_with_device(admins.app, "alice", PASSWORD_A, new_device).status_code == 200


def test_turning_two_factor_off_on_the_account_page_keeps_this_browser_trusted(admins):
    _secret, codes = _turn_on_two_factor(admins.a)
    nonce, token = csrf_pair(secret=SECRET, session=session_of(admins.a))
    from web import csrf
    admins.a.set_cookie(csrf.COOKIE_NAME, nonce, domain="localhost")
    response = admins.a.post("/account/two-factor", data={
        "action": "disable", "current_password": PASSWORD_A, "code": codes[0], csrf.FIELD_NAME: token,
    })
    assert response.status_code == 200
    assert "Two-factor sign-in is off." in response.get_data(as_text=True)
    device = _device_cookie(response)
    assert device
    conn = _control(admins)
    try:
        assert adminAuth.device_username(conn, device) == "alice"
    finally:
        conn.close()


def test_a_code_used_at_the_api_cannot_be_used_again_on_the_login_page(admins):
    secret, _codes = _turn_on_two_factor(admins.a)
    code = _now_code(secret)
    api_browser = admins.app.test_client()
    pending = _login(api_browser, "alice", PASSWORD_A).get_json()["pending"]
    assert api_browser.post("/api/auth/login/totp", json={"pending": pending, "code": code},
                            environ_base=_address()).status_code == 200

    from web import csrf
    page_browser = admins.app.test_client()
    nonce, token = csrf_pair(secret=SECRET)
    page_browser.set_cookie(csrf.COOKIE_NAME, nonce, domain="localhost")
    first = page_browser.post("/login", data={"username": "alice", "password": PASSWORD_A, csrf.FIELD_NAME: token},
                              environ_base=_address())
    assert first.status_code == 200
    page_pending = re.search(r'name="pending" value="([^"]+)"', first.get_data(as_text=True)).group(1)
    second = page_browser.post("/login/code", data={
        "pending": page_pending.replace("&#34;", '"').replace("&amp;", "&"), "code": code,
        "username": "alice", csrf.FIELD_NAME: token,
    }, environ_base=_address())
    assert second.status_code == 401
    assert "code isn" in second.get_data(as_text=True)
    assert page_browser.get("/api/auth/me").status_code == 401

    # And the other way round: a code used on the page is spent at the API.
    later = _now_code(secret, 1)
    page_pending = re.search(r'name="pending" value="([^"]+)"', page_browser.post(
        "/login", data={"username": "alice", "password": PASSWORD_A, csrf.FIELD_NAME: token},
        environ_base=_address()).get_data(as_text=True)).group(1)
    assert page_browser.post("/login/code", data={"pending": page_pending, "code": later, "username": "alice",
                                                  csrf.FIELD_NAME: token},
                             environ_base=_address()).status_code == 303
    pending = _login(api_browser, "alice", PASSWORD_A).get_json()["pending"]
    assert api_browser.post("/api/auth/login/totp", json={"pending": pending, "code": later},
                            environ_base=_address()).status_code == 401


def test_a_pending_login_that_expired_is_refused(admins, monkeypatch):
    secret, _codes = _turn_on_two_factor(admins.a)
    browser = admins.app.test_client()
    pending = _login(browser, "alice", PASSWORD_A).get_json()["pending"]
    # Any pending login is now older than the allowed age.
    monkeypatch.setattr(auth_routes, "PENDING_LOGIN_SECONDS", -1)
    response = browser.post("/api/auth/login/totp", json={"pending": pending, "code": _now_code(secret)},
                            environ_base=_address())
    assert response.status_code == 401
    assert "took too long" in response.get_json()["error"]
    assert browser.get("/api/auth/me").status_code == 401
    # The code wasn't spent by the refused attempt.
    monkeypatch.setattr(auth_routes, "PENDING_LOGIN_SECONDS", 300)
    pending = _login(browser, "alice", PASSWORD_A).get_json()["pending"]
    assert browser.post("/api/auth/login/totp", json={"pending": pending, "code": _now_code(secret)},
                        environ_base=_address()).status_code == 200


def test_a_pending_login_for_another_secret_key_is_refused(admins):
    secret, _codes = _turn_on_two_factor(admins.a)
    from itsdangerous import URLSafeTimedSerializer
    conn = _control(admins)
    try:
        row = conn.execute("SELECT id, password_hash FROM admin_users WHERE username = 'alice'").fetchone()
    finally:
        conn.close()
    forged = URLSafeTimedSerializer("some-other-secret", salt="planetgen-totp-login").dumps(
        {"id": row["id"], "h": auth_routes._hash_fingerprint(row["password_hash"])})
    response = admins.app.test_client().post("/api/auth/login/totp", json={"pending": forged,
                                                                           "code": _now_code(secret)},
                                             environ_base=_address())
    assert response.status_code == 401
    assert "sign in with your password first" in response.get_json()["error"]
