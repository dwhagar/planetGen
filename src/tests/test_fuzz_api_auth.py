# tests/test_fuzz_api_auth.py

"""
Brute-force/property-based tests for admin authentication: the JSON
login (`POST /api/auth/login`), the `/login` page, API keys
(`Authorization: Bearer ...`), the session cookie, the one-shot flash
cookie and credential changes.

Invariants: hostile usernames/passwords never 5xx; a wrong credential is
always refused and never starts a session; a garbage or tampered API key
or session cookie is never accepted, on the API (401) or on the pages
(302 to `/login`); the real ones keep working throughout.

Shares `test_fuzz_web_routes.py`'s seeded-database/app fixtures (a
separate database for this module) and its `check_response` invariants.
Skipped, not failed, without a reachable MySQL test server.
"""

import json

import pytest
from hypothesis import assume, example, given, settings
from hypothesis import strategies as st
from itsdangerous import URLSafeTimedSerializer

from api.authz import SESSION_COOKIE_NAME
from stellarObjects import _db, adminAuth

from web import csrf  # noqa: E402
from web.admin_pages import FLASH_COOKIE  # noqa: E402

from tests.fuzz_support import any_float, hostile_text, scaled
from tests.test_fuzz_web_routes import (  # noqa: F401 -- fixtures
    ADMIN_PASSWORD, ADMIN_USERNAME, SECRET_KEY, admin_client, app, check_response, fuzz_db, login_api, with_csrf,
)

_ME = "/api/auth/me"
_PROTECTED_API = [("GET", _ME), ("GET", "/api/auth/api-keys"), ("GET", "/api/admin/stats"),
                  ("POST", "/api/auth/api-keys"), ("POST", "/api/sectors"), ("POST", "/api/auth/logout")]
_PROTECTED_PAGES = ["/admin", "/admin/stats", "/account", "/admin/generate"]


def _session_cookies(response):
    """The non-empty session cookies a response sets."""
    return [h for h in response.headers.getlist("Set-Cookie")
            if h.startswith(SESSION_COOKIE_NAME + "=") and not h.startswith(SESSION_COOKIE_NAME + "=;")]


def _session_token(response):
    (header,) = _session_cookies(response)
    return header.split(";", 1)[0].split("=", 1)[1]


def raw_client(app):
    """A test client that sends the `Cookie` header exactly as given
    (Flask's default client replaces it with its own cookie jar)."""
    return app.test_client(use_cookies=False)


def _is_real_password(password):
    """The real password, or it followed by NULs: scrypt/PBKDF2-HMAC
    zero-pads the key, so werkzeug's `check_password_hash` treats
    `pw` and `pw\x00` as the same password (a property of HMAC, not a
    bypass -- the password itself is still required)."""
    return isinstance(password, str) and password.rstrip("\x00") == ADMIN_PASSWORD


def _header_safe(value):
    """What an HTTP header can carry: bytes, read back as latin-1 by the
    WSGI server; no CR/LF/NUL (a server rejects the request before the
    app sees it)."""
    raw = value.encode("utf-8").decode("latin-1")
    return "".join(ch for ch in raw if ch not in "\r\n\x00")


@pytest.fixture(scope="module")
def real_session(app):
    """A real, live session token for the admin (its own client, so the
    token-tampering tests below never log the shared admin_client out)."""
    return _session_token(login_api(app.test_client()))


@pytest.fixture(scope="module")
def real_api_key(admin_client):
    response = admin_client.post("/api/auth/api-keys", json={"label": "fuzz key"})
    assert response.status_code == 201
    return response.get_json()["key"]


def _assert_refused_everywhere(app, headers=None, cookie=None):
    """No protected API route or admin page accepts these credentials."""
    client = raw_client(app)
    headers = dict(headers or {})
    if cookie is not None:
        headers["Cookie"] = cookie
    for method, path in _PROTECTED_API:
        response = client.open(path, method=method, headers=headers, json={"label": "x", "name": "x", "edge_ly": 1})
        check_response(response, path)
        assert response.status_code == 401, f"{method} {path} accepted {headers!r} -> {response.status_code}"
    for path in _PROTECTED_PAGES:
        response = client.get(path, headers=headers)
        check_response(response, path)
        assert response.status_code == 302, f"{path} accepted {headers!r} -> {response.status_code}"


# ---------------------------------------------------------------------
# Login: hostile usernames and passwords
# ---------------------------------------------------------------------

json_scalar = st.one_of(hostile_text, st.integers(min_value=-(2 ** 70), max_value=2 ** 70), any_float,
                        st.booleans(), st.none(), st.lists(hostile_text, max_size=2),
                        st.dictionaries(hostile_text, hostile_text, max_size=2))


@settings(max_examples=scaled(60))
@given(username=hostile_text, password=hostile_text)
@example(username=ADMIN_USERNAME, password="")
@example(username=ADMIN_USERNAME, password=ADMIN_PASSWORD[:-1])
@example(username=ADMIN_USERNAME, password="\x00" + ADMIN_PASSWORD)
@example(username=ADMIN_USERNAME, password=ADMIN_PASSWORD.upper())
@example(username=ADMIN_USERNAME, password=" " + ADMIN_PASSWORD)
@example(username=ADMIN_USERNAME, password="' OR '1'='1")
@example(username="' OR '1'='1' -- ", password="' OR '1'='1")
@example(username="admin", password="password")  # the old published default (#39: now random, and rotated here)
@example(username=ADMIN_USERNAME + "' -- ", password="x")
@example(username="%", password="x")
@example(username="\x00", password="\x00")
@example(username="a" * 10000, password="b" * 10000)
def test_login_hostile_credentials_are_refused(app, username, password):
    assume(not _is_real_password(password))
    response = app.test_client().post("/api/auth/login", json={"username": username, "password": password})
    check_response(response, "/api/auth/login", sent=[username, password])
    assert response.status_code in (400, 401), (username, password, response.status_code)
    assert not _session_cookies(response)
    error = response.get_json()["error"]
    # One generic answer: never says whether the username exists.
    assert error in ("invalid username or password", "username and password are required"), error


@settings(max_examples=scaled(30))
@given(username=hostile_text, password=hostile_text)
@example(username="no-such-admin", password="x")
@example(username="admin", password="password")  # #39: the old default user no longer exists here
@example(username=ADMIN_USERNAME, password="wrong-password")
@example(username=ADMIN_USERNAME + "x", password=ADMIN_PASSWORD)
def test_login_costs_one_password_hash_check_whether_or_not_the_user_exists(app, username, password):
    """Security #44: every refused login runs exactly one password hash
    check -- a real user's hash, or a dummy one for an unknown username --
    so the response time doesn't reveal which usernames exist."""
    assume(not _is_real_password(password))
    checked = []
    real_verify = adminAuth.verify_password
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(adminAuth, "verify_password",
                   lambda pw, pw_hash: checked.append(pw_hash) or real_verify(pw, pw_hash))
        response = app.test_client().post("/api/auth/login", json={"username": username, "password": password})
    check_response(response, "/api/auth/login", sent=[username, password])
    assert response.status_code in (400, 401), response.status_code
    if response.status_code == 401:
        assert len(checked) == 1, checked
    else:  # refused before any lookup (empty username or password)
        assert checked == []


@settings(max_examples=scaled(60))
@given(username=json_scalar, password=json_scalar, extra=st.dictionaries(hostile_text, json_scalar, max_size=2))
@example(username=[ADMIN_USERNAME], password=ADMIN_PASSWORD, extra={})
@example(username=ADMIN_USERNAME, password=[ADMIN_PASSWORD], extra={})
@example(username={"$ne": None}, password={"$ne": None}, extra={})
@example(username=None, password=None, extra={})
@example(username=0, password=0, extra={})
@example(username=1, password=ADMIN_PASSWORD, extra={})  # A1
@example(username=1.5, password="x", extra={})  # A1
def test_login_wrongly_typed_fields_are_refused(app, username, password, extra):
    """`username`/`password` of any JSON type: never logs in unless both
    are the real strings, never a 5xx."""
    body = {**extra, "username": username, "password": password}
    response = app.test_client().post("/api/auth/login", json=body)
    check_response(response, "/api/auth/login")
    assert response.status_code in (400, 401)
    assert not _session_cookies(response)


# A1 (fixed): api/auth.py login() calls .strip() on a non-string username (int, list, dict) -> AttributeError -> 500
@pytest.mark.parametrize("username", [1, 1.5, [ADMIN_USERNAME], {"a": 1}, True])
def test_regression_login_non_string_username_is_a_400(app, username):
    response = app.test_client().post("/api/auth/login", json={"username": username, "password": ADMIN_PASSWORD})
    assert response.status_code == 400


@pytest.mark.parametrize("data,content_type", [
    (b"", "application/json"), (b"{", "application/json"), (b"[]", "application/json"), (b"null", "application/json"),
    (b'"admin"', "application/json"), (b"\xff\xfe\x00", "application/json"),
    (json.dumps({"username": ADMIN_USERNAME, "password": ADMIN_PASSWORD}).encode(), "text/plain"),
    (f"username={ADMIN_USERNAME}&password={ADMIN_PASSWORD}".encode(), "application/x-www-form-urlencoded"),
    (json.dumps({"username": ADMIN_USERNAME, "password": ADMIN_PASSWORD}).encode(), ""),
    (b'{"username": "' + ADMIN_USERNAME.encode() + b'", "password": NaN}', "application/json"),
    (b'{"username": "a", "password": "' + b"p" * 1_000_000 + b'"}', "application/json"),
])
def test_login_rejects_non_json_bodies(app, data, content_type):
    """Only a JSON object logs in -- a form-encoded or text/plain body
    (which a cross-site form *can* send) never does."""
    response = app.test_client().post("/api/auth/login", data=data, content_type=content_type)
    check_response(response, "/api/auth/login")
    assert response.status_code in (400, 401), response.status_code
    assert not _session_cookies(response)


def test_login_rejects_a_duplicate_key_trick(app):
    """`{"password": <real>, "password": <wrong>}` -- JSON's last key wins
    (Python's json), so the wrong one is checked; and the other order
    logs in. Never a 5xx either way."""
    client = app.test_client()
    raw = '{"username": "%s", "password": "%s", "password": "wrong"}' % (ADMIN_USERNAME, ADMIN_PASSWORD)
    assert client.post("/api/auth/login", data=raw, content_type="application/json").status_code == 401
    raw = '{"username": "%s", "password": "wrong", "password": "%s"}' % (ADMIN_USERNAME, ADMIN_PASSWORD)
    assert client.post("/api/auth/login", data=raw, content_type="application/json").status_code == 200


_PASSWORD_MUTATIONS = [
    ADMIN_PASSWORD.upper(), ADMIN_PASSWORD.capitalize(), ADMIN_PASSWORD[:-1], ADMIN_PASSWORD[1:],
    ADMIN_PASSWORD + " ", " " + ADMIN_PASSWORD, "\x00" + ADMIN_PASSWORD, ADMIN_PASSWORD + "\n",
    ADMIN_PASSWORD.replace("-", "‐"), ADMIN_PASSWORD.replace("a", "а"),  # hyphen / Cyrillic a lookalikes
    ADMIN_PASSWORD * 2, ADMIN_PASSWORD[:8], "", "password", "%", "*",
]


@pytest.mark.parametrize("password", _PASSWORD_MUTATIONS)
@pytest.mark.parametrize("username", [ADMIN_USERNAME, ADMIN_USERNAME.upper(), " " + ADMIN_USERNAME + " "])
def test_wrong_password_for_the_real_user_is_always_refused(app, username, password):
    response = app.test_client().post("/api/auth/login", json={"username": username, "password": password})
    check_response(response, "/api/auth/login")
    assert response.status_code in (400, 401)
    assert not _session_cookies(response)


@pytest.mark.parametrize("username,accepted", [
    (ADMIN_USERNAME, True),
    ("  " + ADMIN_USERNAME + "\t", True),  # stripped by the login route
    (ADMIN_USERNAME + "\x00", None),  # NUL is ignorable in utf8mb4_unicode_ci: may match
    (ADMIN_USERNAME.upper(), None),  # case-insensitive collation: may match
    (ADMIN_USERNAME.replace("-", "‐"), False),  # a different character, not a case/accent variant
    (ADMIN_USERNAME + "x", False),
    (ADMIN_USERNAME[:-1], False),
    ("%", False), ("_" * len(ADMIN_USERNAME), False), ("fuzz%", False),  # LIKE wildcards are not wildcards
])
def test_right_password_only_opens_the_right_username(app, username, accepted):
    """The real password with another username never logs in. (MySQL's
    `utf8mb4_unicode_ci` does make a case or accent variant of the
    username match -- `FUZZ-ADMIN`, `fúzz-admin` -- which is not an
    authentication bypass, since the password must still be right.)"""
    response = app.test_client().post("/api/auth/login", json={"username": username, "password": ADMIN_PASSWORD})
    check_response(response, "/api/auth/login")
    if accepted:
        assert response.status_code == 200 and _session_cookies(response)
    elif accepted is None:
        assert response.status_code in (200, 401)
    else:
        assert response.status_code == 401 and not _session_cookies(response), username


# ---------------------------------------------------------------------
# /login page
# ---------------------------------------------------------------------

@settings(max_examples=scaled(40))
@given(username=hostile_text, password=hostile_text, next_url=st.one_of(st.none(), hostile_text))
@example(username="<script>alert(1)</script>", password="x", next_url="<script>alert(1)</script>")
@example(username=ADMIN_USERNAME, password=ADMIN_PASSWORD.upper(), next_url="//evil.com")
@example(username=" ", password=" ", next_url=None)
def test_login_page_refuses_hostile_credentials(app, username, password, next_url):
    assume(not _is_real_password(password))
    client = app.test_client()
    data = {"username": username, "password": password, csrf.FIELD_NAME: with_csrf(client)}
    if next_url is not None:
        data["next"] = next_url
    response = client.post("/login", data=data)
    body = check_response(response, "/login", sent=[username, password, next_url or ""])
    assert response.status_code == 200, response.status_code
    assert not _session_cookies(response)
    assert "Invalid username or password." in body or "Enter a username and password." in body
    assert client.get(_ME).status_code == 401


def test_login_page_right_credentials_log_in(app):
    """Sanity for the refusals above."""
    client = app.test_client()
    response = client.post("/login", data={"username": ADMIN_USERNAME, "password": ADMIN_PASSWORD,
                                           csrf.FIELD_NAME: with_csrf(client)})
    assert response.status_code == 303 and _session_cookies(response)
    assert client.get(_ME).status_code == 200


# ---------------------------------------------------------------------
# API keys
# ---------------------------------------------------------------------

def test_real_api_key_works(app, real_api_key):
    response = app.test_client().get(_ME, headers={"Authorization": f"Bearer {real_api_key}"})
    assert response.status_code == 200
    assert response.get_json()["username"] == ADMIN_USERNAME


@settings(max_examples=scaled(60))
@given(value=hostile_text, prefix=st.sampled_from(["Bearer ", "bearer ", "BEARER ", "Bearer", "Bearer  ", "Basic ",
                                                  "Token ", "", " Bearer ", "Bearer\t"]))
@example(value="", prefix="Bearer ")
@example(value="pg_", prefix="Bearer ")
@example(value="pg_" + "A" * 43, prefix="Bearer ")
@example(value="' OR '1'='1", prefix="Bearer ")
@example(value="%", prefix="Bearer ")
def test_garbage_api_keys_are_refused(app, value, prefix):
    header = _header_safe(prefix + value)
    client = app.test_client()
    for method, path in _PROTECTED_API[:3]:
        response = client.open(path, method=method, headers={"Authorization": header})
        check_response(response, path, sent=[value])
        assert response.status_code == 401, (header, path, response.status_code)


def _mutations(secret):
    """Near misses of a real secret: every one must be refused."""
    swapped = secret.swapcase()
    out = [
        secret[:-1], secret[1:], secret + "A", "A" + secret, secret[: len(secret) // 2], swapped,
        secret[:-1] + ("A" if secret[-1] != "A" else "B"), secret[::-1], secret + " x", secret.upper(),
        secret.lower(), secret + "\t", secret * 2, secret.replace("_", "-"),
    ]
    for i in (0, len(secret) // 2, len(secret) - 1):
        ch = secret[i]
        flipped = "b" if ch == "a" else "a"
        out.append(secret[:i] + flipped + secret[i + 1:])
    return [value for value in dict.fromkeys(out) if value != secret and value.strip() != secret]


def test_near_miss_api_keys_are_refused(app, real_api_key):
    for value in _mutations(real_api_key):
        _assert_refused_everywhere(app, headers={"Authorization": f"Bearer {value}"})


def test_api_key_header_wins_over_a_valid_cookie(app, real_session):
    """A garbage `Authorization` header is checked instead of (not
    besides) a valid session cookie -- so it's refused, not silently
    upgraded by the cookie."""
    response = raw_client(app).get(_ME, headers={"Authorization": "Bearer nope",
                                                   "Cookie": f"{SESSION_COOKIE_NAME}={real_session}"})
    assert response.status_code == 401


def test_revoked_api_key_is_refused(app, admin_client):
    created = admin_client.post("/api/auth/api-keys", json={"label": "to revoke"}).get_json()
    headers = {"Authorization": f"Bearer {created['key']}"}
    assert app.test_client().get(_ME, headers=headers).status_code == 200
    assert admin_client.delete(f"/api/auth/api-keys/{created['id']}").status_code == 200
    _assert_refused_everywhere(app, headers=headers)
    # Revoking again, or someone else's/no key, is a clean 404.
    for key_id in (created["id"], 0, 2 ** 63 - 1, 2 ** 64):
        response = admin_client.delete(f"/api/auth/api-keys/{key_id}")
        check_response(response, "/api/auth/api-keys")
        assert response.status_code == 404


def test_session_cookie_is_not_an_api_key(app, real_session, real_api_key):
    """Each credential only works in its own slot."""
    assert app.test_client().get(_ME, headers={"Authorization": f"Bearer {real_session}"}).status_code == 401
    assert raw_client(app).get(_ME, headers={"Cookie": f"{SESSION_COOKIE_NAME}={real_api_key}"}).status_code == 401


# ---------------------------------------------------------------------
# Session cookie tampering
# ---------------------------------------------------------------------

def test_real_session_cookie_works(app, real_session):
    client = raw_client(app)
    assert client.get(_ME, headers={"Cookie": f"{SESSION_COOKIE_NAME}={real_session}"}).status_code == 200
    assert client.get("/admin", headers={"Cookie": f"{SESSION_COOKIE_NAME}={real_session}"}).status_code == 200


def test_near_miss_session_cookies_are_refused(app, real_session):
    for value in _mutations(real_session):
        _assert_refused_everywhere(app, cookie=f"{SESSION_COOKIE_NAME}={value}")


cookie_text = st.text(alphabet=st.characters(min_codepoint=0x21, max_codepoint=0x7e, blacklist_characters=';,'),
                      max_size=80)


@settings(max_examples=scaled(60))
@given(value=st.one_of(cookie_text, hostile_text.map(_header_safe)))
@example(value="")
@example(value='""')
@example(value="' OR '1'='1")
@example(value="%00")
@example(value="a" * 5000)
def test_garbage_session_cookies_are_refused(app, value):
    value = value.replace(";", "")
    client = raw_client(app)
    for path in (_ME, "/admin"):
        response = client.get(path, headers={"Cookie": f"{SESSION_COOKIE_NAME}={value}"})
        check_response(response, path, sent=[value])
        assert response.status_code in (401, 302), (value, path, response.status_code)


def test_session_cookie_with_junk_around_it(app, real_session):
    """Other cookies, duplicates and odd spacing around the real session
    cookie: never a 5xx; with the real token absent never accepted."""
    client = raw_client(app)
    for cookie, may_accept in [
        (f"{SESSION_COOKIE_NAME}=bad; {SESSION_COOKIE_NAME}={real_session}", True),
        (f"{SESSION_COOKIE_NAME}={real_session}; {SESSION_COOKIE_NAME}=bad", True),
        (f"x=1; {SESSION_COOKIE_NAME}={real_session}; y=2", True),
        (f"{SESSION_COOKIE_NAME.upper()}={real_session}", False),
        (f"{SESSION_COOKIE_NAME}_x={real_session}", False),
        (f"x{SESSION_COOKIE_NAME}={real_session}", False),
        (f'{SESSION_COOKIE_NAME}="{real_session}"', True),
        (f"{SESSION_COOKIE_NAME}", False), (f"{SESSION_COOKIE_NAME}=", False), (";;;", False), ("=", False),
    ]:
        response = client.get(_ME, headers={"Cookie": cookie})
        check_response(response, _ME)
        assert response.status_code in ((200, 401) if may_accept else (401,)), (cookie, response.status_code)


def test_logged_out_and_expired_sessions_are_refused(app, fuzz_db):
    client = app.test_client()
    token = _session_token(login_api(client))
    assert client.get(_ME).status_code == 200
    assert client.post("/api/auth/logout").status_code == 200
    _assert_refused_everywhere(app, cookie=f"{SESSION_COOKIE_NAME}={token}")

    token = _session_token(login_api(app.test_client()))
    conn = _db.get_control_connection(fuzz_db["config"], ensure_schema=False)
    try:
        conn.execute("UPDATE admin_sessions SET expires_at = CURRENT_TIMESTAMP - INTERVAL 1 SECOND "
                     "WHERE token_hash = ?", (adminAuth._hash_token(token),))
        conn.commit()
    finally:
        conn.close()
    _assert_refused_everywhere(app, cookie=f"{SESSION_COOKIE_NAME}={token}")


# ---------------------------------------------------------------------
# The flash cookie and the CSRF cookie
# ---------------------------------------------------------------------

@settings(max_examples=scaled(30))
@given(value=st.one_of(cookie_text, hostile_text.map(_header_safe)))
def test_garbage_flash_cookie_is_ignored(admin_client, value):
    value = value.replace(";", "")
    admin_client.set_cookie(FLASH_COOKIE, value, domain="localhost")
    try:
        response = admin_client.get("/admin")
        body = check_response(response, "/admin", sent=[value])
        assert response.status_code == 200
        assert "new-key" not in body or 'id="new-key"' not in body
    finally:
        admin_client.delete_cookie(FLASH_COOKIE, domain="localhost")


@pytest.mark.parametrize("secret,salt,accepted", [
    (SECRET_KEY, "planetgen-web-flash", True),
    ("some-other-secret", "planetgen-web-flash", False),
    (SECRET_KEY, "another-salt", False),
])
def test_forged_flash_cookie_is_ignored(admin_client, secret, salt, accepted):
    """A flash cookie is only shown when signed with the app's secret
    (so nobody can make the admin page display a key or message of
    their choosing)."""
    marker = "pg_FORGED_KEY_MARKER"
    forged = URLSafeTimedSerializer(secret, salt=salt).dumps({"new_key": {"label": "forged", "key": marker},
                                                              "message": "<script>alert(1)</script>"})
    admin_client.set_cookie(FLASH_COOKIE, forged, domain="localhost")
    try:
        response = admin_client.get("/admin")
        body = check_response(response, "/admin", sent=["<script>alert(1)</script>"])
        assert (marker in body) is accepted
    finally:
        admin_client.delete_cookie(FLASH_COOKIE, domain="localhost")


def test_flash_cookie_of_the_wrong_shape_is_ignored(admin_client):
    for data in ([1, 2], "text", 5, None, {"new_key": "not-a-dict"}, {"new_key": {"label": 1}}):
        admin_client.set_cookie(FLASH_COOKIE, URLSafeTimedSerializer(SECRET_KEY, salt="planetgen-web-flash").dumps(data),
                                domain="localhost")
        try:
            response = admin_client.get("/admin")
            check_response(response, "/admin")
            assert response.status_code == 200, data
        finally:
            admin_client.delete_cookie(FLASH_COOKIE, domain="localhost")


@settings(max_examples=scaled(30))
@given(nonce=cookie_text)
def test_csrf_cookie_tampering_is_refused(app, nonce):
    """The form token is bound to the cookie's nonce: swapping the
    cookie for anything else invalidates it."""
    client = app.test_client()
    token = with_csrf(client)
    assume(nonce != "n" * 40)
    client.set_cookie(csrf.COOKIE_NAME, nonce, domain="localhost")
    response = client.post("/login", data={"username": ADMIN_USERNAME, "password": ADMIN_PASSWORD,
                                           csrf.FIELD_NAME: token})
    check_response(response, "/login")
    assert response.status_code == 400
    assert not _session_cookies(response)


# ---------------------------------------------------------------------
# Changing credentials
# ---------------------------------------------------------------------

@settings(max_examples=scaled(30))
@given(current=st.one_of(hostile_text, st.integers(), st.lists(hostile_text, max_size=2)),
       new_username=st.one_of(hostile_text, st.none(), st.integers()),
       new_password=st.one_of(hostile_text, st.none()))
@example(current="", new_username=ADMIN_USERNAME, new_password="another-strong-password")
@example(current=ADMIN_PASSWORD.upper(), new_username="x", new_password="another-strong-password")
def test_change_credentials_with_a_wrong_current_password(app, current, new_username, new_password):
    """Anything but the right current password is a 400 (never a 5xx)
    and changes nothing."""
    assume(not _is_real_password(current))
    client = app.test_client()
    login_api(client)
    response = client.post("/api/auth/change-credentials", json={
        "current_password": current, "new_username": new_username, "new_password": new_password})
    check_response(response, "/api/auth/change-credentials")
    assert response.status_code == 400
    assert not _session_cookies(response)
    assert client.get(_ME).get_json()["username"] == ADMIN_USERNAME
    login_api(app.test_client())  # the old credentials still work


@pytest.mark.parametrize("new_username,new_password", [
    ("", "another-strong-password"), ("   ", "another-strong-password"), ("x", ""), ("x", "short"),
    ("x", "password"), ("samesame-samesame", "SAMESAME-samesame"), ("x", None), (None, "another-strong-password"),
])
def test_change_credentials_policy_violations(app, new_username, new_password):
    client = app.test_client()
    login_api(client)
    response = client.post("/api/auth/change-credentials", json={
        "current_password": ADMIN_PASSWORD, "new_username": new_username, "new_password": new_password})
    check_response(response, "/api/auth/change-credentials")
    assert response.status_code == 400
    login_api(app.test_client())


def test_change_credentials_requires_a_session(app):
    response = app.test_client().post("/api/auth/change-credentials", json={
        "current_password": ADMIN_PASSWORD, "new_username": "hijack", "new_password": "hijack-password-123"})
    check_response(response, "/api/auth/change-credentials")
    assert response.status_code == 401
    login_api(app.test_client())


# A2 (fixed): change-credentials never checks new_username's length; past admin_users.username's VARCHAR(64) MySQL raises DataError -> 500
def test_regression_change_credentials_overlong_username_is_a_400(app):
    client = app.test_client()
    login_api(client)
    response = client.post("/api/auth/change-credentials", json={
        "current_password": ADMIN_PASSWORD, "new_username": "u" * 65, "new_password": "another-strong-password"})
    try:
        assert response.status_code == 400
    finally:
        login_api(app.test_client())  # unchanged either way


# A1 (fixed): change-credentials calls .strip() on a non-string new_username -> AttributeError -> 500
@pytest.mark.parametrize("new_username", [5, ["x"], {"a": 1}])
def test_regression_change_credentials_non_string_username_is_a_400(app, new_username):
    client = app.test_client()
    login_api(client)
    response = client.post("/api/auth/change-credentials", json={
        "current_password": ADMIN_PASSWORD, "new_username": new_username, "new_password": "another-strong-password"})
    assert response.status_code == 400


# ---------------------------------------------------------------------
# API key labels
# ---------------------------------------------------------------------

@settings(max_examples=scaled(30))
@given(label=st.one_of(hostile_text, st.integers(), st.none(), st.lists(hostile_text, max_size=2)))
@example(label="")
@example(label="   ")
@example(label="<script>alert(1)</script>")
@example(label="l" * 128)
@example(label="l" * 129)  # B9
@example(label=5)  # A1
def test_api_key_labels(admin_client, label):
    response = admin_client.post("/api/auth/api-keys", json={"label": label})
    check_response(response, "/api/auth/api-keys", sent=[str(label)])
    if isinstance(label, str) and label.strip() and len(label.strip()) <= 128:
        assert response.status_code == 201
        assert response.get_json()["key"].startswith("pg_")
        listing = admin_client.get("/admin")
        check_response(listing, "/admin", sent=[label])
    else:
        assert response.status_code == 400


# A1 (fixed): POST /api/auth/api-keys calls .strip() on a non-string label -> AttributeError -> 500
@pytest.mark.parametrize("label", [5, ["x"], {"a": 1}, True])
def test_regression_api_key_non_string_label_is_a_400(admin_client, label):
    response = admin_client.post("/api/auth/api-keys", json={"label": label})
    assert response.status_code == 400


def test_shared_admin_client_still_logged_in(admin_client):
    """Guards the module's other tests: nothing above logged the shared
    client out or rotated its credentials."""
    assert admin_client.get(_ME).get_json() == {"username": ADMIN_USERNAME, "must_change_credentials": False}
