# tests/test_web_admin.py

"""
The admin pages on the Flask app (`src/planetgen/web/admin_pages.py`):
`/login`, `/logout`, `/account`, `/admin` and `/admin/stats`, which
replaced `login.py`, `logout.py`, `changecreds.py`, `admin.py` and
`adminstats.py` (those old URLs now redirect here).

Most tests fake the data layer by monkeypatching the `apiclient`
functions the views call (the pattern of `test_web_pages.py`). The tests
at the bottom run the real login flow against a throwaway database
through the in-process transport, and are skipped without a MySQL test
server.
"""

import re
from urllib.parse import parse_qs, urlsplit

import pytest

from planetgen.web.app import create_app
from planetgen.api.authz import SESSION_COOKIE_NAME
from planetgen.api.config import Config

from planetgen.web.lib import apiclient  # noqa: E402
from planetgen.db import store as _db  # noqa: E402
from planetgen.admin import auth as adminAuth
from planetgen.web import admin_pages, csrf  # noqa: E402

DB = "planetgen_web_test"
NONCE = "n" * 43
SESSION_SET = f"{SESSION_COOKIE_NAME}=new-token; Path=/; HttpOnly; SameSite=Strict; Max-Age=43200"
SESSION_CLEAR = (f"{SESSION_COOKIE_NAME}=; Expires=Thu, 01 Jan 1970 00:00:00 GMT; Max-Age=0; "
                 f"Path=/")


class _FakeConfig(Config):
    WEB_DATABASE = DB
    SESSION_COOKIE_SECURE = False
    SECRET_KEY = "test-secret"


def _key(i, revoked=False):
    return {"id": i, "label": f"key {i:03d}", "created_at": "2026-01-01 00:00:00",
            "last_used_at": None, "revoked_at": "2026-02-01 00:00:00" if revoked else None}


def _stats(reachable=True):
    database = {
        "name": DB, "reachable": reachable, "schema_current": True, "schema_version": 30,
        "schema_expected": 30, "size_bytes": 5 * 1024 * 1024,
        "counts": {"sectors": 12, "star_systems": 3456},
        "name_collisions": {"distinct_base_names": 2, "sector": 1, "system": 1},
        "timestamps": [{"table": "sectors", "newest_created_at": "2026-09-01 10:00:00",
                        "last_modified_at": "2026-09-02 10:00:00"},
                       {"table": "star_systems", "newest_created_at": None, "last_modified_at": None}],
        "tables": [{"name": "bright_stars", "approx_rows": 61234, "data_bytes": 4096, "index_bytes": 1024},
                   {"name": "star_systems", "approx_rows": 3456, "data_bytes": 2048, "index_bytes": 1024}],
    }
    if not reachable:
        database = {"name": DB, "reachable": False, "detail": "connection <refused>"}
    return {
        "api": {"version": "5.53.1", "uptime_seconds": 3700, "python_version": "3.12.1",
                "memory": {"available_bytes": 2 ** 30, "total_bytes": 2 ** 32}, "load_average": [0.1, 0.2, 0.3]},
        "mysql": {"version": "10.11", "uptime_seconds": 90000, "threads_connected": 4},
        "query_ms": 3.2,
        "database": database,
    }


class FakeAuth:
    """Stands in for the `apiclient` auth/admin functions; records calls."""

    def __init__(self):
        self.admin = None
        self.calls = []
        self.database_extra = {}
        self.keys = [_key(1), _key(2, revoked=True)]
        self.login_error = None
        self.failures = []
        self.lockouts = []
        self.change_error = None
        self.names = [{"base_name": "Vega", "levels": ["system"], "rows": [
            {"kind": "sector", "id": 5, "name": "Vega Alpha"},
            {"kind": "system", "id": 77, "name": "Vega <b>Beta</b>"},
        ]}]

    def auth_me(self, cookie_header):
        self.calls.append(("auth_me", cookie_header))
        return self.admin

    def auth_login(self, username, password, cookie_header=None):
        self.calls.append(("auth_login", username, password))
        if self.login_error:
            raise self.login_error
        return {"username": username, "must_change_credentials": password == "password"}, [SESSION_SET]

    def auth_logout(self, cookie_header):
        self.calls.append(("auth_logout", cookie_header))
        return [SESSION_CLEAR]

    def auth_change_credentials(self, cookie_header, current_password, new_username, new_password):
        self.calls.append(("auth_change_credentials", cookie_header, current_password, new_username, new_password))
        if self.change_error:
            raise self.change_error
        return {"username": new_username, "must_change_credentials": False}, [SESSION_SET]

    def auth_list_api_keys(self, cookie_header):
        self.calls.append(("auth_list_api_keys", cookie_header))
        return self.keys

    def auth_create_api_key(self, cookie_header, label):
        self.calls.append(("auth_create_api_key", cookie_header, label))
        return {"id": 99, "label": label, "key": "pgk_secret_value_123"}

    def auth_revoke_api_key(self, cookie_header, key_id):
        self.calls.append(("auth_revoke_api_key", cookie_header, key_id))

    def admin_set_sector_wiki_url(self, cookie_header, db, sector_id, wiki_url):
        self.calls.append(("admin_set_sector_wiki_url", db, sector_id, wiki_url))

    def admin_stats(self, cookie_header, db):
        self.calls.append(("admin_stats", db))
        stats = _stats()
        stats["database"].update(self.database_extra)
        return stats

    def admin_duplicate_names(self, cookie_header, db, limit=None, offset=None):
        self.calls.append(("admin_duplicate_names", db, limit, offset))
        return {"items": self.names[offset:offset + limit], "total": len(self.names),
                "limit": limit, "offset": offset}

    generation = {"available": True, "sizes": {}, "buckets": []}

    def admin_generation_stats(self, cookie_header):
        self.calls.append(("admin_generation_stats",))
        return self.generation

    naming = {"database": "planetgen", "key": "0A1B2C3D", "codec_version": 1, "current_codec_version": 1,
              "drawn_at": "2026-10-01T10:00:00Z", "changed_at": None, "changed_by": None}

    def admin_naming_key(self, cookie_header, db):
        self.calls.append(("admin_naming_key", db))
        return self.naming

    settings_files = [{"name": "AB-CD-20261009-010000Z.json", "seed": "AB", "version_key": "CD",
                       "created_at": "2026-10-09T01:00:00Z", "size": 2048, "current": True}]

    def admin_galaxy_settings(self, cookie_header):
        return {"items": self.settings_files}

    def admin_galaxy_settings_file(self, cookie_header, name):
        if name != self.settings_files[0]["name"]:
            raise apiclient.ApiError("No such settings file.", status_code=404)
        return {"format": 1, "seed": "AB"}

    def admin_set_naming_key(self, cookie_header, db, key=None, draw=False):
        self.calls.append(("admin_set_naming_key", db, key, draw))

    def admin_lockouts(self, cookie_header):
        self.calls.append(("admin_lockouts",))
        return {"items": self.lockouts, "proxy_warning": False}

    def admin_lift_lockout(self, cookie_header, scope=None, subject=None, lift_all=False):
        self.calls.append(("admin_lift_lockout", scope, subject, lift_all))
        return 1

    def admin_login_failures(self, cookie_header, limit=None):
        self.calls.append(("admin_login_failures",))
        return self.failures

    def auth_totp_status(self, cookie_header):
        self.calls.append(("auth_totp_status",))
        return {"enabled": False, "recovery_codes_left": 0}

    def called(self, name):
        return [call for call in self.calls if call[0] == name]


_FAKED = ("auth_me", "auth_login", "auth_logout", "auth_change_credentials", "auth_list_api_keys",
          "auth_create_api_key", "auth_revoke_api_key", "admin_set_sector_wiki_url", "admin_stats",
          "admin_duplicate_names", "admin_login_failures", "admin_lockouts", "admin_generation_stats",
          "admin_lift_lockout", "auth_totp_status", "admin_naming_key", "admin_set_naming_key",
          "admin_galaxy_settings", "admin_galaxy_settings_file")


@pytest.fixture
def fake(monkeypatch):
    data = FakeAuth()
    for name in _FAKED:
        monkeypatch.setattr(apiclient, name, getattr(data, name))
    return data


@pytest.fixture
def app():
    application = create_app(_FakeConfig)
    application.testing = True
    return application


@pytest.fixture
def client(app):
    test_client = app.test_client()
    test_client.set_cookie(csrf.COOKIE_NAME, NONCE)
    return test_client


@pytest.fixture
def token(app):
    """A CSRF token for a visitor who isn't logged in (the login form)."""
    with app.app_context():
        return csrf._sign(NONCE)


LOGGED_IN_SESSION = "token-value"


@pytest.fixture
def admin_token(app):
    """A CSRF token for the session `_logged_in` sets: tokens are bound
    to the login session (security #49)."""
    with app.app_context():
        return csrf._sign(NONCE, LOGGED_IN_SESSION)


def _logged_in(client, fake, must_change=False):
    fake.admin = {"username": "boss", "must_change_credentials": must_change}
    client.set_cookie(SESSION_COOKIE_NAME, LOGGED_IN_SESSION)


_FORM_TOKEN = re.compile(rf'name="{csrf.FIELD_NAME}" value="([^"]+)"')


def _form_token(client, path):
    """The CSRF token the page at `path` renders into its forms, as a
    browser would submit it."""
    return _FORM_TOKEN.search(client.get(path).get_data(as_text=True)).group(1)


def _set_cookies(resp):
    return resp.headers.getlist("Set-Cookie")


def _redirect(resp):
    """`(path, next)` of a redirect's Location."""
    parts = urlsplit(resp.headers["Location"])
    return parts.path, parse_qs(parts.query).get("next", [None])[0]


# --- Routing -------------------------------------------------------------------------

def test_admin_pages_are_routes(app):
    for name in ("login", "logout", "account", "admin", "admin_stats"):
        assert f"web.{name}" in app.view_functions


def test_header_links_point_at_new_urls(client, fake):
    html = client.get("/login").get_data(as_text=True)
    assert '<a href="/login" aria-current="page">Login</a>' in html
    _logged_in(client, fake)
    html = client.get("/admin").get_data(as_text=True)
    assert '<a href="/admin" aria-current="page">Admin</a>' in html
    assert '<a href="/admin/stats">Stats</a>' in html  # the admin tab row
    assert '<a href="/logout">Logout</a>' in html


# --- Login gate and ?next= -------------------------------------------------------------

@pytest.mark.parametrize("path", ["/admin", "/admin/stats?names_page=2", "/account"])
def test_unauthenticated_admin_pages_redirect_to_login(client, fake, path):
    resp = client.get(path)
    assert resp.status_code == 302
    assert _redirect(resp) == ("/login", path)
    assert "boss" not in resp.get_data(as_text=True)
    assert not fake.called("auth_list_api_keys") and not fake.called("admin_stats")


def test_unauthenticated_post_to_admin_changes_nothing(client, fake, token):
    resp = client.post("/admin", data={csrf.FIELD_NAME: token, "action": "create_key", "label": "x"})
    assert resp.status_code == 302
    assert resp.headers["Location"].startswith("/login?next=")
    assert not fake.called("auth_create_api_key")


def test_default_credentials_go_to_account_first(client, fake):
    _logged_in(client, fake, must_change=True)
    resp = client.get("/admin/stats")
    assert resp.status_code == 302
    assert _redirect(resp) == ("/account", "/admin/stats")
    assert client.get("/account").status_code == 200


@pytest.mark.parametrize("value,expected", [
    ("/admin/stats?names_page=2", "/admin/stats?names_page=2"),
    ("/", "/"),
    (None, "/admin"),
    ("", "/admin"),
    ("https://evil.example/", "/admin"),
    ("//evil.example/", "/admin"),
    ("/\\evil.example/", "/admin"),
    ("\\\\evil.example", "/admin"),
    ("javascript:alert(1)", "/admin"),
    ("admin", "/admin"),
    ("/ /evil.example", "/admin"),
    ("/\t/evil.example", "/admin"),
    ("/%0d%0aSet-Cookie:x", "/%0d%0aSet-Cookie:x"),
    ("/" + "a" * 3000, "/admin"),
])
def test_safe_next(app, value, expected):
    with app.test_request_context("/"):
        assert admin_pages.safe_next(value) == expected


# --- /login --------------------------------------------------------------------------------

def test_login_form(client, fake):
    resp = client.get("/login?next=/admin/stats")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert resp.headers["Cache-Control"] == "no-store"
    assert '<form method="post" action="/login"' in html
    assert f'name="{csrf.FIELD_NAME}"' in html
    assert '<input type="hidden" name="next" value="/admin/stats">' in html
    assert 'autocomplete="current-password"' in html


def test_login_form_drops_unsafe_next(client, fake):
    html = client.get("/login?next=//evil.example").get_data(as_text=True)
    assert '<input type="hidden" name="next" value="/admin">' in html
    assert "evil.example" not in html


def test_login_requires_csrf_token(client, fake):
    resp = client.post("/login", data={"username": "a", "password": "b"})
    assert resp.status_code == 400
    assert not fake.called("auth_login")


def test_login_success_relays_cookie_and_redirects_to_next(client, fake, token):
    resp = client.post("/login", data={csrf.FIELD_NAME: token, "username": "boss", "password": "long-password-1",
                                       "next": "/admin/stats?names_page=2"})
    assert resp.status_code == 303
    assert resp.headers["Location"] == "/admin/stats?names_page=2"
    assert SESSION_SET in _set_cookies(resp)
    assert fake.called("auth_login") == [("auth_login", "boss", "long-password-1")]


def test_login_never_redirects_off_site(client, fake, token):
    resp = client.post("/login", data={csrf.FIELD_NAME: token, "username": "boss", "password": "long-password-1",
                                       "next": "//evil.example/"})
    assert resp.headers["Location"] == "/admin"


def test_login_with_default_credentials_goes_to_account(client, fake, token):
    resp = client.post("/login", data={csrf.FIELD_NAME: token, "username": "admin", "password": "password",
                                       "next": "/admin/stats"})
    assert resp.status_code == 303
    assert _redirect(resp) == ("/account", "/admin/stats")
    assert SESSION_SET in _set_cookies(resp)


def test_login_wrong_password_shows_form_again(client, fake, token):
    fake.login_error = apiclient.ApiError("planetGen API error (401): invalid credentials", status_code=401)
    resp = client.post("/login", data={csrf.FIELD_NAME: token, "username": "<b>x</b>", "password": "nope"})
    html = resp.get_data(as_text=True)
    assert resp.status_code == 401
    assert "Invalid username or password." in html
    assert 'value="&lt;b&gt;x&lt;/b&gt;"' in html
    assert "nope" not in html
    assert not any(SESSION_COOKIE_NAME in header for header in _set_cookies(resp))


def test_login_rate_limited_message(client, fake, token):
    fake.login_error = apiclient.ApiError("planetGen API error (429): too many", status_code=429)
    resp = client.post("/login", data={csrf.FIELD_NAME: token, "username": "a", "password": "b"})
    assert resp.status_code == 429
    assert "Too many login attempts" in resp.get_data(as_text=True)


def test_login_missing_fields(client, fake, token):
    resp = client.post("/login", data={csrf.FIELD_NAME: token, "username": " ", "password": ""})
    assert resp.status_code == 200
    assert not fake.called("auth_login")


def test_login_api_down_is_502(client, fake, token):
    fake.login_error = apiclient.ApiError("planetGen API error (500): SELECT secret", status_code=500)
    resp = client.post("/login", data={csrf.FIELD_NAME: token, "username": "a", "password": "b"})
    assert resp.status_code == 502
    assert "SELECT secret" not in resp.get_data(as_text=True)


def test_login_page_when_already_logged_in_goes_to_next(client, fake):
    _logged_in(client, fake)
    resp = client.get("/login?next=/admin/stats")
    assert resp.status_code == 302
    assert resp.headers["Location"] == "/admin/stats"


@pytest.mark.parametrize("path,form", [
    ("/login", {"username": "boss", "password": "long-password-1"}),
    ("/login/code", {"username": "boss", "pending": "signed-pending", "code": "123456"}),
])
def test_a_sign_in_form_sent_again_after_signing_in_never_shows_form_expired(client, fake, token, path, form):
    # SEC.31: the login (or two-step code) form, sent a second time after
    # the first sign-in set the session cookie, carries the token minted
    # before sign-in. It goes back to the login page, which forwards the
    # signed-in admin to `next`; no 400 and no second sign-in.
    _logged_in(client, fake)
    resp = client.post(path, data={csrf.FIELD_NAME: token, "next": "/admin/stats", **form})
    assert resp.status_code == 303
    assert _redirect(resp) == ("/login", "/admin/stats")
    assert not fake.called("auth_login") and not fake.called("auth_login_totp")
    follow = client.get(resp.headers["Location"])
    assert follow.status_code == 302
    assert follow.headers["Location"] == "/admin/stats"


def test_a_stale_sign_in_form_without_a_session_is_still_refused(client, fake):
    resp = client.post("/login/code", data={csrf.FIELD_NAME: "wrong", "pending": "p", "code": "1"})
    assert resp.status_code == 400
    assert not fake.called("auth_login_totp")


# --- /logout -------------------------------------------------------------------------------

def test_logout_get_changes_nothing(client, fake):
    _logged_in(client, fake)
    resp = client.get("/logout")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert not fake.called("auth_logout")
    assert not any(SESSION_COOKIE_NAME in header for header in _set_cookies(resp))
    assert '<form method="post" action="/logout"' in html
    assert f'name="{csrf.FIELD_NAME}"' in html
    assert "boss" in html


def test_logout_get_when_logged_out(client, fake):
    html = client.get("/logout").get_data(as_text=True)
    assert "You are not signed in." in html
    assert '<form method="post" action="/logout"' not in html


def test_logout_post_requires_csrf(client, fake):
    _logged_in(client, fake)
    assert client.post("/logout").status_code == 400
    assert not fake.called("auth_logout")


def test_logout_post_clears_cookie(client, fake, admin_token):
    _logged_in(client, fake)
    resp = client.post("/logout", data={csrf.FIELD_NAME: admin_token})
    assert resp.status_code == 303
    assert resp.headers["Location"] == "/"
    assert SESSION_CLEAR in _set_cookies(resp)
    assert f"{SESSION_COOKIE_NAME}=token-value" in fake.called("auth_logout")[0][1]


# --- /account ------------------------------------------------------------------------------

def test_account_form(client, fake):
    _logged_in(client, fake, must_change=True)
    html = client.get("/account?next=/admin/stats").get_data(as_text=True)
    assert "still uses the first password the installer printed" in html
    assert 'name="new_username" value="boss"' in html
    assert '<input type="hidden" name="next" value="/admin/stats">' in html
    assert f'name="{csrf.FIELD_NAME}"' in html


def test_account_change_relays_new_cookie(client, fake, admin_token):
    _logged_in(client, fake, must_change=True)
    resp = client.post("/account", data={csrf.FIELD_NAME: admin_token, "current_password": "password",
                                         "new_username": "boss2", "new_password": "a-long-new-password",
                                         "next": "/admin/stats"})
    assert resp.status_code == 303
    assert resp.headers["Location"] == "/admin/stats"
    assert SESSION_SET in _set_cookies(resp)
    call = fake.called("auth_change_credentials")[0]
    assert call[2:] == ("password", "boss2", "a-long-new-password")


def test_account_change_to_admin_flashes_a_message(client, fake, admin_token):
    _logged_in(client, fake)
    resp = client.post("/account", data={csrf.FIELD_NAME: admin_token, "current_password": "x",
                                         "new_username": "boss", "new_password": "a-long-new-password"})
    assert resp.headers["Location"] == "/admin"
    html = client.get("/admin").get_data(as_text=True)
    assert "Your username and password were changed." in html


def test_account_policy_error_shows_inline(client, fake, admin_token):
    _logged_in(client, fake)
    fake.change_error = apiclient.ApiError("planetGen API error (400): password is too short", status_code=400)
    resp = client.post("/account", data={csrf.FIELD_NAME: admin_token, "current_password": "x",
                                         "new_username": "<i>new</i>", "new_password": "short"})
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert '<p class="error" role="alert">password is too short</p>' in html
    assert 'value="&lt;i&gt;new&lt;/i&gt;"' in html
    assert not any(SESSION_COOKIE_NAME in header for header in _set_cookies(resp))


def test_account_requires_login(client, fake, token):
    resp = client.post("/account", data={csrf.FIELD_NAME: token, "current_password": "x"})
    assert resp.status_code == 302
    assert not fake.called("auth_change_credentials")


# --- /admin --------------------------------------------------------------------------------

def test_admin_lists_keys(client, fake):
    _logged_in(client, fake)
    resp = client.get("/admin")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert resp.headers["Cache-Control"] == "no-store"
    assert "key 001" in html and "key 002" in html
    assert "revoked 2026-02-01 00:00 UTC" in html
    # One revoke form (the active key), carrying the CSRF token.
    assert html.count('value="revoke_key"') == 1
    revoke = re.search(r'<form method="post" action="/admin" class="table-form">.*?</form>', html, re.S).group(0)
    assert f'name="{csrf.FIELD_NAME}"' in revoke and 'name="key_id" value="1"' in revoke
    assert 'name="db"' not in html


def test_admin_keys_are_a_data_table_paged_with_get_links(client, fake):
    _logged_in(client, fake)
    fake.keys = [_key(i) for i in range(1, 61)]
    html = client.get("/admin?keys_page=2").get_data(as_text=True)
    assert "key 051" in html and "key 050" not in html
    assert 'href="/admin?keys_page=1#api-keys"' in html


def test_admin_keys_table_route_is_for_admins_only_and_carries_the_csrf_token(client, fake):
    assert client.get("/table/api-keys").status_code == 403
    _logged_in(client, fake)
    fake.keys = [_key(1), _key(2, revoked=True)]
    data = client.get("/table/api-keys?sort=label&facets=1").get_json()
    assert data["total"] == 2 and data["rows"][0][0]["text"] == "key 001"
    form = data["rows"][0][4]["form"]
    assert form["action"] == "/admin" and ["key_id", 1] in form["fields"]
    assert [csrf.FIELD_NAME] == [name for name, _ in form["fields"] if name == csrf.FIELD_NAME]
    assert "form" not in data["rows"][1][4]
    only = client.get("/table/api-keys?keys_status=revoked").get_json()
    assert [row[0]["text"] for row in only["rows"]] == ["key 002"]


def test_admin_create_key_is_post_redirect_get_and_shown_once(client, fake, admin_token):
    _logged_in(client, fake)
    resp = client.post("/admin", data={csrf.FIELD_NAME: admin_token, "action": "create_key", "label": "my script"})
    assert resp.status_code == 303
    assert resp.headers["Location"] == "/admin#new-key"
    assert "pgk_secret_value_123" not in resp.headers["Location"]
    assert fake.called("auth_create_api_key")[0][2] == "my script"
    flash = [h for h in _set_cookies(resp) if h.startswith(f"{admin_pages.FLASH_COOKIE}=")][0]
    assert "HttpOnly" in flash and "SameSite=Strict" in flash

    page = client.get("/admin")
    html = page.get_data(as_text=True)
    assert "pgk_secret_value_123" in html
    assert "New API key: my script" in html
    assert any(h.startswith(f"{admin_pages.FLASH_COOKIE}=;") for h in _set_cookies(page))
    assert "pgk_secret_value_123" not in client.get("/admin").get_data(as_text=True)


def test_new_key_flash_is_only_sent_back_to_admin(client, fake, admin_token):
    """Security #50: the flash cookie carrying a new key's raw value is
    scoped to `/admin` (the page that shows it), so the browser doesn't
    send it with requests to the rest of the site."""
    _logged_in(client, fake)
    resp = client.post("/admin", data={csrf.FIELD_NAME: admin_token, "action": "create_key", "label": "k"})
    flash = [h for h in _set_cookies(resp) if h.startswith(f"{admin_pages.FLASH_COOKIE}=")][0]
    assert "Path=/admin" in flash and "Path=/;" not in flash
    for path in ("/logout", "/login", "/account"):
        sent = client.get(path).request.headers.get("Cookie", "")
        assert admin_pages.FLASH_COOKIE not in sent, path
    page = client.get("/admin")
    assert admin_pages.FLASH_COOKIE in page.request.headers.get("Cookie", "")
    assert "pgk_secret_value_123" in page.get_data(as_text=True)
    cleared = [h for h in _set_cookies(page) if h.startswith(f"{admin_pages.FLASH_COOKIE}=;")][0]
    assert "Path=/admin" in cleared  # deleting must match the path it was set with
    assert client.get_cookie(admin_pages.FLASH_COOKIE, path="/admin") is None


def test_account_change_flash_is_scoped_to_admin(client, fake, admin_token):
    _logged_in(client, fake)
    resp = client.post("/account", data={csrf.FIELD_NAME: admin_token, "current_password": "x",
                                         "new_username": "boss", "new_password": "a-long-new-password",
                                         "next": "/admin"})
    flash = [h for h in _set_cookies(resp) if h.startswith(f"{admin_pages.FLASH_COOKIE}=")][0]
    assert "Path=/admin" in flash


# --- CSRF tokens are bound to the login session (security #49) ---------------------------------

def test_csrf_token_for_one_session_fails_for_another(app, client, fake):
    _logged_in(client, fake)
    with app.app_context():
        other_session = csrf._sign(NONCE, "someone-elses-session")
        anonymous = csrf._sign(NONCE)
        mine = csrf._sign(NONCE, LOGGED_IN_SESSION)
    for bad in (other_session, anonymous):
        resp = client.post("/admin", data={csrf.FIELD_NAME: bad, "action": "revoke_key", "key_id": "1"})
        assert resp.status_code == 400
    assert not fake.called("auth_revoke_api_key")
    resp = client.post("/admin", data={csrf.FIELD_NAME: mine, "action": "revoke_key", "key_id": "1"})
    assert resp.status_code == 303
    assert fake.called("auth_revoke_api_key")


def test_logged_in_token_fails_without_the_session(client, fake, admin_token):
    """A token rendered for a session can't be replayed with the session
    cookie dropped (or before logging in)."""
    resp = client.post("/login", data={csrf.FIELD_NAME: admin_token, "username": "boss",
                                       "password": "long-password-1"})
    assert resp.status_code == 400
    assert not fake.called("auth_login")


def test_rendered_token_follows_the_session(client, fake):
    """Tokens are computed per render: the login form's token is for no
    session, and after logging in the pages render one for the new
    session."""
    anonymous = _form_token(client, "/login")
    _logged_in(client, fake)
    logged_in = _form_token(client, "/admin")
    assert anonymous != logged_in
    resp = client.post("/admin", data={csrf.FIELD_NAME: logged_in, "action": "revoke_key", "key_id": "1"})
    assert resp.status_code == 303


def test_admin_forged_flash_is_ignored(client, fake):
    _logged_in(client, fake)
    client.set_cookie(admin_pages.FLASH_COOKIE, "eyJuZXdfa2V5Ijp7ImtleSI6ImZvcmdlZCJ9fQ.bad.sig")
    html = client.get("/admin").get_data(as_text=True)
    assert "forged" not in html and "New API key" not in html


def test_admin_create_key_needs_label(client, fake, admin_token):
    _logged_in(client, fake)
    resp = client.post("/admin", data={csrf.FIELD_NAME: admin_token, "action": "create_key", "label": "  "})
    assert resp.headers["Location"] == "/admin#api-keys"
    assert not fake.called("auth_create_api_key")
    assert "Label is required." in client.get("/admin").get_data(as_text=True)


def test_admin_revoke_key_returns_to_the_keys(client, fake, admin_token):
    _logged_in(client, fake)
    resp = client.post("/admin", data={csrf.FIELD_NAME: admin_token, "action": "revoke_key", "key_id": "1"})
    assert resp.status_code == 303
    assert resp.headers["Location"] == "/admin#api-keys"
    assert fake.called("auth_revoke_api_key")[0][2] == 1


def test_admin_revoke_bad_id(client, fake, admin_token):
    _logged_in(client, fake)
    client.post("/admin", data={csrf.FIELD_NAME: admin_token, "action": "revoke_key", "key_id": "x"})
    assert not fake.called("auth_revoke_api_key")
    assert "Invalid key id." in client.get("/admin").get_data(as_text=True)


def test_admin_api_error_is_flashed(client, fake, admin_token, monkeypatch):
    _logged_in(client, fake)

    def refuse(cookie_header, key_id):
        raise apiclient.ApiError("planetGen API error (403): <credentials> stale", status_code=403)
    monkeypatch.setattr(apiclient, "auth_revoke_api_key", refuse)
    client.post("/admin", data={csrf.FIELD_NAME: admin_token, "action": "revoke_key", "key_id": "1"})
    html = client.get("/admin").get_data(as_text=True)
    assert "&lt;credentials&gt; stale" in html


def test_admin_post_requires_csrf(client, fake):
    _logged_in(client, fake)
    resp = client.post("/admin", data={"action": "create_key", "label": "x"})
    assert resp.status_code == 400
    assert not fake.called("auth_create_api_key")


def test_admin_page_no_longer_has_the_sector_wiki_form(client, fake):
    """UX.73: the form lives in the sector page's Admin menu."""
    _logged_in(client, fake)
    html = client.get("/admin").get_data(as_text=True)
    assert "Sector wiki link" not in html and 'value="set_sector_wiki_url"' not in html


def test_admin_hub_lists_the_admin_pages_without_a_signed_in_card(client, fake):
    """UX.72: the hub links the pages under it, and the tab row repeats them."""
    _logged_in(client, fake)
    html = client.get("/admin").get_data(as_text=True)
    assert "Signed in" not in html
    hub = re.search(r'<ul class="admin-hub">.*?</ul>', html, re.S).group(0)
    for path in ("/admin/generate", "/admin/queue", "/admin/stats", "/account"):
        assert f'href="{path}"' in hub
    tabs = re.search(r'<nav class="admin-tabs".*?</nav>', html, re.S).group(0)
    assert '<a href="/admin" aria-current="page">Overview</a>' in tabs
    assert all(f">{name}</a>" in tabs for name in ("Generate", "Queue", "Stats"))


def test_gear_menu_has_only_account_admin_and_logout(client, fake):
    _logged_in(client, fake)
    html = client.get("/admin").get_data(as_text=True)
    gear = re.search(r'<nav class="site-account".*?</nav>', html, re.S).group(0)
    assert re.findall(r">([A-Za-z: ]+)</(?:a|button)>", gear) == ["Theme: System", "Account", "Admin", "Logout"]


def test_admin_unknown_action(client, fake, admin_token):
    _logged_in(client, fake)
    resp = client.post("/admin", data={csrf.FIELD_NAME: admin_token, "action": "drop_tables"})
    assert resp.headers["Location"] == "/admin"
    assert "Unrecognized form action." in client.get("/admin").get_data(as_text=True)


# --- /admin/stats --------------------------------------------------------------------------

def test_admin_stats_renders(client, fake):
    _logged_in(client, fake)
    resp = client.get("/admin/stats?db=someone_elses_db")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert resp.headers["Cache-Control"] == "no-store"
    assert fake.called("admin_stats") == [("admin_stats", DB)]
    assert "someone_elses_db" not in html
    assert '<span class="badge">healthy</span>' in html
    assert "3,456" in html and "5.0 MB" in html and "1.03 hours" in html
    assert "0.10 / 0.20 / 0.30" in html
    # Duplicate names: links to the sector/system pages, escaped names.
    assert "Vega &lt;b&gt;Beta&lt;/b&gt;" in html
    assert 'href="/sector.py?db=' in html or 'href="/sector/5"' in html
    assert "Planet/moon collisions" not in html
    assert '<a href="/admin/stats" aria-current="page">Stats</a>' in html
    assert "/admin/stats/lockouts" not in re.search(r'<main.*</main>', html, re.S).group(0)


def test_admin_stats_shows_generation_speed(client, fake):
    """PERF.10: the measured speed per density bucket, and this galaxy's
    size per system."""
    _logged_in(client, fake)
    html = client.get("/admin/stats").get_data(as_text=True)
    assert "Nothing measured yet" in html
    fake.generation = {"available": True, "sizes": {DB: {"bytes_per_system": 62000.0, "systems": 1200,
                                                         "total_bytes": 74_400_000}},
                       "buckets": [{"kind": "sector", "bucket": 4, "density_low": 1.0, "density_high": 3.1623,
                                    "samples": 76, "seconds_per_task": 1.269, "seconds_per_system": 0.0551,
                                    "systems_per_task": 23.8, "stars_per_system": 1.3, "max_density": 3.0}]}
    html = client.get("/admin/stats").get_data(as_text=True)
    assert "60.5 KB per star system (1,200 systems, 71.0 MB)" in html
    assert "<td>Sector fill</td><td>1 to 3.16</td>" in html and "1.27 s" in html and "55.1 ms" in html
    fake.generation = {"available": False, "sizes": {}, "buckets": []}
    assert "run <code>update.sh</code>" in client.get("/admin/stats").get_data(as_text=True)


def test_admin_stats_counts_bright_stars(client, fake):
    _logged_in(client, fake)
    html = client.get("/admin/stats").get_data(as_text=True)
    assert re.search(r'<th scope="row">Bright stars</th><td>about 61,234 placed in all \(estimate\)</td>', html)


def test_admin_stats_exact_bright_star_counts(client, fake):
    _logged_in(client, fake)
    fake.database_extra = {"bright_stars": {"placed": 61234, "filled": 1234, "unfilled": 60000}}
    html = client.get("/admin/stats").get_data(as_text=True)
    assert re.search(r'<th scope="row">Bright stars</th><td>61,234 placed: 1,234 built into systems, '
                     r'60,000 waiting for their sectors</td>', html)


@pytest.mark.parametrize("database, text", [
    ({"tables": []}, "not tracked by this schema"),
    ({"tables": [{"name": "bright_stars", "approx_rows": 0}]}, "none pre-placed"),
    ({"tables": [], "bright_stars": {"placed": 0, "filled": 0, "unfilled": 0}}, "none pre-placed"),
    ({"tables": [], "bright_stars": None}, "not tracked by this schema"),
])
def test_bright_star_text_without_a_scatter(database, text):
    assert admin_pages.bright_star_text(database) == text


def test_admin_stats_sector_density_row(client, fake):
    # PERF.11: the decaying average of actual against expected systems.
    _logged_in(client, fake)
    fake.database_extra = {"sector_stats": {"measured": 12, "backfilled": 300, "ratio": 0.9731, "fills": 12}}
    html = client.get("/admin/stats").get_data(as_text=True)
    assert re.search(r'<th scope="row">Sector density</th><td>0.97 of the expected systems \(average of 12 '
                     r'fills\); 12 sectors measured, 300 sectors backfilled to their own level</td>', html)


@pytest.mark.parametrize("database, text", [
    ({}, "not tracked by this schema"),
    ({"sector_stats": {"measured": 0, "backfilled": 0, "ratio": None, "fills": 0}},
     "no sector filled yet; 0 sectors backfilled to their own level"),
])
def test_density_text_without_fills(database, text):
    assert admin_pages.density_text(database) == text


def test_admin_stats_names_paged(client, fake):
    _logged_in(client, fake)
    fake.names = [{"base_name": f"Name {i:03d}", "levels": ["sector"], "rows": []} for i in range(60)]
    html = client.get("/admin/stats?names_page=2").get_data(as_text=True)
    assert ("admin_duplicate_names", DB, 50, 50) in fake.calls
    assert "Name 050" in html and "Name 049" not in html
    assert 'href="/admin/stats?names_page=1#duplicate-names"' in html
    assert "no rows left" in html


def test_admin_names_table_route_pages_and_links_the_names(client, fake):
    assert client.get("/table/duplicate-names").status_code == 403
    _logged_in(client, fake)
    fake.names = [{"base_name": f"Name {i:03d}", "levels": ["sector", "system"],
                   "rows": [{"kind": "sector", "id": i, "name": f"Alpha Name {i:03d}"},
                            {"kind": "system", "id": i, "name": f"Beta Name {i:03d}"}]} for i in range(60)]
    data = client.get("/table/duplicate-names?offset=50&limit=50").get_json()
    assert data["total"] == 60 and len(data["rows"]) == 10
    base, levels, named = data["rows"][0]
    assert base["text"] == "Name 050" and "parts" in named
    links = [part for part in named["parts"] if isinstance(part, dict)]
    assert [link["text"] for link in links] == ["Alpha Name 050", "Beta Name 050"]
    assert links[0]["href"] == "/sector/50" and links[1]["href"] == "/system/50"


def test_admin_stats_database_unreachable(client, fake, monkeypatch):
    _logged_in(client, fake)
    monkeypatch.setattr(apiclient, "admin_stats", lambda cookie, db: _stats(reachable=False))
    html = client.get("/admin/stats").get_data(as_text=True)
    assert "database unreachable" in html
    assert "connection &lt;refused&gt;" in html
    assert "Names made unique" not in html


# --- Old CGI URLs ----------------------------------------------------------------------------

@pytest.mark.parametrize("url,location", [
    ("/login.py?db=x", "/login"),
    ("/logout.py", "/logout"),
    ("/changecreds.py", "/account"),
    ("/admin.py?keys_page=3", "/admin?keys_page=3"),
    ("/adminstats.py?db=x&names_page=2", "/admin/stats?names_page=2"),
])
def test_old_admin_urls_redirect(client, url, location):
    result = client.get(url)
    assert result.status_code == 301
    assert result.headers["Location"] == location


def test_old_logout_url_does_not_log_out(client):
    """A GET of the old logout.py only redirects; it never reaches the API."""
    client.set_cookie(SESSION_COOKIE_NAME, "abc")
    result = client.get("/logout.py")
    assert result.status_code == 301
    assert "Set-Cookie" not in result.headers


# --- Real database, in-process ----------------------------------------------------------------

@pytest.fixture
def db_app(mysql_config):
    class RealConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "test-secret"

    _username, first_password = adminAuth.bootstrap_control_schema(mysql_config)
    _db.get_connection(mysql_config).close()
    application = create_app(RealConfig)
    application.testing = True
    application.first_admin_password = first_password  # printed once by migrateDb on a real install
    return application


def _csrf_client(app):
    test_client = app.test_client()
    test_client.set_cookie(csrf.COOKIE_NAME, NONCE)
    with app.app_context():
        return test_client, csrf._sign(NONCE)  # not logged in yet: bound to no session


def test_real_login_account_admin_logout_flow(db_app, monkeypatch):
    def no_http(*args, **kwargs):
        raise AssertionError("HTTP transport used inside a Flask request")
    monkeypatch.setattr(apiclient, "_http_transport", no_http)
    client, token = _csrf_client(db_app)
    first = db_app.first_admin_password

    assert _redirect(client.get("/admin")) == ("/login", "/admin")

    resp = client.post("/login", data={csrf.FIELD_NAME: token, "username": "admin", "password": first,
                                       "next": "/admin/stats"})
    assert resp.status_code == 303
    assert _redirect(resp) == ("/account", "/admin/stats")
    session_cookie = [h for h in _set_cookies(resp) if h.startswith(f"{SESSION_COOKIE_NAME}=")][0]
    assert "HttpOnly" in session_cookie and "SameSite=Strict" in session_cookie and "Path=/" in session_cookie

    # Still on the first credentials: admin pages send us to /account.
    assert _redirect(client.get("/admin")) == ("/account", "/admin")
    # The login changed the session, so the pre-login token is refused and
    # the page now renders one for the new session (security #49).
    assert client.post("/account", data={csrf.FIELD_NAME: token, "current_password": first,
                                         "new_username": "boss", "new_password": "short"}).status_code == 400
    token = _form_token(client, "/account")
    bad = client.post("/account", data={csrf.FIELD_NAME: token, "current_password": first,
                                        "new_username": "boss", "new_password": "short"})
    assert bad.status_code == 200
    resp = client.post("/account", data={csrf.FIELD_NAME: token, "current_password": first,
                                         "new_username": "boss", "new_password": "a-much-longer-password",
                                         "next": "/admin/stats"})
    assert resp.status_code == 303
    assert resp.headers["Location"] == "/admin/stats"
    assert any(h.startswith(f"{SESSION_COOKIE_NAME}=") for h in _set_cookies(resp))

    stats = client.get("/admin/stats")
    assert stats.status_code == 200
    assert "Server health" in stats.get_data(as_text=True)

    token = _form_token(client, "/admin")  # the change re-issued the session
    resp = client.post("/admin", data={csrf.FIELD_NAME: token, "action": "create_key", "label": "ci"})
    assert resp.status_code == 303
    flash = [h for h in _set_cookies(resp) if h.startswith(f"{admin_pages.FLASH_COOKIE}=")][0]
    assert "Path=/admin" in flash  # the new key is never sent to the rest of the site
    assert client.get_cookie(admin_pages.FLASH_COOKIE, path="/admin") is not None
    html = client.get("/admin").get_data(as_text=True)
    key = re.search(r'<p class="api-key-value">([^<]+)</p>', html).group(1)
    assert client.get("/api/auth/me", headers={"Authorization": f"Bearer {key}"}).status_code == 200

    assert client.get("/logout").status_code == 200
    assert client.get("/admin").status_code == 200  # GET /logout changed nothing
    resp = client.post("/logout", data={csrf.FIELD_NAME: token})
    assert resp.status_code == 303
    assert any(h.startswith(f"{SESSION_COOKIE_NAME}=;") for h in _set_cookies(resp))
    assert client.get("/admin").status_code == 302


def test_real_login_keeps_rate_limit(db_app):
    from planetgen.api.auth import LOGIN_RATE_LIMIT
    from planetgen.api.limiter import limiter
    per_minute = int(LOGIN_RATE_LIMIT.split()[0])
    client, token = _csrf_client(db_app)
    limiter.reset()
    try:
        statuses = [client.post("/login", data={csrf.FIELD_NAME: token, "username": "admin",
                                                "password": "wrong"}).status_code
                    for _ in range(per_minute + 1)]
        assert statuses[:per_minute] == [401] * per_minute  # the form again, "Invalid username or password."
        assert statuses[-1] == 429
        last = client.post("/login", data={csrf.FIELD_NAME: token, "username": "admin",
                                           "password": db_app.first_admin_password})
        assert last.status_code == 429
        assert "Too many login attempts" in last.get_data(as_text=True)
    finally:
        limiter.reset()


def test_admin_stats_lists_and_lifts_lockouts(client, fake, admin_token):
    _logged_in(client, fake)
    fake.lockouts = [{"scope": "ip", "subject": "93.184.216.34", "retry_after": 290,
                      "locked_until": "2026-10-01T10:00:00Z", "level": 1}]
    fake.failures = [{"action": "login.failed", "username": "<b>x</b>", "ip": "93.184.216.34",
                      "created_at": "2026-10-01T09:55:00Z"}]
    html = client.get("/admin/stats").get_data(as_text=True)
    assert "<code>93.184.216.34</code>" in html and "4.83 minutes" in html
    assert "&lt;b&gt;x&lt;/b&gt;" in html and "Wrong username or password" in html
    resp = client.post("/admin/stats/lockouts", data={csrf.FIELD_NAME: admin_token, "scope": "ip",
                                                      "subject": "93.184.216.34"})
    assert resp.status_code == 303 and resp.headers["Location"] == "/admin/stats#lockouts"
    assert fake.called("admin_lift_lockout") == [("admin_lift_lockout", "ip", "93.184.216.34", False)]
    client.post("/admin/stats/lockouts", data={csrf.FIELD_NAME: admin_token, "all": "1"})
    assert fake.called("admin_lift_lockout")[-1] == ("admin_lift_lockout", None, None, True)


def test_admin_stats_shows_and_changes_the_naming_key(client, fake, admin_token):
    _logged_in(client, fake)
    html = client.get("/admin/stats").get_data(as_text=True)
    assert "<code>0A1B2C3D</code>" in html and "Naming key" in html
    resp = client.post("/admin/stats/naming-key", data={csrf.FIELD_NAME: admin_token, "key": "ffeeddcc"})
    assert resp.status_code == 303 and resp.headers["Location"] == "/admin/stats#naming-key"
    assert fake.called("admin_set_naming_key")[-1][2:] == ("ffeeddcc", False)
    client.post("/admin/stats/naming-key", data={csrf.FIELD_NAME: admin_token, "draw": "1"})
    assert fake.called("admin_set_naming_key")[-1][2:] == (None, True)


def test_naming_key_form_needs_an_admin(client, fake, token):
    resp = client.post("/admin/stats/naming-key", data={csrf.FIELD_NAME: token, "draw": "1"})
    assert resp.status_code == 302 and "/login" in resp.headers["Location"]
    assert not fake.called("admin_set_naming_key")


def test_lift_lockout_needs_an_admin(client, fake, token):
    resp = client.post("/admin/stats/lockouts", data={csrf.FIELD_NAME: token, "all": "1"})
    assert resp.status_code == 302 and "/login" in resp.headers["Location"]
    assert not fake.called("admin_lift_lockout")


def test_admin_stats_lists_and_downloads_galaxy_settings(client, fake):
    """ADM.18: the creation-settings file is listed and downloadable."""
    _logged_in(client, fake)
    name = fake.settings_files[0]["name"]
    html = client.get("/admin/stats").get_data(as_text=True)
    assert name in html and "Current" in html
    resp = client.get(f"/admin/stats/galaxy-settings/{name}")
    assert resp.status_code == 200
    assert f'attachment; filename="{name}"' in resp.headers["Content-Disposition"]
    assert resp.get_json()["seed"] == "AB"
    assert client.get("/admin/stats/galaxy-settings/nope.json").status_code == 404
