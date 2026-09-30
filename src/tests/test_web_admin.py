# tests/test_web_admin.py

"""
The admin pages on the Flask app (`src/html/web/admin_pages.py`):
`/login`, `/logout`, `/account`, `/admin` and `/admin/stats`, which
replaced `login.py`, `logout.py`, `changecreds.py`, `admin.py` and
`adminstats.py`.

Most tests fake the data layer by monkeypatching the `apiclient`
functions the views call (the pattern of `test_web_pages.py`). The tests
at the bottom run the real login flow against a throwaway database
through the in-process transport, and are skipped without a MySQL test
server.
"""

import re
from urllib.parse import parse_qs, urlsplit

import pytest

from api.app import create_app
from api.authz import SESSION_COOKIE_NAME
from api.config import Config

import web  # noqa: F401 -- puts src/html/lib on sys.path
import apiclient  # noqa: E402
from stellarObjects import _db, adminAuth  # noqa: E402
from web import admin_pages, csrf  # noqa: E402
from web.helpers import LEGACY_PAGES  # noqa: E402

from tests.webpage_support import run_page  # noqa: E402

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
        "name_collisions": {"distinct_base_names": 2, "sector": 1, "system": 1, "body": 0},
        "timestamps": [{"table": "sectors", "newest_created_at": "2026-09-01 10:00:00",
                        "last_modified_at": "2026-09-02 10:00:00"},
                       {"table": "star_systems", "newest_created_at": None, "last_modified_at": None}],
        "tables": [{"name": "star_systems", "approx_rows": 3456, "data_bytes": 2048, "index_bytes": 1024}],
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
        self.keys = [_key(1), _key(2, revoked=True)]
        self.login_error = None
        self.change_error = None
        self.names = [{"base_name": "Vega", "levels": ["system"], "rows": [
            {"kind": "sector", "id": 5, "name": "Vega Alpha"},
            {"kind": "system", "id": 77, "name": "Vega <b>Beta</b>"},
            {"kind": "planet", "id": 9, "name": "Vega Kin", "star_system_id": 77, "system_name": "Vega Beta"},
        ]}]

    def auth_me(self, cookie_header):
        self.calls.append(("auth_me", cookie_header))
        return self.admin

    def auth_login(self, username, password):
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
        return _stats()

    def admin_duplicate_names(self, cookie_header, db, limit=None, offset=None):
        self.calls.append(("admin_duplicate_names", db, limit, offset))
        return {"items": self.names[offset:offset + limit], "total": len(self.names),
                "limit": limit, "offset": offset}

    def called(self, name):
        return [call for call in self.calls if call[0] == name]


_FAKED = ("auth_me", "auth_login", "auth_logout", "auth_change_credentials", "auth_list_api_keys",
          "auth_create_api_key", "auth_revoke_api_key", "admin_set_sector_wiki_url", "admin_stats",
          "admin_duplicate_names")


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
    with app.app_context():
        return csrf._sign(NONCE)


def _logged_in(client, fake, must_change=False):
    fake.admin = {"username": "boss", "must_change_credentials": must_change}
    client.set_cookie(SESSION_COOKIE_NAME, "token-value")


def _set_cookies(resp):
    return resp.headers.getlist("Set-Cookie")


def _redirect(resp):
    """`(path, next)` of a redirect's Location."""
    parts = urlsplit(resp.headers["Location"])
    return parts.path, parse_qs(parts.query).get("next", [None])[0]


# --- Routing -------------------------------------------------------------------------

def test_admin_pages_are_routes_not_legacy(app):
    for name in ("login", "logout", "account", "admin", "admin_stats"):
        assert f"web.{name}" in app.view_functions
        assert name not in LEGACY_PAGES


def test_header_links_point_at_new_urls(client, fake):
    html = client.get("/login").get_data(as_text=True)
    assert '<a href="/login" aria-current="page">Login</a>' in html
    _logged_in(client, fake)
    html = client.get("/admin").get_data(as_text=True)
    assert '<a href="/admin" aria-current="page">Admin</a>' in html
    assert '<a href="/admin/stats">Stats</a>' in html
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
    assert resp.status_code == 200
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


def test_logout_post_clears_cookie(client, fake, token):
    _logged_in(client, fake)
    resp = client.post("/logout", data={csrf.FIELD_NAME: token})
    assert resp.status_code == 303
    assert resp.headers["Location"] == "/"
    assert SESSION_CLEAR in _set_cookies(resp)
    assert f"{SESSION_COOKIE_NAME}=token-value" in fake.called("auth_logout")[0][1]


# --- /account ------------------------------------------------------------------------------

def test_account_form(client, fake):
    _logged_in(client, fake, must_change=True)
    html = client.get("/account?next=/admin/stats").get_data(as_text=True)
    assert "still uses the default username and password" in html
    assert 'name="new_username" value="boss"' in html
    assert '<input type="hidden" name="next" value="/admin/stats">' in html
    assert f'name="{csrf.FIELD_NAME}"' in html


def test_account_change_relays_new_cookie(client, fake, token):
    _logged_in(client, fake, must_change=True)
    resp = client.post("/account", data={csrf.FIELD_NAME: token, "current_password": "password",
                                         "new_username": "boss2", "new_password": "a-long-new-password",
                                         "next": "/admin/stats"})
    assert resp.status_code == 303
    assert resp.headers["Location"] == "/admin/stats"
    assert SESSION_SET in _set_cookies(resp)
    call = fake.called("auth_change_credentials")[0]
    assert call[2:] == ("password", "boss2", "a-long-new-password")


def test_account_change_to_admin_flashes_a_message(client, fake, token):
    _logged_in(client, fake)
    resp = client.post("/account", data={csrf.FIELD_NAME: token, "current_password": "x",
                                         "new_username": "boss", "new_password": "a-long-new-password"})
    assert resp.headers["Location"] == "/admin"
    html = client.get("/admin").get_data(as_text=True)
    assert "Your username and password were changed." in html


def test_account_policy_error_shows_inline(client, fake, token):
    _logged_in(client, fake)
    fake.change_error = apiclient.ApiError("planetGen API error (400): password is too short", status_code=400)
    resp = client.post("/account", data={csrf.FIELD_NAME: token, "current_password": "x",
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
    assert "revoked 2026-02-01 00:00:00" in html
    # One revoke form (the active key), carrying the CSRF token.
    assert html.count('value="revoke_key"') == 1
    revoke = re.search(r'<form method="post" action="/admin" class="table-form">.*?</form>', html, re.S).group(0)
    assert f'name="{csrf.FIELD_NAME}"' in revoke and 'name="key_id" value="1"' in revoke
    assert f"<code>{DB}</code>" in html
    assert 'name="db"' not in html


def test_admin_keys_are_paged_with_get_links(client, fake):
    _logged_in(client, fake)
    fake.keys = [_key(i) for i in range(1, 61)]
    html = client.get("/admin?keys_page=2").get_data(as_text=True)
    assert "key 051" in html and "key 050" not in html
    assert 'href="/admin?keys_page=1#api-keys"' in html
    assert '<input type="hidden" name="keys_page" value="2">' in html


def test_admin_create_key_is_post_redirect_get_and_shown_once(client, fake, token):
    _logged_in(client, fake)
    resp = client.post("/admin", data={csrf.FIELD_NAME: token, "action": "create_key", "label": "my script"})
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


def test_admin_forged_flash_is_ignored(client, fake):
    _logged_in(client, fake)
    client.set_cookie(admin_pages.FLASH_COOKIE, "eyJuZXdfa2V5Ijp7ImtleSI6ImZvcmdlZCJ9fQ.bad.sig")
    html = client.get("/admin").get_data(as_text=True)
    assert "forged" not in html and "New API key" not in html


def test_admin_create_key_needs_label(client, fake, token):
    _logged_in(client, fake)
    resp = client.post("/admin", data={csrf.FIELD_NAME: token, "action": "create_key", "label": "  "})
    assert resp.headers["Location"] == "/admin#api-keys"
    assert not fake.called("auth_create_api_key")
    assert "Label is required." in client.get("/admin").get_data(as_text=True)


def test_admin_revoke_key_returns_to_same_page(client, fake, token):
    _logged_in(client, fake)
    resp = client.post("/admin", data={csrf.FIELD_NAME: token, "action": "revoke_key", "key_id": "1",
                                       "keys_page": "2"})
    assert resp.status_code == 303
    assert resp.headers["Location"] == "/admin?keys_page=2#api-keys"
    assert fake.called("auth_revoke_api_key")[0][2] == 1


def test_admin_revoke_bad_id(client, fake, token):
    _logged_in(client, fake)
    client.post("/admin", data={csrf.FIELD_NAME: token, "action": "revoke_key", "key_id": "x"})
    assert not fake.called("auth_revoke_api_key")
    assert "Invalid key id." in client.get("/admin").get_data(as_text=True)


def test_admin_api_error_is_flashed(client, fake, token, monkeypatch):
    _logged_in(client, fake)

    def refuse(cookie_header, key_id):
        raise apiclient.ApiError("planetGen API error (403): <credentials> stale", status_code=403)
    monkeypatch.setattr(apiclient, "auth_revoke_api_key", refuse)
    client.post("/admin", data={csrf.FIELD_NAME: token, "action": "revoke_key", "key_id": "1"})
    html = client.get("/admin").get_data(as_text=True)
    assert "&lt;credentials&gt; stale" in html


def test_admin_post_requires_csrf(client, fake):
    _logged_in(client, fake)
    resp = client.post("/admin", data={"action": "create_key", "label": "x"})
    assert resp.status_code == 400
    assert not fake.called("auth_create_api_key")


def test_admin_sets_wiki_url_in_configured_database(client, fake, token):
    _logged_in(client, fake)
    resp = client.post("/admin", data={csrf.FIELD_NAME: token, "action": "set_sector_wiki_url", "sector_id": "5",
                                       "wiki_url": "https://wiki.example/S5", "db": "someone_elses_db"})
    assert resp.headers["Location"] == "/admin#sector-wiki-link"
    assert fake.called("admin_set_sector_wiki_url") == [
        ("admin_set_sector_wiki_url", DB, 5, "https://wiki.example/S5")]
    assert "Wiki link for sector 5 set to https://wiki.example/S5." in client.get("/admin").get_data(as_text=True)


def test_admin_wiki_url_unknown_sector(client, fake, token, monkeypatch):
    _logged_in(client, fake)

    def missing(*args):
        raise apiclient.NotFoundError("Unknown sector 404")
    monkeypatch.setattr(apiclient, "admin_set_sector_wiki_url", missing)
    client.post("/admin", data={csrf.FIELD_NAME: token, "action": "set_sector_wiki_url", "sector_id": "404"})
    assert "Unknown sector 404" in client.get("/admin").get_data(as_text=True)


def test_admin_unknown_action(client, fake, token):
    _logged_in(client, fake)
    resp = client.post("/admin", data={csrf.FIELD_NAME: token, "action": "drop_tables"})
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
    assert "3,456" in html and "5.0 MB" in html and "1h 1m" in html
    assert "0.10 / 0.20 / 0.30" in html
    # Duplicate names: links to the sector/system pages, escaped names.
    assert "Vega &lt;b&gt;Beta&lt;/b&gt;" in html
    assert 'href="/sector.py?db=' in html or 'href="/sector/5"' in html
    assert "planet in" in html
    assert '<a href="/admin/stats" aria-current="page">Stats</a>' in html
    assert "<form" not in re.search(r'<main.*</main>', html, re.S).group(0)


def test_admin_stats_names_paged(client, fake):
    _logged_in(client, fake)
    fake.names = [{"base_name": f"Name {i:03d}", "levels": ["sector"], "rows": []} for i in range(60)]
    html = client.get("/admin/stats?names_page=2").get_data(as_text=True)
    assert ("admin_duplicate_names", DB, 50, 50) in fake.calls
    assert "Name 050" in html and "Name 049" not in html
    assert 'href="/admin/stats?names_page=1#duplicate-names"' in html
    assert "no rows left" in html


def test_admin_stats_database_unreachable(client, fake, monkeypatch):
    _logged_in(client, fake)
    monkeypatch.setattr(apiclient, "admin_stats", lambda cookie, db: _stats(reachable=False))
    html = client.get("/admin/stats").get_data(as_text=True)
    assert "database unreachable" in html
    assert "connection &lt;refused&gt;" in html
    assert "Names made unique" not in html


# --- CGI shims ---------------------------------------------------------------------------

@pytest.mark.parametrize("script,query,location", [
    ("login.py", {"db": "x"}, "/login"),
    ("logout.py", {}, "/logout"),
    ("changecreds.py", {}, "/account"),
    ("admin.py", {"keys_page": "3"}, "/admin?keys_page=3"),
    ("adminstats.py", {"db": "x", "names_page": "2"}, "/admin/stats?names_page=2"),
])
def test_cgi_shims_redirect(script, query, location):
    result = run_page("http://127.0.0.1:9/api", script, query=query)
    assert result.status_code == 301
    assert result.headers["Location"] == location


def test_cgi_logout_shim_does_not_log_out_on_post():
    """A POST to the old logout.py only redirects; it never reaches the API."""
    result = run_page("http://127.0.0.1:9/api", "logout.py", method="POST", body={"x": "1"},
                      cookie=f"{SESSION_COOKIE_NAME}=abc")
    assert result.status_code == 301
    assert result.headers["Location"] == "/logout"
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

    adminAuth.bootstrap_control_schema(mysql_config)
    _db.get_connection(mysql_config).close()
    application = create_app(RealConfig)
    application.testing = True
    return application


def _csrf_client(app):
    test_client = app.test_client()
    test_client.set_cookie(csrf.COOKIE_NAME, NONCE)
    with app.app_context():
        return test_client, csrf._sign(NONCE)


def test_real_login_account_admin_logout_flow(db_app, monkeypatch):
    def no_http(*args, **kwargs):
        raise AssertionError("HTTP transport used inside a Flask request")
    monkeypatch.setattr(apiclient, "_http_transport", no_http)
    client, token = _csrf_client(db_app)

    assert _redirect(client.get("/admin")) == ("/login", "/admin")

    resp = client.post("/login", data={csrf.FIELD_NAME: token, "username": "admin", "password": "password",
                                       "next": "/admin/stats"})
    assert resp.status_code == 303
    assert _redirect(resp) == ("/account", "/admin/stats")
    session_cookie = [h for h in _set_cookies(resp) if h.startswith(f"{SESSION_COOKIE_NAME}=")][0]
    assert "HttpOnly" in session_cookie and "SameSite=Strict" in session_cookie and "Path=/" in session_cookie

    # Still on the default credentials: admin pages send us to /account.
    assert _redirect(client.get("/admin")) == ("/account", "/admin")
    bad = client.post("/account", data={csrf.FIELD_NAME: token, "current_password": "password",
                                        "new_username": "boss", "new_password": "short"})
    assert bad.status_code == 200
    resp = client.post("/account", data={csrf.FIELD_NAME: token, "current_password": "password",
                                         "new_username": "boss", "new_password": "a-much-longer-password",
                                         "next": "/admin/stats"})
    assert resp.status_code == 303
    assert resp.headers["Location"] == "/admin/stats"
    assert any(h.startswith(f"{SESSION_COOKIE_NAME}=") for h in _set_cookies(resp))

    stats = client.get("/admin/stats")
    assert stats.status_code == 200
    assert "Server health" in stats.get_data(as_text=True)

    resp = client.post("/admin", data={csrf.FIELD_NAME: token, "action": "create_key", "label": "ci"})
    assert resp.status_code == 303
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
    from api.auth import LOGIN_RATE_LIMIT
    from api.limiter import limiter
    per_minute = int(LOGIN_RATE_LIMIT.split()[0])
    client, token = _csrf_client(db_app)
    limiter.reset()
    try:
        statuses = [client.post("/login", data={csrf.FIELD_NAME: token, "username": "admin",
                                                "password": "wrong"}).status_code
                    for _ in range(per_minute + 1)]
        assert statuses[:per_minute] == [200] * per_minute  # the form again, "Invalid username or password."
        assert statuses[-1] == 429
        last = client.post("/login", data={csrf.FIELD_NAME: token, "username": "admin", "password": "password"})
        assert last.status_code == 429
        assert "Too many login attempts" in last.get_data(as_text=True)
    finally:
        limiter.reset()
