# tests/test_web_pages.py

"""
The Flask-served HTML pages (`src/html/web/`).

Most tests fake the data layer by monkeypatching the `apiclient`
functions the views call (the same functions the CGI pages use), so they
need no database. The tests at the bottom (marked by the `mysql_config`
fixture) run against a real throwaway database through the in-process
transport, and are skipped without a MySQL test server, like the rest of
the suite.
"""

import re

import pytest

from api.app import create_app
from api.authz import SESSION_COOKIE_NAME
from api.config import Config

import web  # noqa: F401 -- puts src/html/lib on sys.path
import apiclient  # noqa: E402
from stellarObjects import _db, adminAuth  # noqa: E402
from stellarObjects.spaceSector import SpaceSector  # noqa: E402
from web import csrf  # noqa: E402
from web.helpers import LEGACY_PAGES, page_url  # noqa: E402

from tests.webpage_support import run_page  # noqa: E402

DB = "planetgen_web_test"


class _FakeConfig(Config):
    WEB_DATABASE = DB
    SESSION_COOKIE_SECURE = False
    SECRET_KEY = "test-secret"


def _sector(i, name=None, placed=False):
    return {
        "id": i, "name": name or f"Sector {i:03d}", "system_count": i % 5, "edge_ly": 10.0,
        "placed": placed, "center_x_pc": 100.0 if placed else None, "center_y_pc": 50.0 if placed else None,
        "galactic_radius_ly": 1234.5 if placed else None,
    }


def _system(i, name=None):
    return {"id": 1000 + i, "name": name or f"System {i:03d}", "is_binary": i % 2, "star_summary": "G2V"}


class FakeData:
    """Stands in for the `apiclient` read functions; records calls."""

    def __init__(self, sectors=3, systems=2):
        self.sectors = [_sector(i) for i in range(sectors)]
        self.systems = [_system(i) for i in range(systems)]
        self.calls = []
        self.admin = None

    def get_sectors(self, db, limit=None, offset=None):
        self.calls.append(("get_sectors", db, limit, offset))
        return {"items": self.sectors[offset:offset + limit], "total": len(self.sectors),
                "limit": limit, "offset": offset}

    def get_systems(self, db, star_type=None, sector_id=None, limit=None, offset=None):
        self.calls.append(("get_systems", db, sector_id, limit, offset))
        return {"items": self.systems[offset:offset + limit], "total": len(self.systems),
                "limit": limit, "offset": offset}

    def auth_me(self, cookie_header):
        self.calls.append(("auth_me", cookie_header))
        return self.admin


@pytest.fixture
def fake(monkeypatch):
    data = FakeData()
    for name in ("get_sectors", "get_systems", "auth_me"):
        monkeypatch.setattr(apiclient, name, getattr(data, name))
    return data


@pytest.fixture
def app():
    application = create_app(_FakeConfig)
    application.testing = True
    return application


@pytest.fixture
def client(app):
    return app.test_client()


def _section_link(html, label):
    match = re.search(r'<nav class="site-sections" aria-label="Main">.*?</nav>', html, re.S)
    link = re.search(rf'<a href="[^"]*"( aria-current="page")?>{label}</a>', match.group(0))
    return link


# --- Rendering ------------------------------------------------------------------

def test_home_renders_both_tables_and_shell(client, fake):
    resp = client.get("/")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert resp.mimetype == "text/html"
    assert '<a class="skip-link" href="#main">' in html
    assert '<main id="main"' in html
    assert '<meta name="viewport" content="width=device-width, initial-scale=1">' in html
    assert '<meta name="description"' in html
    assert 'name="theme-color"' in html
    assert re.search(r'href="/static/favicon.svg\?v=[^"]+"', html)
    assert re.search(r'<script src="/static/theme.js\?v=[^"]+"></script>', html)
    assert re.search(r'<link rel="stylesheet" href="/static/style.css\?v=[^"]+">', html)
    assert "data-theme-toggle hidden" in html
    assert '<form class="site-search" role="search" method="get" action="/search">' in html
    assert "<details class=\"site-menu\">" in html
    assert "Sector 002" in html and "System 001" in html
    assert "3 sectors" in html and "2 standalone systems" in html
    # Not-yet-moved pages are plain GET links to the CGI scripts.
    assert f'href="/sector.py?db={DB}&amp;id=1"' in html
    assert f'href="/system.py?db={DB}&amp;id=1001"' in html
    # The home page is no section; nothing in the main nav is current.
    assert 'aria-current="page"' not in re.search(
        r'<nav class="site-sections".*?</nav>', html, re.S).group(0)
    # No breadcrumb trail on the home page.
    assert 'aria-label="Breadcrumb"' not in html


def test_home_security_headers(client, fake):
    resp = client.get("/")
    csp = resp.headers["Content-Security-Policy"]
    for directive in ("default-src 'self'", "base-uri 'self'", "form-action 'self'",
                      "frame-ancestors 'none'", "object-src 'none'"):
        assert directive in csp
    assert resp.headers["X-Frame-Options"] == "DENY"
    assert resp.headers["X-Content-Type-Options"] == "nosniff"
    assert resp.headers["Referrer-Policy"] == "no-referrer"


def test_api_keeps_its_own_csp(client):
    resp = client.get("/api/wiki-config")
    assert resp.status_code == 200
    assert resp.headers["Content-Security-Policy"] == "default-src 'none'"


@pytest.mark.parametrize("path,label", [("/sectors", "Sectors"), ("/systems", "Systems")])
def test_section_pages_mark_current_section_and_breadcrumbs(client, fake, path, label):
    html = client.get(path).get_data(as_text=True)
    assert _section_link(html, label).group(1) == ' aria-current="page"'
    for other in ("Galaxy", "Phenomena", "Nav", "Sectors", "Systems"):
        if other != label:
            assert _section_link(html, other).group(1) is None
    crumbs = re.search(r'<nav class="breadcrumbs" aria-label="Breadcrumb">(.*?)</nav>', html, re.S).group(1)
    assert '<li><a href="/">Home</a></li>' in crumbs
    assert f'<li><span aria-current="page">{label}</span></li>' in crumbs


def test_sections_link_to_moved_and_legacy_pages(client, fake):
    html = client.get("/").get_data(as_text=True)
    assert _section_link(html, "Sectors").group(0).startswith('<a href="/sectors"')
    assert _section_link(html, "Systems").group(0).startswith('<a href="/systems"')
    assert _section_link(html, "Galaxy").group(0).startswith('<a href="/galaxy"')
    assert _section_link(html, "Nav").group(0).startswith(f'<a href="/nav.py?db={DB}"')


def test_database_comes_from_config_never_the_url(client, fake):
    html = client.get("/?db=someone_elses_db").get_data(as_text=True)
    assert {call[1] for call in fake.calls if call[0].startswith("get_")} == {DB}
    assert "someone_elses_db" not in html
    # No moved-page URL carries the database.
    assert not re.search(r'href="/(sectors|systems)?\?[^"]*db=', html)


def test_database_defaults_to_mysql_config(fake):
    class DefaultDb(_FakeConfig):
        WEB_DATABASE = ""
    app = create_app(DefaultDb)
    app.testing = True
    app.test_client().get("/")
    assert {call[1] for call in fake.calls if call[0].startswith("get_")} == {Config.MYSQL_CONFIG.database}


def test_names_are_escaped(client, fake):
    fake.sectors = [_sector(1, name='<script>alert("x")</script>')]
    fake.systems = [_system(1, name="Tom & <b>Jerry</b>")]
    html = client.get("/").get_data(as_text=True)
    assert "<script>alert" not in html
    assert "&lt;script&gt;alert(&#34;x&#34;)&lt;/script&gt;" in html
    assert "Tom &amp; &lt;b&gt;Jerry&lt;/b&gt;" in html


def test_placed_sector_links_its_quadrant(client, fake):
    fake.sectors = [_sector(7, placed=True)]
    html = client.get("/").get_data(as_text=True)
    assert 'href="/galaxy?quadrant=I">Quadrant I</a>' in html
    assert "1,234.5 ly" in html


def test_pagination_uses_get_links_and_keeps_other_table_page(client, fake):
    fake.sectors = [_sector(i) for i in range(120)]
    fake.systems = [_system(i) for i in range(120)]
    html = client.get("/?sectors_page=2&standalone_page=3").get_data(as_text=True)
    assert ("get_sectors", DB, 50, 50) in fake.calls
    assert ("get_systems", DB, "none", 50, 100) in fake.calls
    assert "Sector 050" in html and "Sector 049" not in html
    assert "System 100" in html
    # The sectors pager keeps standalone_page=3 and vice versa; no forms.
    assert 'href="/?standalone_page=3&amp;sectors_page=3#sectors"' in html
    assert 'href="/?sectors_page=2&amp;standalone_page=2#standalone-systems"' in html
    assert '<form method="post"' not in html


def test_page_past_the_end_shows_last_page(client, fake):
    fake.sectors = [_sector(i) for i in range(60)]
    html = client.get("/sectors?sectors_page=99").get_data(as_text=True)
    assert "Showing 51&ndash;60 of 60" in html


def test_search_box_forwards_to_cgi_search(client, fake):
    resp = client.get("/search?q=Kepler 42")
    assert resp.status_code == 302
    assert resp.headers["Location"] == f"/search.py?db={DB}&system_q=Kepler+42"


# --- Account menu, one login check per request -------------------------------------

def test_login_link_and_no_auth_lookup_without_cookie(client, fake):
    html = client.get("/").get_data(as_text=True)
    assert f'<a href="/login.py?db={DB}"' in html
    assert ">Admin</a>" not in html
    assert not [call for call in fake.calls if call[0] == "auth_me"]


def test_admin_menu_with_one_lookup_per_request(client, fake):
    fake.admin = {"username": "admin", "must_change_credentials": False}
    client.set_cookie(SESSION_COOKIE_NAME, "token-value")
    html = client.get("/").get_data(as_text=True)
    assert ">Admin</a>" in html and ">Stats</a>" in html and ">Logout</a>" in html
    assert ">Login</a>" not in html
    lookups = [call for call in fake.calls if call[0] == "auth_me"]
    assert len(lookups) == 1  # both menus (wide and narrow) share it
    assert f"{SESSION_COOKIE_NAME}=token-value" in lookups[0][1]


def test_failed_login_lookup_counts_as_logged_out(client, fake, monkeypatch):
    def broken(cookie_header):
        raise apiclient.ApiError("down", status_code=503)
    monkeypatch.setattr(apiclient, "auth_me", broken)
    client.set_cookie(SESSION_COOKIE_NAME, "token-value")
    resp = client.get("/")
    assert resp.status_code == 200
    assert ">Login</a>" in resp.get_data(as_text=True)


# --- Errors -------------------------------------------------------------------------

def test_unknown_page_is_html_404(client):
    resp = client.get("/no/such/page")
    assert resp.status_code == 404
    assert resp.mimetype == "text/html"
    assert "There is no page at this address." in resp.get_data(as_text=True)


def test_unknown_api_path_stays_json_404(client):
    resp = client.get("/api/no-such-endpoint")
    assert resp.status_code == 404
    assert resp.get_json() == {"error": "not found"}


def test_not_found_from_data_layer_is_404_page(client, fake, monkeypatch):
    def missing(*args, **kwargs):
        raise apiclient.NotFoundError("Unknown database <x>")
    monkeypatch.setattr(apiclient, "get_sectors", missing)
    resp = client.get("/")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 404
    assert "Unknown database &lt;x&gt;" in html


def test_api_error_is_502_without_detail(client, fake, monkeypatch):
    def failing(*args, **kwargs):
        raise apiclient.ApiError("planetGen API error (500): SELECT secret FROM /srv/app.py")
    monkeypatch.setattr(apiclient, "get_sectors", failing)
    resp = client.get("/")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 502
    assert "could not be loaded" in html
    assert "SELECT secret" not in html and "/srv/app.py" not in html


def test_unexpected_exception_is_500_without_traceback(app, fake, monkeypatch):
    app.testing = False  # let the error handler run instead of re-raising
    monkeypatch.delenv("PLANETGEN_DEBUG", raising=False)

    def boom(*args, **kwargs):
        raise RuntimeError("secret internal detail at /srv/planetgen/src/x.py")
    monkeypatch.setattr(apiclient, "get_sectors", boom)
    resp = app.test_client().get("/")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 500
    assert resp.mimetype == "text/html"
    assert "An unexpected error occurred." in html
    assert "Traceback" not in html and "secret internal detail" not in html


# --- CSRF ---------------------------------------------------------------------------

def _token_for(app, nonce):
    with app.app_context():
        return csrf._sign(nonce)


def test_post_without_token_is_rejected(client):
    resp = client.post("/", data={"x": "1"})
    assert resp.status_code == 400
    assert "expired or was not sent from this site" in resp.get_data(as_text=True)


def test_post_with_wrong_token_is_rejected(app, client):
    nonce = "n" * 43
    client.set_cookie(csrf.COOKIE_NAME, nonce)
    resp = client.post("/", data={csrf.FIELD_NAME: _token_for(app, "other" * 9)})
    assert resp.status_code == 400


def test_post_with_token_but_no_cookie_is_rejected(app, client):
    resp = client.post("/", data={csrf.FIELD_NAME: _token_for(app, "n" * 43)})
    assert resp.status_code == 400


def test_post_with_valid_token_passes_the_check(app, client):
    nonce = "n" * 43
    client.set_cookie(csrf.COOKIE_NAME, nonce)
    resp = client.post("/", data={csrf.FIELD_NAME: _token_for(app, nonce)})
    # Past the CSRF check; "/" itself only answers GET.
    assert resp.status_code == 405
    assert resp.mimetype == "text/html"


def test_api_posts_are_not_subject_to_the_form_token(client):
    resp = client.post("/api/no-such-endpoint", json={})
    assert resp.status_code == 404
    assert resp.get_json() == {"error": "not found"}


def test_csrf_field_sets_cookie_once(app):
    with app.test_request_context("/"):
        field = str(csrf.csrf_field())
        token = re.search(r'value="([0-9a-f]{64})"', field).group(1)
        response = app.make_response("ok")
        response = csrf.set_cookie(response)
        cookie = response.headers["Set-Cookie"]
        assert cookie.startswith(f"{csrf.COOKIE_NAME}=")
        assert "HttpOnly" in cookie and "SameSite=Strict" in cookie
        nonce = cookie.split(";")[0].split("=", 1)[1]
        assert csrf.valid(token, nonce)
    with app.test_request_context("/", headers={"Cookie": f"{csrf.COOKIE_NAME}={nonce}"}):
        assert str(csrf.csrf_field()) == field
        assert "Set-Cookie" not in csrf.set_cookie(app.make_response("ok")).headers


# --- Helpers ------------------------------------------------------------------------

def test_page_url_resolves_moved_and_legacy_pages(app):
    with app.test_request_context("/"):
        assert page_url("index") == "/"
        assert page_url("sectors") == "/sectors"
        assert page_url("sector", sector_id=5) == f"/sector.py?db={DB}&id=5"
        assert page_url("phenomenon", phenomenon_type="nebula", phenomenon_id=3) == \
            f"/phenomenon.py?db={DB}&type=nebula&id=3"
        assert page_url("galaxy", _anchor="map") == "/galaxy#map"
        with pytest.raises(KeyError):
            page_url("no_such_page")


def test_legacy_pages_are_not_also_routes(app):
    """A page PR that adds `web.<name>` must delete its LEGACY_PAGES entry."""
    for name in LEGACY_PAGES:
        assert f"web.{name}" not in app.view_functions, name


# --- In-process transport ---------------------------------------------------------------

def test_in_process_transport_skips_http_inside_a_request(app, monkeypatch):
    def no_http(*args, **kwargs):
        raise AssertionError("HTTP transport used inside a Flask request")
    monkeypatch.setattr(apiclient, "_http_transport", no_http)
    with app.test_request_context("/"):
        assert set(apiclient.get_wiki_config()) == {"wikijs", "mediawiki"}
        with pytest.raises(apiclient.NotFoundError):
            apiclient._request("/no-such-endpoint")


def test_in_process_transport_declines_outside_a_request(app, monkeypatch):
    def unreachable(*args, **kwargs):
        raise apiclient.TransportUnreachable("no server in tests")
    monkeypatch.setattr(apiclient, "_http_transport", unreachable)
    with pytest.raises(apiclient.ApiError, match="Could not reach"):
        apiclient.get_wiki_config()


# --- CGI shims ------------------------------------------------------------------------

@pytest.mark.parametrize("script", ["index.py", "browse.py"])
def test_cgi_shim_redirects_get(script):
    result = run_page("http://127.0.0.1:9/api", script,
                      query={"db": "x", "sectors_page": "2", "standalone_page": "junk"})
    assert result.status_code == 301
    assert result.headers["Location"] == "/?sectors_page=2"


@pytest.mark.parametrize("script", ["index.py", "browse.py"])
def test_cgi_shim_redirects_post(script):
    result = run_page("http://127.0.0.1:9/api", script, method="POST",
                      body={"db": "x", "standalone_page": "4"})
    assert result.status_code == 301
    assert result.headers["Location"] == "/?standalone_page=4"


def test_cgi_shim_without_params_goes_home():
    result = run_page("http://127.0.0.1:9/api", "browse.py")
    assert result.status_code == 301
    assert result.headers["Location"] == "/"


# --- Real database, in-process ----------------------------------------------------------

@pytest.fixture
def db_client(mysql_config):
    class RealConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "test-secret"

    application = create_app(RealConfig)
    application.testing = True
    return application.test_client()


def test_real_database_home_and_paging(db_client, mysql_config, monkeypatch):
    def no_http(*args, **kwargs):
        raise AssertionError("HTTP transport used inside a Flask request")
    monkeypatch.setattr(apiclient, "_http_transport", no_http)
    for i in range(55):
        _db.save_sector(SpaceSector(f"Web <i>Sector</i> {i:03d}", edge_ly=10.0), config=mysql_config)

    page1 = db_client.get("/")
    html = page1.get_data(as_text=True)
    assert page1.status_code == 200
    assert "55 sectors" in html
    assert "Web &lt;i&gt;Sector&lt;/i&gt; 000" in html
    assert "Web &lt;i&gt;Sector&lt;/i&gt; 050" not in html
    assert f"db={mysql_config.database}" in html  # legacy links only

    page2 = db_client.get("/sectors?sectors_page=2").get_data(as_text=True)
    assert "Web &lt;i&gt;Sector&lt;/i&gt; 054" in page2
    assert "Showing 51&ndash;55 of 55" in page2


def test_real_login_session_shows_admin_menu(db_client, mysql_config):
    adminAuth.bootstrap_control_schema(mysql_config)
    _db.get_connection(mysql_config).close()
    resp = db_client.post("/api/auth/login", json={"username": "admin", "password": "password"})
    assert resp.status_code == 200
    html = db_client.get("/").get_data(as_text=True)
    assert ">Admin</a>" in html and ">Logout</a>" in html

    db_client.delete_cookie(SESSION_COOKIE_NAME)
    assert ">Login</a>" in db_client.get("/").get_data(as_text=True)


def test_page_views_are_not_rate_limited(db_client, mysql_config):
    """
    Each home page view makes two in-process API calls. 60 views (120
    calls) is well past the default 50/hour, so any of them counting
    against the default limit would turn a page into a 502/429. A
    direct API call from the same address still counts.
    """
    from api.limiter import limiter
    _db.get_connection(mysql_config).close()  # lay down the (empty) content schema
    limiter.reset()
    try:
        statuses = {db_client.get("/").status_code for _ in range(60)}
        assert statuses == {200}
        direct = [db_client.get("/api/sectors").status_code for _ in range(51)]
        assert direct[-1] == 429
    finally:
        limiter.reset()
