# tests/test_web_pages.py

"""
The Flask-served HTML pages (`src/planetgen/web/`).

Most tests fake the data layer by monkeypatching the `apiclient`
functions the views call, so they need no database. The tests at the bottom (marked by the `mysql_config`
fixture) run against a real throwaway database through the in-process
transport, and are skipped without a MySQL test server, like the rest of
the suite.
"""

import re

import pytest

from planetgen.web.app import create_app
from planetgen.api.authz import SESSION_COOKIE_NAME
from planetgen.api.config import Config

from planetgen.web.lib import apiclient  # noqa: E402
from planetgen.db import store  # noqa: E402
from planetgen.admin import auth as adminAuth
from planetgen.galaxy.sector import SpaceSector  # noqa: E402
from planetgen.web import csrf  # noqa: E402
from planetgen.web.helpers import page_url  # noqa: E402
from werkzeug.routing import BuildError  # noqa: E402
from planetgen.web.old_urls import OLD_PAGES  # noqa: E402

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
        self.uncharted = []
        self.calls = []
        self.asked = []
        self.admin = None

    def get_sectors(self, db, limit=None, offset=None, sort=None, descending=False, quadrants=(), facets=False):
        self.calls.append(("get_sectors", db, limit, offset))
        self.asked.append(("sectors", sort, descending, list(quadrants)))
        body = {"items": self.sectors[offset:offset + limit], "total": len(self.sectors),
                "limit": limit, "offset": offset}
        if facets:
            body["facets"] = {"quadrant": [{"value": "I", "count": len(self.sectors)}]}
        return body

    def get_systems(self, db, star_type=None, sector_id=None, limit=None, offset=None, sort=None,
                    descending=False, binary=None, placement=None, octants=(), facets=False):
        self.calls.append(("get_systems", db, sector_id, limit, offset))
        self.asked.append(("systems", sort, descending, binary, placement, list(octants)))
        body = {"items": self.systems[offset:offset + limit], "total": len(self.systems),
                "limit": limit, "offset": offset}
        if facets:
            body["facets"] = {"placement": [{"value": "sector", "count": 1}],
                              "binary": [{"value": "yes", "count": 1}, {"value": "no", "count": 1}],
                              "octant": [{"value": "I", "count": 2}]}
        return body

    def get_uncharted_systems(self, db, limit=None, offset=None, sort=None, descending=False):
        self.calls.append(("get_uncharted_systems", db, limit, offset))
        self.asked.append(("uncharted", sort, descending))
        return {"items": self.uncharted[offset or 0:(offset or 0) + (limit or 50)], "total": len(self.uncharted),
                "limit": limit, "offset": offset}

    def get_phenomena(self, db, limit=None, offset=None, **kwargs):
        self.calls.append(("get_phenomena", db, limit, offset))
        return {"items": [], "total": 4, "limit": limit, "offset": offset}

    def auth_me(self, cookie_header):
        self.calls.append(("auth_me", cookie_header))
        return self.admin


@pytest.fixture
def fake(monkeypatch):
    data = FakeData()
    for name in ("get_sectors", "get_systems", "get_uncharted_systems", "get_phenomena", "auth_me"):
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

def test_home_renders_the_front_door_and_shell(client, fake):
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
    assert '<sl-dropdown class="site-menu"' in html
    assert re.search(r'<script type="module" src="/static/components.js\?v=[^"]+"></script>', html)
    # UX.55: the counts are links to their lists, and the tools are entry points.
    assert '<a href="/sectors"><strong>3</strong> sectors</a>' in html
    assert '<a href="/systems"><strong>2</strong> systems</a>' in html
    assert '<a href="/phenomena"><strong>4</strong> phenomena</a>' in html
    assert re.search(r'<a href="/systems\?placement=standalone#all-systems"><strong>2</strong> standalone systems</a>', html)
    for path in ("/galaxy", "/nav", "/search"):
        assert f'<a href="{path}">' in html
    assert "Sector 002" not in html
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


def test_static_scripts_get_the_page_csp(client):
    """A script can be started as a Web Worker, which runs under its own
    response's CSP (the Galaxy Map's galaxyblocks.js imports
    galaxyprisms.js), so scripts get the pages' policy, not 'none'."""
    from planetgen.web import CONTENT_SECURITY_POLICY

    resp = client.get("/static/galaxyblocks.js")
    assert resp.status_code == 200
    assert resp.headers["Content-Security-Policy"] == CONTENT_SECURITY_POLICY
    resp.close()
    resp = client.get("/static/style.css")
    assert resp.headers["Content-Security-Policy"] == "default-src 'none'"
    resp.close()


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
    assert _section_link(html, "Nav").group(0).startswith('<a href="/nav"')


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
    html = client.get("/sectors").get_data(as_text=True)
    assert "<script>alert" not in html
    assert "&lt;script&gt;alert(&#34;x&#34;)&lt;/script&gt;" in html
    assert "Tom &amp; &lt;b&gt;Jerry&lt;/b&gt;" in client.get("/systems").get_data(as_text=True)


def test_placed_sector_links_its_quadrant(client, fake):
    fake.sectors = [_sector(7, placed=True)]
    html = client.get("/sectors").get_data(as_text=True)
    assert 'href="/galaxy?quadrant=I">Quadrant I</a>' in html
    assert "378 pc (1,235 ly)" in html


def test_pagination_uses_get_links(client, fake):
    fake.sectors = [_sector(i) for i in range(120)]
    html = client.get("/sectors?sectors_page=2").get_data(as_text=True)
    assert ("get_sectors", DB, 50, 50) in fake.calls
    assert "Sector 050" in html and "Sector 049" not in html
    assert 'href="/sectors?sectors_page=3#sectors"' in html
    assert '<form method="post"' not in html


def test_standalone_systems_are_a_filter_on_the_systems_list(client, fake):
    """UX.55: no standalone card; the badge links to the filtered list."""
    html = client.get("/systems").get_data(as_text=True)
    assert 'id="standalone-systems"' not in html
    assert '<a href="/systems?placement=standalone#all-systems">2 standalone</a>' in html
    assert ("get_systems", DB, "none", 1, 0) in fake.calls


def test_page_past_the_end_shows_last_page(client, fake):
    fake.sectors = [_sector(i) for i in range(60)]
    html = client.get("/sectors?sectors_page=99").get_data(as_text=True)
    assert "Showing 51&ndash;60 of 60" in html


# --- Account menu, one login check per request -------------------------------------

def test_login_link_and_no_auth_lookup_without_cookie(client, fake):
    html = client.get("/").get_data(as_text=True)
    assert '<a href="/login"' in html
    assert ">Admin</a>" not in html
    assert not [call for call in fake.calls if call[0] == "auth_me"]


def test_admin_menu_with_one_lookup_per_request(client, fake):
    fake.admin = {"username": "admin", "must_change_credentials": False}
    client.set_cookie(SESSION_COOKIE_NAME, "token-value")
    html = client.get("/").get_data(as_text=True)
    assert ">Admin</a>" in html and ">Account</a>" in html and ">Logout</a>" in html
    assert ">Stats</a>" not in html  # UX.72: it is under Admin
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
        assert page_url("sector", sector_id=5) == "/sector/5"
        assert page_url("phenomenon", phenomenon_type="nebula", phenomenon_id=3) == "/phenomenon/nebula/3"
        assert page_url("galaxy", _anchor="map") == "/galaxy#map"
        with pytest.raises(BuildError):
            page_url("no_such_page")


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


# --- Old CGI URLs ----------------------------------------------------------------------

@pytest.mark.parametrize("script", ["index.py", "browse.py"])
def test_old_home_urls_redirect_keeping_page_numbers(client, script):
    result = client.get(f"/{script}?db=x&sectors_page=2&standalone_page=4")
    assert result.status_code == 301
    assert result.headers["Location"] == "/?sectors_page=2&standalone_page=4"


def test_old_url_without_params_goes_home(client):
    result = client.get("/browse.py")
    assert result.status_code == 301
    assert result.headers["Location"] == "/"


@pytest.mark.parametrize("script,location", [
    ("sectors.py", None), ("wsgi.py", None), ("nope.py", None),
    ("galaxy3d.py", "/galaxy"), ("galaxy_view.py", "/galaxy"),
])
def test_old_url_names(client, script, location):
    result = client.get(f"/{script}")
    if location is None:
        assert result.status_code == 404
        assert "text/html" in result.headers["Content-Type"]
    else:
        assert result.status_code == 301
        assert result.headers["Location"] == location


def test_every_old_page_points_at_a_route(app):
    for name, endpoint in OLD_PAGES.items():
        assert endpoint in app.view_functions, name


def test_old_url_post_is_not_replayed(client):
    """An old form POSTed to a `.py` URL is refused like any stale form."""
    assert client.post("/login.py", data={"username": "a", "password": "b"}).status_code in (400, 405)


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
        store.save_sector(SpaceSector(f"Web <i>Sector</i> {i:03d}", edge_ly=10.0), config=mysql_config)

    home = db_client.get("/")
    assert home.status_code == 200
    assert "<strong>55</strong> sectors" in home.get_data(as_text=True)
    page1 = db_client.get("/sectors")
    html = page1.get_data(as_text=True)
    assert page1.status_code == 200
    assert "Web &lt;i&gt;Sector&lt;/i&gt; 000" in html
    assert "Web &lt;i&gt;Sector&lt;/i&gt; 050" not in html
    assert f"db={mysql_config.database}" not in html  # every page is on Flask now

    page2 = db_client.get("/sectors?sectors_page=2").get_data(as_text=True)
    assert "Web &lt;i&gt;Sector&lt;/i&gt; 054" in page2
    assert "Showing 51&ndash;55 of 55" in page2


def test_real_login_session_shows_admin_menu(db_client, mysql_config):
    _username, first_password = adminAuth.bootstrap_control_schema(mysql_config)
    store.get_connection(mysql_config).close()
    resp = db_client.post("/api/auth/login", json={"username": "admin", "password": first_password})
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
    from planetgen.api.limiter import limiter
    store.get_connection(mysql_config).close()  # lay down the (empty) content schema
    limiter.reset()
    try:
        statuses = {db_client.get("/").status_code for _ in range(60)}
        assert statuses == {200}
        direct = [db_client.get("/api/sectors").status_code for _ in range(51)]
        assert direct[-1] == 429
    finally:
        limiter.reset()


def test_systems_page_lists_every_system_with_sector_and_octant(client, fake):
    in_sector = dict(_system(1, name="Inner"), sector_id=7, sector_name="Home <Sector>", quadrant="III")
    fake.systems = [in_sector] + [_system(i) for i in range(2, 60)]
    html = client.get("/systems").get_data(as_text=True)
    assert 'id="all-systems"' in html and 'id="standalone-systems"' not in html
    assert '<a href="/sector/7">Home &lt;Sector&gt;</a>' in html and "<td>III</td>" in html
    assert "Standalone" in html
    assert ("get_systems", DB, None, 50, 0) in fake.calls
    client.get("/systems?systems_page=2")
    assert ("get_systems", DB, None, 50, 50) in fake.calls


# --- The Sectors and Systems data tables (UX.41) ---------------------------------

def test_sectors_table_asks_for_the_sort_and_quadrant_filter(client, fake):
    html = client.get("/sectors?sectors_sort=systems&sectors_order=desc&quadrant=I&quadrant=unplaced").get_data(as_text=True)
    assert ("sectors", "systems", True, ["I", "unplaced"]) in fake.asked
    assert '<th scope="col" data-col="systems" aria-sort="descending">' in html
    # The default sort is by distance from the core.
    assert '<th scope="col" data-col="distance" aria-sort="none">' in html
    assert 'href="/sectors?quadrant=I&amp;quadrant=unplaced#sectors"' in html
    assert 'name="quadrant" value="I" checked' in html
    assert '<span class="datatable-option-label">Quadrant I</span>' in html


def test_sectors_table_default_sort_is_distance(client, fake):
    html = client.get("/sectors").get_data(as_text=True)
    assert ("sectors", "distance", False, []) in fake.asked
    assert '<th scope="col" data-col="distance" aria-sort="ascending">' in html


def test_systems_table_takes_its_sort_and_filters(client, fake):
    html = client.get("/systems?systems_sort=sector&binary=yes&placement=sector&octant=I").get_data(as_text=True)
    assert ("systems", "sector", False, True, "sector", ["I"]) in fake.asked
    assert 'name="systems_sort" value="sector"' in html


def test_both_binary_choices_filter_nothing(client, fake):
    client.get("/systems?binary=yes&binary=no")
    assert ("systems", "name", False, None, None, []) in fake.asked


def test_unrelated_address_parameters_are_not_carried_into_links(client, fake):
    html = client.get("/sectors?sectors_sort=name&stray=1&db=someone_elses_db").get_data(as_text=True)
    assert "stray" not in html and "someone_elses_db" not in html


def test_table_route_serves_the_sectors_and_systems_tables(client, fake):
    fake.sectors = [_sector(1, placed=True), _sector(2, placed=False)]
    body = client.get("/table/sectors?sort=systems&order=desc&facets=1").get_json()
    assert body["total"] == 2
    assert body["rows"][0][0] == {"text": "Sector 001", "href": "/sector/1"}
    assert body["rows"][0][3] == {"text": "Quadrant I", "href": "/galaxy?quadrant=I"}
    assert body["rows"][1][3] == {"text": "Unplaced", "muted": True}
    assert body["facets"]["quadrant"] == [{"value": "I", "label": "Quadrant I", "count": 2}]
    assert ("sectors", "systems", True, []) in fake.asked

    body = client.get("/table/systems?facets=1&placement=standalone").get_json()
    assert [cell["text"] for cell in body["rows"][0]] == ["System 000", "Standalone", "–", "No", "G2V"]
    assert {o["value"]: o["label"] for o in body["facets"]["binary"]} == {"yes": "Binary", "no": "Single star"}
    assert client.get("/table/standalone-systems").status_code == 404


# --- Uncharted stars on the Systems page (UX.87) ------------------------------------

def _uncharted_star(i):
    return {"id": 50 + i, "star_type": "B2V", "yerkes_class": "V", "population": "young", "mass_solar": 9.0,
            "radius_solar": 5.0, "temperature_k": 22000.0, "luminosity_sol": 6000.0 - i, "age_gy": 0.02,
            "ring_index": 3, "layer_index": 0, "ring_slot_index": 1, "designation": "6000100001",
            "x": 12.0, "y": 3.0, "z": 0.5, "local_x": 0.4, "local_y": -0.2, "local_z": 0.1, "galactic_radius_pc": 12.4}


def test_the_systems_page_offers_the_uncharted_stars_off_by_default(client, fake):
    fake.uncharted = [_uncharted_star(i) for i in range(3)]
    html = client.get("/systems").get_data(as_text=True)
    assert 'id="all-systems"' in html and "3 uncharted stars" in html
    assert "/systems?uncharted=1" in html
    assert "Uncharted B2V star" not in html
    assert not [call for call in fake.calls if call[0] == "get_uncharted_systems" and call[2] != 1]


def test_the_uncharted_filter_lists_every_scattered_star_with_where_it_is(client, fake):
    fake.uncharted = [_uncharted_star(i) for i in range(3)]
    html = client.get("/systems?uncharted=1").get_data(as_text=True)
    assert ("uncharted", "luminosity", True) in fake.asked, "brightest first"
    assert "Uncharted B2V star 50" in html and "6,000" in html and "B2V (Young disk)" in html
    assert "6000100001 (ring 3, layer 0, slot 1)" in html
    assert "/galaxy?sector=6000100001&amp;open=1" in html
    assert "12.00, 3.00, 0.50 pc" in html and "0.40, -0.20, 0.10 pc" in html
    assert "generating the whole sector" in html or "recommended" in html
    assert 'value="Generate"' not in html and ">Generate<" not in html.replace("Generate</th>", ""), "no button for visitors"


def test_an_admin_gets_a_generate_button_per_uncharted_star(client, fake, monkeypatch):
    fake.uncharted = [_uncharted_star(0)]
    fake.admin = {"username": "boss", "must_change_credentials": False}
    client.set_cookie(SESSION_COOKIE_NAME, "token-value")
    html = client.get("/systems?uncharted=1").get_data(as_text=True)
    assert 'action="/systems/uncharted/50/generate"' in html and "Generate" in html


def _csrf_form(app, client):
    nonce = "n" * 43
    client.set_cookie(csrf.COOKIE_NAME, nonce)
    session = client.get_cookie(SESSION_COOKIE_NAME)
    with app.app_context():
        return {csrf.FIELD_NAME: csrf._sign(nonce, session.value if session else "")}


def test_generating_an_uncharted_star_opens_its_system(app, client, fake, monkeypatch):
    fake.admin = {"username": "boss", "must_change_credentials": False}
    client.set_cookie(SESSION_COOKIE_NAME, "token-value")
    calls = []

    def edit(cookie, db, method, path, body=None):
        calls.append((method, path))
        return {"id": 77}

    monkeypatch.setattr(apiclient, "admin_edit", edit)
    response = client.post("/systems/uncharted/50/generate", data=_csrf_form(app, client))
    assert response.status_code == 303 and response.headers["Location"].endswith("/system/77")
    assert calls == [("POST", "/uncharted-systems/50/generate")]


def test_a_visitor_cannot_generate_an_uncharted_star(app, client, fake, monkeypatch):
    monkeypatch.setattr(apiclient, "admin_edit", lambda *args, **kwargs: pytest.fail("must not call"))
    response = client.post("/systems/uncharted/50/generate", data=_csrf_form(app, client))
    assert response.status_code == 303 and "uncharted=1" in response.headers["Location"]


def test_a_star_that_is_no_longer_waiting_flashes_a_message(app, client, fake, monkeypatch):
    fake.admin = {"username": "boss", "must_change_credentials": False}
    client.set_cookie(SESSION_COOKIE_NAME, "token-value")

    def gone(*args, **kwargs):
        raise apiclient.NotFoundError("gone")

    monkeypatch.setattr(apiclient, "admin_edit", gone)
    response = client.post("/systems/uncharted/50/generate", data=_csrf_form(app, client), follow_redirects=True)
    assert "no longer waiting" in response.get_data(as_text=True)


def test_an_api_refusal_is_shown_on_the_uncharted_list(app, client, fake, monkeypatch):
    fake.admin = {"username": "boss", "must_change_credentials": False}
    client.set_cookie(SESSION_COOKIE_NAME, "token-value")

    def refuse(*args, **kwargs):
        raise apiclient.ApiError("planetGen API error (409): That star already has a system.", 409)

    monkeypatch.setattr(apiclient, "admin_edit", refuse)
    response = client.post("/systems/uncharted/50/generate", data=_csrf_form(app, client), follow_redirects=True)
    assert "That star already has a system." in response.get_data(as_text=True)
