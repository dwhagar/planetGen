# tests/test_web_system_phen.py

"""
The Flask-served system and phenomenon pages (`web/system_pages.py`):
`/system/<id>`, `/phenomena` and `/phenomenon/<type>/<id>`, plus the old
CGI URLs that redirect there.

Most tests fake the data layer (the `apiclient` functions the views
call), like `test_web_pages.py`; the ones at the bottom run against a
real throwaway database through the in-process transport and are skipped
without a MySQL test server.
"""

import re

import pytest

from api.app import create_app
from api.authz import SESSION_COOKIE_NAME
from api.config import Config

import web  # noqa: F401 -- puts src/html/lib on sys.path
import apiclient  # noqa: E402
from fmt import linkify_location, nearest_neighbors_location  # noqa: E402
from stellarObjects import _db, adminAuth  # noqa: E402
from stellarObjects.config import SystemConfig  # noqa: E402
from stellarObjects.spaceSector import SpaceSector  # noqa: E402
from stellarObjects.systemData import StarSystem  # noqa: E402
from web import csrf  # noqa: E402

DB = "planetgen_web_test"


class _FakeConfig(Config):
    WEB_DATABASE = DB
    SESSION_COOKIE_SECURE = False
    SECRET_KEY = "test-secret"


def _star(i, name="Sol", role="single"):
    return {
        "id": i, "role": role, "name": name, "star_type": "G2V", "mass_kg": 1.989e30,
        "radius_km": 696_000.0, "temperature_k": 5778.0, "luminosity_w": 3.828e26,
    }


def _system_detail(system_id=5, name="Kepler <b>42</b>", sector_id=7):
    return {
        "id": system_id, "name": name, "sector_id": sector_id, "quadrant": "III",
        "location": "Voranthis & Kelmoor -- nearest: Alpha (4.2 ly), Beta (5.1 ly)" if sector_id else None,
        "is_binary": 0, "binary_type": None, "binary_configuration": None,
        "binary_mutual_position_x_km": None, "binary_mutual_position_y_km": None,
        "binary_mutual_position_z_km": None,
        "wikijs_url": None, "mediawiki_url": None,
        "stars": [_star(1)], "planets": [], "belts": [], "comets": [],
        "sector_siblings": [{"id": 8, "name": "Alpha"}, {"id": 9, "name": "Beta"}],
        "nearest_neighbors": [{"id": 8, "name": "Alpha", "distance_ly": 4.2},
                              {"id": 9, "name": "Be<ta>", "distance_ly": 5.1}],
    }


def _phenomenon_row(i, kind="nebula", sector_id=None):
    return {
        "id": i, "type": kind, "name": f"Phenomenon {i:03d}", "descriptor": "emission_nebula",
        "radius_ly": 12.5 if kind == "nebula" else None, "sector_id": sector_id,
        "sector_name": "Home <Sector>" if sector_id else None, "placed": i % 2,
    }


class FakeData:
    """Stands in for the `apiclient` functions these pages call."""

    def __init__(self):
        self.system = _system_detail()
        self.sections = {"overview": "An *overview*.", "stars": {"1": "A yellow star."},
                         "planets": {}, "moons": {}, "belts": {}, "comets": {}}
        self.phenomena = [_phenomenon_row(i, sector_id=3 if i % 3 == 0 else None) for i in range(3)]
        self.phenomenon = {
            "id": 4, "name": "Crab <Nebula>", "nebula_type": "supernova_remnant", "radius_ly": 5.5,
            "composition": "Hydrogen & helium", "formation_cause": "stellar death",
            "galactic_orbital_speed_kms": 220.0, "galactic_orbital_period_gy": 0.25,
            "galactic_radius_pc": 1000.0, "sector_id": 3, "sector_name": "Crab <Sector>",
        }
        self.wiki_config = {"wikijs": True, "mediawiki": True}
        self.upload_error = None
        self.admin = None
        self.calls = []

    def get_system(self, db, system_id):
        self.calls.append(("get_system", db, system_id))
        if int(system_id) != self.system["id"]:
            raise apiclient.NotFoundError(f"no such system: {system_id}")
        return self.system

    def get_system_sections(self, db, system_id):
        self.calls.append(("get_system_sections", db, system_id))
        return self.sections

    def get_system_facilities(self, db, system_id):
        self.calls.append(("get_system_facilities", db, system_id))
        return []

    def get_system_text(self, db, system_id, fmt):
        self.calls.append(("get_system_text", db, system_id, fmt))
        return {"id": system_id, "format": fmt, "content": f"== {fmt} <page> ==\nline"}

    def get_wiki_config(self):
        self.calls.append(("get_wiki_config",))
        return self.wiki_config

    def upload_system_to_wiki(self, cookie_header, db, system_id, backend, path=None):
        self.calls.append(("upload_system_to_wiki", cookie_header, db, system_id, backend, path))
        if self.upload_error:
            raise self.upload_error
        return {"id": 1, "path": path, "title": "x", "url": "https://wiki.example/x"}

    def get_phenomena(self, db, limit=None, offset=None):
        self.calls.append(("get_phenomena", db, limit, offset))
        return {"items": self.phenomena[offset:offset + limit], "total": len(self.phenomena),
                "limit": limit, "offset": offset}

    def get_phenomenon(self, db, phenomenon_type, phenomenon_id):
        self.calls.append(("get_phenomenon", db, phenomenon_type, phenomenon_id))
        return self.phenomenon

    def auth_me(self, cookie_header):
        self.calls.append(("auth_me", cookie_header))
        return self.admin


_FAKED = ("get_system", "get_system_sections", "get_system_facilities", "get_system_text", "get_wiki_config", "upload_system_to_wiki",
          "get_phenomena", "get_phenomenon", "auth_me")


@pytest.fixture
def fake(monkeypatch):
    data = FakeData()
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
    return app.test_client()


def _as_admin(client, fake):
    fake.admin = {"username": "admin", "must_change_credentials": False}
    client.set_cookie(SESSION_COOKIE_NAME, "token-value")


def _csrf(app, client):
    nonce = "n" * 43
    client.set_cookie(csrf.COOKIE_NAME, nonce)
    session = client.get_cookie(SESSION_COOKIE_NAME)
    with app.app_context():
        return csrf._sign(nonce, session.value if session else "")  # bound to the login session


def _section_current(html, label):
    nav = re.search(r'<nav class="site-sections" aria-label="Main">.*?</nav>', html, re.S).group(0)
    return re.search(rf'<a href="[^"]*" aria-current="page">{label}</a>', nav) is not None


# --- /system/<id> -------------------------------------------------------------------

def _url(app, name, **params):
    """`page_url` as it appears in the HTML (escaped). Used for links to
    pages other groups are moving, whose URLs change when they land."""
    from markupsafe import escape
    from web.helpers import page_url
    with app.test_request_context("/"):
        return str(escape(page_url(name, **params)))


def _nav_links(app, kind, entity_id):
    """The expected Navigate links (`/nav?from=`/`?to=<kind>:<id>`),
    HTML-escaped."""
    from markupsafe import escape
    from web.system_pages import nav_links
    with app.test_request_context("/"):
        return {key: str(escape(url)) for key, url in nav_links(kind, entity_id).items()}


def test_system_page_renders(app, client, fake):
    resp = client.get("/system/5")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert '<h1 class="page-title">Kepler &lt;b&gt;42&lt;/b&gt;</h1>' in html
    assert "<title>Kepler &lt;b&gt;42&lt;/b&gt; - " in html
    # Breadcrumbs: Home > Sectors > <sector, from the location prefix> > system.
    crumbs = re.search(r'<nav class="breadcrumbs" aria-label="Breadcrumb">(.*?)</nav>', html, re.S).group(1)
    assert '<a href="/sectors">Sectors</a>' in crumbs
    assert f'<a href="{_url(app, "sector", sector_id=7)}">Voranthis &amp; Kelmoor</a>' in crumbs
    assert '<span aria-current="page">Kepler &lt;b&gt;42&lt;/b&gt;</span>' in crumbs
    assert _section_current(html, "Sectors")
    # Badges, nav links, location links: all plain GET links.
    assert '<span class="badge">Octant III</span>' in html and "Single star" in html
    links = _nav_links(app, "system", 5)
    assert links["from"].startswith("/nav?from=system") and links["to"].startswith("/nav?to=system")
    assert f'href="{links["from"]}">Navigate from here</a>' in html
    assert f'href="{links["to"]}">Navigate to here</a>' in html
    assert '<a href="/system/8">Alpha</a> (4.2 ly)' in html
    assert '<a href="/system/9">Be&lt;ta&gt;</a> (5.1 ly)' in html
    # The map, the body list, the stars table, and the page's scripts.
    assert 'id="sysmap-' in html
    assert 'class="system-list system-list-root"' in html
    assert "<p>An *overview*.</p>" in html and "A yellow star." in html
    assert "<h2>Stars</h2>" in html
    assert re.search(r'<script type="module" src="/static/systemmap.js\?v=[^"]+"></script>', html)
    assert re.search(r'<script type="module" src="/static/copycode.js\?v=[^"]+"></script>', html)
    # No POST forms, no db in any moved-page URL, no code box by default.
    assert "<form method=\"post\"" not in html
    assert 'id="system-code"' not in html
    assert not re.search(r'href="/system/[^"]*db=', html)
    assert {call[1] for call in fake.calls if call[0].startswith("get_system")} == {DB}
    # Logged out: no wiki lookup and no upload form.
    assert not [c for c in fake.calls if c[0] in ("get_wiki_config", "auth_me")]
    assert "Upload to Wiki" not in html


def test_standalone_system_breadcrumbs_and_no_nav(client, fake):
    fake.system = _system_detail(sector_id=None)
    html = client.get("/system/5").get_data(as_text=True)
    crumbs = re.search(r'<nav class="breadcrumbs" aria-label="Breadcrumb">(.*?)</nav>', html, re.S).group(1)
    assert '<a href="/systems">Systems</a>' in crumbs
    assert _section_current(html, "Systems")
    assert "Navigate from here" not in html
    assert "Location:" not in html


def test_system_page_code_views(client, fake):
    html = client.get("/system/5?code=wikitext").get_data(as_text=True)
    assert ("get_system_text", DB, 5, "wikitext") in fake.calls
    assert 'id="system-code"' in html and 'data-copy-target="system-code"' in html
    assert "== wikitext &lt;page&gt; ==" in html
    assert 'href="/system/5#system-panel"' in html and ">Hide Wikitext</a>" in html
    assert 'href="/system/5?code=markdown#system-panel"' in html
    html = client.get("/system/5?code=bogus").get_data(as_text=True)
    assert 'id="system-code"' not in html
    assert 'href="/system/5?code=wikitext#system-panel">Wikitext</a>' in html


def test_system_page_wiki_links(client, fake):
    fake.system["wikijs_url"] = "https://wiki.example/a?x=1&y=2"
    html = client.get("/system/5").get_data(as_text=True)
    assert 'href="https://wiki.example/a?x=1&amp;y=2" target="_blank" rel="noopener noreferrer">View on Wiki.js' in html
    assert "View on MediaWiki" not in html


def test_unknown_system_is_404(client, fake):
    resp = client.get("/system/6")
    assert resp.status_code == 404
    assert resp.mimetype == "text/html"
    assert client.get("/system/abc").status_code == 404


def test_admin_sees_upload_form_with_remaining_backends(client, fake):
    _as_admin(client, fake)
    fake.system["mediawiki_url"] = "https://wiki.example/m"
    html = client.get("/system/5").get_data(as_text=True)
    form = re.search(r'<form method="post" action="/system/5".*?</form>', html, re.S).group(0)
    assert 'name="csrf_token"' in form
    assert 'value="wikijs" checked' in form
    assert 'value="mediawiki"' not in form
    assert 'placeholder="e.g. systems/Kepler &lt;b&gt;42&lt;/b&gt;"' in form


def test_admin_upload_form_hidden_when_nothing_left(client, fake):
    _as_admin(client, fake)
    fake.wiki_config = {"wikijs": False, "mediawiki": False}
    html = client.get("/system/5").get_data(as_text=True)
    assert 'id="wiki-upload"' not in html


def test_upload_post_redirects_to_get(app, client, fake):
    _as_admin(client, fake)
    token = _csrf(app, client)
    resp = client.post("/system/5", data={csrf.FIELD_NAME: token, "backend": "wikijs", "path": " systems/k42 "})
    assert resp.status_code == 303
    assert resp.headers["Location"] == "/system/5?wiki=uploaded#wiki-upload"
    upload = [c for c in fake.calls if c[0] == "upload_system_to_wiki"][0]
    assert upload[2:] == (DB, 5, "wikijs", "systems/k42")
    assert f"{SESSION_COOKIE_NAME}=token-value" in upload[1]
    html = client.get("/system/5?wiki=uploaded").get_data(as_text=True)
    assert "Uploaded to the wiki." in html


@pytest.mark.parametrize("status,code,message", [
    (409, "exists", "A page already exists at that location."),
    (400, "invalid", "The upload was rejected."),
    (501, "unconfigured", "That wiki is not configured"),
    (502, "failed", "The upload failed."),
])
def test_upload_errors_become_fixed_messages(app, client, fake, status, code, message):
    _as_admin(client, fake)
    fake.upload_error = apiclient.ApiError("detail <x>", status_code=status)
    resp = client.post("/system/5", data={csrf.FIELD_NAME: _csrf(app, client), "backend": "mediawiki"})
    assert resp.headers["Location"] == f"/system/5?wiki={code}#wiki-upload"
    html = client.get(f"/system/5?wiki={code}").get_data(as_text=True)
    assert message in html and "detail &lt;x&gt;" not in html


def test_upload_rejects_unknown_backend_without_calling_api(app, client, fake):
    _as_admin(client, fake)
    resp = client.post("/system/5", data={csrf.FIELD_NAME: _csrf(app, client), "backend": "evil"})
    assert resp.headers["Location"] == "/system/5?wiki=invalid#wiki-upload"
    assert not [c for c in fake.calls if c[0] == "upload_system_to_wiki"]


def test_upload_needs_admin_and_csrf(app, client, fake):
    assert client.post("/system/5", data={"backend": "wikijs"}).status_code == 400
    resp = client.post("/system/5", data={csrf.FIELD_NAME: _csrf(app, client), "backend": "wikijs"})
    assert resp.status_code == 403
    assert not [c for c in fake.calls if c[0] == "upload_system_to_wiki"]
    # A visitor can't show the outcome message either.
    assert "Uploaded to the wiki." not in client.get("/system/5?wiki=uploaded").get_data(as_text=True)


# --- /phenomena and /phenomenon/<type>/<id> ------------------------------------------

def test_phenomena_list(app, client, fake):
    resp = client.get("/phenomena")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert _section_current(html, "Phenomena")
    assert "3 phenomena" in html
    assert '<a href="/phenomenon/nebula/0">Phenomenon 000</a>' in html
    assert f'<a href="{_url(app, "sector", sector_id=3)}">Home &lt;Sector&gt;</a>' in html
    assert "Emission nebula" in html and "3.83 pc (12.5 ly)" in html
    assert '<form method="post"' not in html
    assert ("get_phenomena", DB, 50, 0) in fake.calls


def test_phenomena_pages_with_get_links(client, fake):
    fake.phenomena = [_phenomenon_row(i) for i in range(120)]
    html = client.get("/phenomena?page=2").get_data(as_text=True)
    assert ("get_phenomena", DB, 50, 50) in fake.calls
    assert "Phenomenon 050" in html and "Phenomenon 049" not in html
    assert 'href="/phenomena?page=3#phenomena-list"' in html
    assert "Showing 51&ndash;100 of 120" in html


def test_phenomenon_detail(app, client, fake):
    resp = client.get("/phenomenon/nebula/4")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert ("get_phenomenon", DB, "nebula", 4) in fake.calls
    assert '<h1 class="page-title">Crab &lt;Nebula&gt;</h1>' in html
    assert '<a href="/phenomena">Phenomena</a>' in html
    assert _section_current(html, "Phenomena")
    assert '<span class="badge">Nebula</span>' in html
    assert "1 kpc (3,262 ly) from Galactic Center" in html
    assert f'<a href="{_url(app, "sector", sector_id=3)}">Crab &lt;Sector&gt;</a>' in html
    links = _nav_links(app, "nebula", 4)
    assert f'href="{links["from"]}">Navigate from here</a>' in html
    assert f'href="{links["to"]}">Navigate to here</a>' in html
    assert "4" in links["from"] and "nebula" in links["from"]
    assert '<th scope="row">Composition</th><td>Hydrogen &amp; helium</td>' in html
    assert "<td>Supernova remnant</td>" in html
    assert 'id="phenomenonmap-svg"' in html
    # The diagram's scripts come from the template with absolute URLs.
    assert re.search(r'<script src="/static/mapzoom.js\?v=[^"]+" defer></script>', html)
    assert re.search(r'<script type="module" src="/static/phenomenonmap.js\?v=[^"]+"></script>', html)
    assert 'src="static/' not in html


def test_quasar_detail(client, fake):
    fake.phenomenon = {"id": 1, "name": "Core Q", "black_hole_mass_solar": 2.5e9, "eddington_ratio": 0.42,
                       "is_radio_loud": 1, "jet_length_ly": 150000.0, "galactic_radius_pc": 0.0,
                       "sector_id": None}
    html = client.get("/phenomenon/quasar/1").get_data(as_text=True)
    assert '<span class="badge">Quasar</span>' in html
    assert "<h2 id=\"phenomenon-data-heading\">Quasar Data</h2>" in html
    assert "2.50e+09 solar masses" in html and "<td>42%</td>" in html and "150,000 ly" in html
    assert "Galactic Orbital" not in html


def test_phenomenon_unknown_type_is_404(client, fake):
    assert client.get("/phenomenon/quasar_x/4").status_code == 404
    assert not [c for c in fake.calls if c[0] == "get_phenomenon"]


def test_black_hole_skips_missing_optional_rows(client, fake):
    fake.phenomenon = {"id": 2, "name": "BH", "mass_solar": 10.0, "has_accretion_disk": 1,
                       "galactic_orbital_speed_kms": None, "sector_id": None}
    html = client.get("/phenomenon/black_hole/2").get_data(as_text=True)
    assert "10.00 solar masses" in html and "<td>Yes</td>" in html
    assert "Galactic Orbital Speed" not in html
    assert "Sector:" not in html


# --- Links elsewhere now point here --------------------------------------------------

def test_moved_pages_are_routes(app):
    for name in ("system", "phenomena", "phenomenon"):
        assert f"web.{name}" in app.view_functions


def test_header_phenomena_section_links_here(client, fake):
    html = client.get("/phenomena").get_data(as_text=True)
    assert '<a href="/phenomena" aria-current="page">Phenomena</a>' in html


def test_location_link_hook():
    url = lambda system_id: f"/system/{system_id}"  # noqa: E731
    html = nearest_neighbors_location("S -- nearest: x", [{"id": 3, "name": "A&B", "distance_ly": 1.0}], url)
    assert html == 'S -- nearest: <a href="/system/3">A&amp;B</a> (1.0 ly)'
    html = linkify_location("S -- nearest: A (1.0 ly), Z (2.0 ly)", {"A": 3}, url)
    assert html == 'S -- nearest: <a href="/system/3">A</a> (1.0 ly), Z (2.0 ly)'


# --- Old CGI URLs ---------------------------------------------------------------------

@pytest.mark.parametrize("url,location", [
    ("/system.py?db=x&id=12", "/system/12"),
    ("/system.py?db=x&id=12&code=markdown", "/system/12?code=markdown"),
    ("/system.py?db=x&id=../admin", "/systems"),
    ("/system.py", "/systems"),
    ("/phenomenon.py?db=x&type=black_hole&id=3", "/phenomenon/black_hole/3"),
    ("/phenomenon.py?db=x&type=../x&id=3", "/phenomena"),
    ("/phenomenon.py?db=x", "/phenomena"),
    ("/phenomena.py?db=x&page=3", "/phenomena?page=3"),
    ("/phenomena.py", "/phenomena"),
])
def test_old_system_and_phenomenon_urls_redirect(client, url, location):
    result = client.get(url)
    assert result.status_code == 301
    assert result.headers["Location"] == location


# --- Real database, in-process ---------------------------------------------------------

@pytest.fixture
def db_app(mysql_config, monkeypatch):
    class RealConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "test-secret"

    def no_http(*args, **kwargs):
        raise AssertionError("HTTP transport used inside a Flask request")
    monkeypatch.setattr(apiclient, "_http_transport", no_http)
    application = create_app(RealConfig)
    application.testing = True
    return application


def _save_system(mysql_config, name=None, moons=True):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = True
    cfg.MOONS = moons
    cfg.BINARY_SYSTEM = False
    if name:
        cfg.NAME = name
    # The body list's Habitable/Inhabited chips belong to planet rows, and
    # about 1 in 20 random G2V systems has no planet (nothing, or only a
    # belt), so draw until one has a planet.
    system = StarSystem(system_config=cfg)
    for _ in range(50):
        if system.planet_count:
            break
        system = StarSystem(system_config=cfg)
    assert system.planet_count
    return _db.save_system(system, cfg, config=mysql_config)


def test_real_system_page_lists_bodies_and_code(db_app, mysql_config):
    client = db_app.test_client()
    system_id = _save_system(mysql_config)
    resp = client.get(f"/system/{system_id}")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert 'class="system-list system-list-root"' in html
    # One ordered list (TODO 48): no Planets & Moons table any more, and
    # the stars table comes first, right under the map.
    assert "Planets &amp; Moons" not in html and "<h2>Asteroid Belts" not in html
    assert '<span class="stat">Class ' in html
    assert html.index("<h2>Stars</h2>") < html.index('id="system-panel"')
    assert 'id="system-code"' not in html
    for fmt, marker in (("wikitext", "[[Category:Star Systems]]"), ("markdown", "| Property | Value |")):
        page = client.get(f"/system/{system_id}?code={fmt}").get_data(as_text=True)
        assert 'id="system-code"' in page
        assert marker in page.replace("&#39;", "'")


def test_real_system_in_sector_links_neighbours(db_app, mysql_config):
    sector = SpaceSector("Web Sector", edge_ly=10.0)
    for position in ((1.0, 1.0, 1.0), (-2.0, 0.5, 3.0)):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "M5V"
        cfg.PLANETS = False
        cfg.BINARY_SYSTEM = False
        sector.add_system(StarSystem(system_config=cfg), position=position, system_config=cfg)
    sector_id = _db.save_sector(sector, config=mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        ids = [r["id"] for r in conn.execute(
            "SELECT id FROM star_systems WHERE sector_id = ? ORDER BY id", (sector_id,)).fetchall()]
    finally:
        conn.close()
    html = db_app.test_client().get(f"/system/{ids[0]}").get_data(as_text=True)
    assert f'<a href="/system/{ids[1]}">' in html
    assert f'<a href="{_url(db_app, "sector", sector_id=sector_id)}">Web Sector</a>' in html
    assert "Navigate from here" in html


def test_real_system_name_is_escaped(db_app, mysql_config):
    system_id = _save_system(mysql_config, name="<script>alert(1)</script>", moons=False)
    html = db_app.test_client().get(f"/system/{system_id}").get_data(as_text=True)
    assert "<script>alert(1)</script>" not in html
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in html


def test_real_unknown_system_is_404(db_app, mysql_config):
    _db.get_connection(mysql_config).close()
    resp = db_app.test_client().get("/system/999999999")
    assert resp.status_code == 404
    assert "Traceback" not in resp.get_data(as_text=True)


def test_real_phenomena_list_and_detail(db_app, mysql_config):
    from stellarObjects.roguePlanetData import InterstellarComet, RoguePlanet

    sector = SpaceSector(name="Phenomena Sector", edge_ly=40.0)
    sector.add_phenomenon(RoguePlanet(SystemConfig()), "rogue-planet")
    sector.add_phenomenon(InterstellarComet(SystemConfig()), "comet")
    _db.save_sector(sector, config=mysql_config)
    client = db_app.test_client()

    html = client.get("/phenomena").get_data(as_text=True)
    assert "2 phenomena" in html
    links = re.findall(r'<a href="(/phenomenon/[a-z_]+/\d+)">', html)
    assert len(links) == 2
    for link in links:
        resp = client.get(link)
        page = resp.get_data(as_text=True)
        assert resp.status_code == 200
        # Rogue planets and comets get a rendered view, built from the
        # real row's fields.
        assert 'id="phenomrender"' in page
        assert "Galactic Orbital Speed" in page


def test_real_admin_upload_without_wiki_config(db_app, mysql_config):
    """With no wiki configured, an admin sees no form, and a hand-made
    POST comes back as a fixed message: first "change your default
    credentials" (the API's own rule), then "not configured"."""
    _username, first_password = adminAuth.bootstrap_control_schema(mysql_config)
    system_id = _save_system(mysql_config, moons=False)
    client = db_app.test_client()
    assert client.post("/api/auth/login", json={"username": "admin", "password": first_password}).status_code == 200
    html = client.get(f"/system/{system_id}").get_data(as_text=True)
    assert 'id="wiki-upload"' not in html
    token = _csrf(db_app, client)
    resp = client.post(f"/system/{system_id}", data={csrf.FIELD_NAME: token, "backend": "mediawiki"})
    assert resp.headers["Location"] == f"/system/{system_id}?wiki=forbidden#wiki-upload"
    assert "Change the admin username and password the installer set first." in client.get(resp.headers["Location"]).get_data(as_text=True)

    assert client.post("/api/auth/change-credentials", json={
        "current_password": first_password, "new_username": "boss", "new_password": "a-long-new-password-1",
    }).status_code == 200
    token = _csrf(db_app, client)  # the change re-issued the session; tokens are bound to it
    resp = client.post(f"/system/{system_id}", data={csrf.FIELD_NAME: token, "backend": "mediawiki"})
    assert resp.status_code == 303
    assert resp.headers["Location"] == f"/system/{system_id}?wiki=unconfigured#wiki-upload"
    html = client.get(resp.headers["Location"]).get_data(as_text=True)
    assert "That wiki is not configured for this site." in html
