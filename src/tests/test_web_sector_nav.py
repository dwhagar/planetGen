# tests/test_web_sector_nav.py

"""
The Flask-served sector page (`/sector/<id>`, `web/sector_page.py`) and
NAV page (`/nav`, `web/nav_page.py`), plus the `sector.py`/`nav.py` CGI
shims that now redirect to them.

Most tests fake the data layer by monkeypatching the `apiclient`
functions the views call (see `test_web_pages.py`). The tests at the
bottom use a real throwaway database through the in-process transport and
are skipped without a MySQL test server.
"""

import html as html_lib
import json
import re

import pytest

from api.app import create_app
from api.authz import SESSION_COOKIE_NAME
from api.config import Config

import web  # noqa: F401 -- puts src/html/lib on sys.path
import apiclient  # noqa: E402
from stellarObjects import _db, adminAuth  # noqa: E402
from stellarObjects.config import SystemConfig  # noqa: E402
from stellarObjects.spaceSector import SpaceSector  # noqa: E402
from stellarObjects.systemData import StarSystem  # noqa: E402
from web import csrf  # noqa: E402
from web.helpers import page_url  # noqa: E402
from web.nav_page import endpoint, nav_url  # noqa: E402

from tests.webpage_support import run_page  # noqa: E402

DB = "planetgen_web_test"


class _FakeConfig(Config):
    WEB_DATABASE = DB
    SESSION_COOKIE_SECURE = False
    SECRET_KEY = "test-secret"


def _star(star_type="G2V"):
    return {"star_type": star_type, "temperature_k": 5772.0, "radius_km": 696000.0, "luminosity_w": 3.8e26}


def _sector_system(system_id, name, distance, x=1.0):
    return {
        "id": system_id, "name": name, "is_binary": 0, "binary_type": None, "stars": [_star()],
        "quadrant": "+x+y+z", "location": f"Fake Sector -- nearest: Other ({distance:.1f} ly)",
        "center_distance_ly": distance,
        "position_x_mpc": x, "position_y_mpc": 0.5, "position_z_mpc": -0.5,
    }


def _sector_detail(sector_id=5, name="Fake Sector", systems=None, placed=True, wiki_url=None):
    return {
        "id": sector_id, "name": name, "edge_mpc": 3066.0, "edge_ly": 10.0,
        "shell_index": 5 if placed else None, "shell_slot_index": 100 if placed else None,
        "placed": placed,
        "center_x_pc": 500.0 if placed else None, "center_y_pc": 200.0 if placed else None,
        "center_z_pc": 10.0 if placed else None,
        "wiki_url": wiki_url,
        "systems": systems if systems is not None else [
            _sector_system(1001, "Alpha", 2.0), _sector_system(1002, "Other", 1.0, x=-1.0),
        ],
        "phenomena": [{
            "id": 3, "type": "nebula", "name": "Veil", "descriptor": "emission", "radius_ly": 2.5,
            "distance_ly": 4.0, "offset_x_ly": 1.0, "offset_y_ly": 0.0, "offset_z_ly": 0.0,
        }],
        "neighbors": [{
            "direction_pc": (1.0, 0.0, 0.0), "shell_index": 5, "shell_slot_index": 101,
            "designation": "ABC", "exists": True, "sector_id": 6, "sector_name": "Next Door",
        }],
    }


class FakeData:
    """Stands in for the `apiclient` functions these pages call."""

    def __init__(self):
        self.sectors = {5: _sector_detail(), 9: _sector_detail(9, "Far Sector", systems=[
            _sector_system(2001, "Distant", 3.0)])}
        self.systems = {
            1001: {"id": 1001, "name": "Alpha", "sector_id": 5},
            1002: {"id": 1002, "name": "Other", "sector_id": 5},
            1500: {"id": 1500, "name": "Waypoint", "sector_id": 5},
            2001: {"id": 2001, "name": "Distant", "sector_id": 9},
            3000: {"id": 3000, "name": "Loner", "sector_id": None},
        }
        self.phenomena = {("nebula", 3): {"id": 3, "name": "Veil", "galactic_radius_pc": 540.0}}
        self.admin = None
        self.wiki_config = {"wikijs": True, "mediawiki": False}
        self.calls = []
        self.nav_error = None
        self.action_error = None

    def get_sector(self, db, sector_id):
        self.calls.append(("get_sector", db, sector_id))
        try:
            return self.sectors[int(sector_id)]
        except (KeyError, ValueError):
            raise apiclient.NotFoundError(f"No such sector: {sector_id}")

    def get_sectors(self, db, limit=None, offset=None):
        self.calls.append(("get_sectors", db, limit))
        items = [{"id": key, "name": value["name"]} for key, value in sorted(self.sectors.items())]
        return {"items": items, "total": len(items), "limit": limit, "offset": 0}

    def get_system(self, db, system_id):
        self.calls.append(("get_system", db, system_id))
        try:
            return self.systems[int(system_id)]
        except KeyError:
            raise apiclient.NotFoundError(f"No such system: {system_id}")

    def get_phenomenon(self, db, phenomenon_type, phenomenon_id):
        self.calls.append(("get_phenomenon", db, phenomenon_type, phenomenon_id))
        try:
            return self.phenomena[(phenomenon_type, int(phenomenon_id))]
        except KeyError:
            raise apiclient.NotFoundError("No such phenomenon")

    def get_nav(self, db, from_id, to_id, from_kind="system", to_kind="system", from_type=None, to_type=None):
        self.calls.append(("get_nav", db, from_id, to_id, from_kind, to_kind, from_type, to_type))
        if self.nav_error:
            raise self.nav_error
        return {
            "scope": "sector",
            "direct": {"distance_ly": 3.25, "azimuth_deg": 45.0, "altitude_deg": -2.5},
            "warp_times": [{"warp_factor": 1, "velocity_multiple_of_c": 1.0, "formatted": "3 years"}],
            "origin_position": (0.0, 0.0, 0.0), "destination_position": (3.0, 1.0, 0.0),
            "route": {"path": [from_id if from_kind == "system" else f"phenomenon:{from_type}:{from_id}",
                               1500, to_id],
                      "distance_ly": 3.5, "positions": {"1500": (1.5, 0.5, 0.0)}},
        }

    def get_wiki_config(self):
        return self.wiki_config

    def auth_me(self, cookie_header):
        self.calls.append(("auth_me", cookie_header))
        return self.admin

    def upload_sector_to_wiki(self, cookie_header, db, sector_id, backend, path=None):
        self.calls.append(("upload", db, sector_id, backend, path))
        if self.action_error:
            raise self.action_error
        self.sectors[sector_id]["wiki_url"] = "https://wiki.example/sectors/fake"
        return {"url": "https://wiki.example/sectors/fake"}

    def generate_sector_neighborhood(self, cookie_header, sector_id, radius_ly=None):
        self.calls.append(("generate", sector_id))
        if self.action_error:
            raise self.action_error
        return {"generated": 4, "already_existed": 2, "candidates": 6}


_FAKED = ("get_sector", "get_sectors", "get_system", "get_phenomenon", "get_nav", "get_wiki_config",
          "auth_me", "upload_sector_to_wiki", "generate_sector_neighborhood")


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


def _log_in(client, fake, must_change=False):
    fake.admin = {"username": "admin", "must_change_credentials": must_change}
    client.set_cookie(SESSION_COOKIE_NAME, "token-value")


def _csrf(app, client):
    nonce = "n" * 43
    client.set_cookie(csrf.COOKIE_NAME, nonce)
    with app.app_context():
        return csrf._sign(nonce)


def _scene(html):
    match = re.search(r'<script type="application/json" id="starmap-data">(.*?)</script>', html, re.S)
    return json.loads(match.group(1))


# --- Sector page ------------------------------------------------------------------------

def test_sector_page_renders_badges_map_and_contents(client, fake):
    resp = client.get("/sector/5")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert "<title>Fake Sector - " in html
    assert '<h1 class="page-title">Fake Sector</h1>' in html
    assert re.search(r'<a href="/sectors" aria-current="page">Sectors</a>', html)
    crumbs = re.search(r'<nav class="breadcrumbs".*?</nav>', html, re.S).group(0)
    assert '<a href="/sectors">Sectors</a>' in crumbs and '<span aria-current="page">Fake Sector</span>' in crumbs
    assert "Cube edge 10.00 ly" in html and "2 systems" in html and "1 phenomenon" in html
    assert f'href="/galaxy.py?db={DB}&amp;quadrant=' in html
    assert re.search(r'<script type="module" src="/static/sectormap.js\?v=[^"]+"></script>', html)
    # Contents: nearest first, systems and phenomena, plain links.
    contents = html[html.index('id="sector-contents"'):]
    assert contents.index("Other") < contents.index("Alpha") < contents.index("Veil")
    assert f'<a href="/system.py?db={DB}&amp;id=1001">Alpha</a>' in contents
    assert f'<a href="/phenomenon.py?db={DB}&amp;type=nebula&amp;id=3">Veil</a>' in contents
    assert "Emission, 2.50 ly radius" in contents
    # Location's neighbor names link too.
    assert f'nearest: <a href="/system.py?db={DB}&amp;id=1002">Other</a> (1.0 ly)' in contents
    # Anonymous visitors get no forms at all.
    assert "<form method=\"post\"" not in html and "data-nav" not in html


def test_sector_map_entries_are_plain_links(client, fake):
    html = client.get("/sector/5").get_data(as_text=True)
    scene = _scene(html)
    assert {star["href"] for star in scene["stars"]} == {
        f"/system.py?db={DB}&id=1001", f"/system.py?db={DB}&id=1002"}
    assert scene["clouds"][0]["href"] == f"/phenomenon.py?db={DB}&type=nebula&id=3"
    assert scene["neighbors"][0]["href"] == "/sector/6"
    assert not any("navTarget" in entry for entry in scene["stars"] + scene["clouds"] + scene["neighbors"])
    assert '<li><a href="/sector/6">Next Door</a></li>' in html  # <noscript> list


def test_sector_contents_pager_uses_get_links(client, fake):
    fake.sectors[5]["systems"] = [_sector_system(4000 + i, f"S{i:03d}", float(i)) for i in range(55)]
    html = client.get("/sector/5").get_data(as_text=True)
    assert "Showing 1&ndash;50 of 56" in html
    assert 'href="/sector/5?contents_page=2#sector-contents"' in html
    page2 = client.get("/sector/5?contents_page=2").get_data(as_text=True)
    assert "Showing 51&ndash;56 of 56" in page2
    contents = page2[page2.index('id="sector-contents"'):]
    assert ">S054</a>" in contents and ">S000</a>" not in contents
    # The map still plots every system.
    assert len(_scene(page2)["stars"]) == 55


def test_sector_page_takes_database_from_config(client, fake):
    client.get("/sector/5?db=someone_elses")
    assert {call[1] for call in fake.calls if call[0] == "get_sector"} == {DB}


def test_unknown_sector_is_404(client, fake):
    resp = client.get("/sector/404")
    assert resp.status_code == 404
    assert "No such sector" in resp.get_data(as_text=True)
    assert client.get("/sector/abc").status_code == 404


def test_sector_name_is_escaped(client, fake):
    fake.sectors[5]["name"] = '"><img src=x onerror=alert(1)>'
    html = client.get("/sector/5").get_data(as_text=True)
    assert "<img src=x" not in html
    assert "&lt;img src=x onerror=alert(1)&gt;" in html


def test_sector_wiki_link_when_set(client, fake):
    fake.sectors[5]["wiki_url"] = "https://wiki.example/s"
    html = client.get("/sector/5").get_data(as_text=True)
    assert '<a href="https://wiki.example/s" target="_blank" rel="noopener noreferrer">View on Wiki</a>' in html


def test_admin_sees_forms_with_csrf_token(client, fake):
    _log_in(client, fake)
    html = client.get("/sector/5").get_data(as_text=True)
    forms = re.findall(r'<form method="post" action="/sector/5".*?</form>', html, re.S)
    assert len(forms) == 2
    for form in forms:
        assert f'name="{csrf.FIELD_NAME}"' in form
        assert 'name="db"' not in form and 'name="id"' not in form
    assert 'value="upload_wiki"' in forms[0] and 'value="wikijs"' in forms[0] and "mediawiki" not in forms[0]
    assert 'value="generate_neighborhood"' in forms[1]


def test_admin_with_stale_credentials_gets_no_generate_form(client, fake):
    _log_in(client, fake, must_change=True)
    html = client.get("/sector/5").get_data(as_text=True)
    assert 'value="generate_neighborhood"' not in html
    assert 'value="upload_wiki"' in html


def test_unplaced_sector_explains_why_it_cannot_generate(client, fake):
    fake.sectors[5] = _sector_detail(placed=False)
    _log_in(client, fake)
    html = client.get("/sector/5").get_data(as_text=True)
    assert 'value="generate_neighborhood"' not in html
    assert "never been placed in a galaxy" in html


def test_admin_post_without_csrf_token_is_rejected(client, fake):
    _log_in(client, fake)
    resp = client.post("/sector/5", data={"action": "generate_neighborhood"})
    assert resp.status_code == 400
    assert not [call for call in fake.calls if call[0] == "generate"]


def test_generate_neighborhood_posts_then_redirects_to_get(app, client, fake):
    _log_in(client, fake)
    token = _csrf(app, client)
    resp = client.post("/sector/5?contents_page=2",
                       data={"action": "generate_neighborhood", csrf.FIELD_NAME: token})
    assert resp.status_code == 303
    assert resp.headers["Location"] == "/sector/5"
    assert ("generate", 5) in fake.calls
    html = client.get("/sector/5").get_data(as_text=True)
    assert "Generated 4 new sector(s) (2 already existed, 6 candidate slot(s) within radius)." in html
    # The message is shown once.
    assert "Generated 4 new sector(s)" not in client.get("/sector/5").get_data(as_text=True)


def test_wiki_upload_posts_then_redirects_to_get(app, client, fake):
    _log_in(client, fake)
    token = _csrf(app, client)
    resp = client.post("/sector/5", data={
        "action": "upload_wiki", "backend": "wikijs", "path": " sectors/fake ", csrf.FIELD_NAME: token})
    assert resp.status_code == 303 and resp.headers["Location"] == "/sector/5"
    assert ("upload", DB, 5, "wikijs", "sectors/fake") in fake.calls
    html = client.get("/sector/5").get_data(as_text=True)
    assert "Uploaded to the wiki: https://wiki.example/sectors/fake" in html
    assert "View on Wiki" in html and 'value="upload_wiki"' not in html


def test_wiki_upload_conflict_is_shown_inline(app, client, fake):
    _log_in(client, fake)
    fake.action_error = apiclient.ApiError("planetGen API error (409): exists", status_code=409)
    resp = client.post("/sector/5", data={
        "action": "upload_wiki", "backend": "wikijs", csrf.FIELD_NAME: _csrf(app, client)})
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert '<p class="error" role="alert">A page already exists at that location.</p>' in html


def test_generate_error_is_shown_inline(app, client, fake):
    _log_in(client, fake)
    fake.action_error = apiclient.NotFoundError("Sector 5 has no galaxy placement")
    resp = client.post("/sector/5", data={"action": "generate_neighborhood", csrf.FIELD_NAME: _csrf(app, client)})
    assert resp.status_code == 200
    assert "Sector 5 has no galaxy placement" in resp.get_data(as_text=True)


def test_post_without_admin_session_changes_nothing(app, client, fake):
    resp = client.post("/sector/5", data={"action": "generate_neighborhood", csrf.FIELD_NAME: _csrf(app, client)})
    assert resp.status_code == 303 and resp.headers["Location"] == "/sector/5"
    assert not [call for call in fake.calls if call[0] in ("generate", "upload")]


def test_stale_credentials_cannot_generate(app, client, fake):
    _log_in(client, fake, must_change=True)
    resp = client.post("/sector/5", data={"action": "generate_neighborhood", csrf.FIELD_NAME: _csrf(app, client)})
    assert resp.status_code == 303
    assert ("generate", 5) not in fake.calls


def test_page_url_for_sector(app):
    with app.test_request_context("/"):
        assert page_url("sector", sector_id=5) == "/sector/5"
        assert page_url("sector", sector_id=5, contents_page=2) == "/sector/5?contents_page=2"


# --- NAV page ----------------------------------------------------------------------------

def _form(html, field):
    return re.search(rf'<form method="get" action="/nav"[^>]*>(?:(?!</form>).)*name="{field}".*?</form>',
                     html, re.S).group(0)


def test_nav_starts_with_a_sector_picker(client, fake):
    resp = client.get("/nav")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert re.search(r'<a href="/nav" aria-current="page">Nav</a>', html)
    form = _form(html, "from_sector")
    assert '<option value="5">Fake Sector</option>' in form and '<option value="9">Far Sector</option>' in form
    assert 'type="hidden"' not in form
    assert "csrf" not in form  # a GET form needs no token


def test_nav_second_step_lists_systems_as_endpoints(client, fake):
    html = client.get("/nav?from_sector=5").get_data(as_text=True)
    form = _form(html, "from")
    assert '<option value="system:1001">Alpha</option>' in form
    assert '<a href="/nav">&larr; Start over</a>' in html


def test_nav_carries_a_preset_destination_through_the_origin_picker(client, fake):
    html = client.get("/nav?to=nebula:3").get_data(as_text=True)
    assert '<input type="hidden" name="to" value="nebula:3">' in _form(html, "from_sector")
    html = client.get("/nav?to=nebula:3&from_sector=5").get_data(as_text=True)
    assert '<input type="hidden" name="to" value="nebula:3">' in _form(html, "from")
    assert '<a href="/nav?to=nebula:3">&larr; Start over</a>' in html


def test_nav_destination_pickers_for_a_system(client, fake):
    html = client.get("/nav?from=system:1001").get_data(as_text=True)
    assert "<title>Nav: Alpha - " in html
    assert f'From: <a href="/system.py?db={DB}&amp;id=1001">Alpha</a>' in html
    same = _form(html, "to")
    assert '<input type="hidden" name="from" value="system:1001">' in same
    assert '<option value="system:1002">Other</option>' in same
    assert '<option value="system:1001">' not in same
    cross = _form(html, "to_sector")
    assert '<option value="9">Far Sector</option>' in cross and 'value="5"' not in cross
    html = client.get("/nav?from=system:1001&to_sector=9").get_data(as_text=True)
    assert '<option value="system:2001">Distant</option>' in html
    assert '<a href="/nav?from=system:1001">&larr; Start over</a>' in html


def test_nav_unplaced_sector_offers_same_sector_only(client, fake):
    fake.sectors[5]["placed"] = False
    html = client.get("/nav?from=system:1001").get_data(as_text=True)
    assert 'name="to_sector"' not in html
    assert "no galaxy placement" in html


def test_nav_system_without_sector_is_unavailable(client, fake):
    html = client.get("/nav?from=system:3000").get_data(as_text=True)
    assert "NAV is not available: this system isn&#39;t assigned to a sector." in html


def test_nav_course_and_route(client, fake):
    resp = client.get("/nav?from=system:1001&to=system:1002")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert ("get_nav", DB, 1001, 1002, "system", "system", None, None) in fake.calls
    assert "3.25 ly" in html and "45.00&deg;" in html and "-2.50&deg;" in html and "3 years" in html
    assert "Same sector" in html
    assert f'To: <a href="/system.py?db={DB}&amp;id=1002">Other</a>' in html
    route = re.search(r'<ol class="nav-route">.*?</ol>', html, re.S).group(0)
    assert [name for name in re.findall(r">([^<]+)</a>", route)] == ["Alpha", "Waypoint", "Other"]
    assert "3 stops, 3.50 ly total." in html
    # The NAV map's points are plain links.
    assert f'<a class="navmap-point navmap-hop" href="/system.py?db={DB}&amp;id=1500">' in html
    assert "data-nav" not in html
    assert '<a href="/nav?from=system:1002&amp;to=system:1001">Reverse course</a>' in html


def test_nav_phenomenon_origin(client, fake):
    html = client.get("/nav?from=nebula:3&to=system:1002").get_data(as_text=True)
    assert ("get_nav", DB, 3, 1002, "phenomenon", "system", "nebula", None) in fake.calls
    assert f'From: <a href="/phenomenon.py?db={DB}&amp;type=nebula&amp;id=3">Veil</a>' in html
    route = re.search(r'<ol class="nav-route">.*?</ol>', html, re.S).group(0)
    assert f'href="/phenomenon.py?db={DB}&amp;type=nebula&amp;id=3">Veil</a>' in route


def test_nav_phenomenon_origin_picks_any_sector(client, fake):
    html = client.get("/nav?from=nebula:3").get_data(as_text=True)
    assert "Same-sector destination" not in html
    assert '<option value="5">Fake Sector</option>' in _form(html, "to_sector")


def test_nav_unplaced_phenomenon_is_unavailable(client, fake):
    fake.phenomena[("nebula", 3)]["galactic_radius_pc"] = None
    html = client.get("/nav?from=nebula:3").get_data(as_text=True)
    assert "this phenomenon has not been placed in the galaxy" in html


def test_nav_unsupported_pair_is_a_message_not_an_error(client, fake):
    fake.nav_error = apiclient.ApiError("planetGen API error (400): <different> sectors", status_code=400)
    resp = client.get("/nav?from=system:1001&to=system:2001")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert '<p class="error" role="alert">&lt;different&gt; sectors</p>' in html


def test_nav_api_failure_is_502(client, fake):
    fake.nav_error = apiclient.ApiError("planetGen API error (500): boom", status_code=500)
    assert client.get("/nav?from=system:1001&to=system:2001").status_code == 502


@pytest.mark.parametrize("query", [
    "from=abc", "from=system:x", "from=bogus:3", "from=system:1001&to=abc", "to=nope", "from_sector=zz",
    "from=system:99999",
])
def test_nav_bad_parameters_are_404(client, fake, query):
    resp = client.get(f"/nav?{query}")
    assert resp.status_code == 404
    assert "Traceback" not in resp.get_data(as_text=True)


def test_nav_bare_number_means_a_system(client, fake):
    html = client.get("/nav?from=1001&to=1002").get_data(as_text=True)
    assert ("get_nav", DB, 1001, 1002, "system", "system", None, None) in fake.calls
    assert "Reverse course" in html


@pytest.mark.parametrize("query,location", [
    ("from_id=1001", "/nav?from=system:1001"),
    ("from_id=1001&to_id=1002", "/nav?from=system:1001&to=system:1002"),
    ("from=3&from_kind=phenomenon&from_type=nebula&to=1002&to_kind=system",
     "/nav?from=nebula:3&to=system:1002"),
    ("to=3&to_kind=phenomenon&to_type=black_hole&from_sector=5", "/nav?to=black_hole:3&from_sector=5"),
])
def test_nav_old_parameter_style_redirects(client, fake, query, location):
    resp = client.get(f"/nav?{query}")
    assert resp.status_code == 301
    assert resp.headers["Location"] == location


def test_nav_url_helpers(app):
    with app.test_request_context("/"):
        assert page_url("nav") == "/nav"
        assert nav_url(endpoint("system", 12)) == "/nav?from=system:12"
        assert nav_url(None, endpoint("black_hole", 3)) == "/nav?to=black_hole:3"
        assert page_url("nav", **{"from": "system:1", "to": "nebula:2"}) == \
            "/nav?from=system:1&to=nebula:2"


# --- CGI shims -------------------------------------------------------------------------------

_NO_API = "http://127.0.0.1:9/api"


def test_sector_shim_redirects_get_and_post():
    result = run_page(_NO_API, "sector.py", query={"db": "x", "id": "5", "contents_page": "2"})
    assert result.status_code == 301
    assert result.headers["Location"] == "/sector/5?contents_page=2"
    result = run_page(_NO_API, "sector.py", method="POST",
                      body={"db": "x", "id": "7", "action": "generate_neighborhood"})
    assert result.status_code == 301
    assert result.headers["Location"] == "/sector/7"


@pytest.mark.parametrize("sector_id", ["", "abc", "0"])
def test_sector_shim_without_valid_id_goes_to_sectors(sector_id):
    result = run_page(_NO_API, "sector.py", query={"db": "x", "id": sector_id})
    assert result.status_code == 301
    assert result.headers["Location"] == "/sectors"


@pytest.mark.parametrize("params,location", [
    ({"db": "x"}, "/nav"),
    ({"db": "x", "from": "12"}, "/nav?from=system:12"),
    ({"db": "x", "from": "12", "to": "3", "to_kind": "phenomenon", "to_type": "neutron_star"},
     "/nav?from=system:12&to=neutron_star:3"),
    ({"db": "x", "to": "4", "from_sector": "9"}, "/nav?to=system:4&from_sector=9"),
    ({"db": "x", "from": "12", "to_sector": "9"}, "/nav?from=system:12&to_sector=9"),
    ({"db": "x", "from": "junk", "from_sector": "<x>"}, "/nav"),
])
def test_nav_shim_translates_parameters(params, location):
    for method in ("GET", "POST"):
        kwargs = {"query": params} if method == "GET" else {"method": "POST", "body": params}
        result = run_page(_NO_API, "nav.py", **kwargs)
        assert result.status_code == 301
        assert result.headers["Location"] == location


# --- Real database, in-process ------------------------------------------------------------

@pytest.fixture
def db_app(mysql_config):
    class RealConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "test-secret"

    application = create_app(RealConfig)
    application.testing = True
    return application


@pytest.fixture
def db_client(db_app):
    return db_app.test_client()


def _system_config(star_type):
    cfg = SystemConfig()
    cfg.STAR_TYPE = star_type
    cfg.PLANETS = False
    cfg.BINARY_SYSTEM = False
    return cfg


def _two_system_sector(mysql_config, name="Test Sector"):
    sector = SpaceSector(name, edge_ly=10.0)
    for star_type, position in (("G2V", (1.0, 1.0, 1.0)), ("M5V", (-2.0, 0.5, 3.0))):
        cfg = _system_config(star_type)
        sector.add_system(StarSystem(system_config=cfg), position=position, system_config=cfg)
    sector_id = _db.save_sector(sector, config=mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        rows = conn.execute("SELECT id, name FROM star_systems WHERE sector_id = ? ORDER BY id",
                            (sector_id,)).fetchall()
    finally:
        conn.close()
    return sector_id, rows


def test_real_sector_page_and_nav(db_client, mysql_config, monkeypatch):
    def no_http(*args, **kwargs):
        raise AssertionError("HTTP transport used inside a Flask request")
    monkeypatch.setattr(apiclient, "_http_transport", no_http)
    sector_id, systems = _two_system_sector(mysql_config)

    resp = db_client.get(f"/sector/{sector_id}")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    for row in systems:
        assert html_lib.escape(row["name"]) in html
    assert len(_scene(html)["stars"]) == 2

    html = db_client.get("/nav").get_data(as_text=True)
    assert f'<option value="{sector_id}">Test Sector</option>' in html
    html = db_client.get(f"/nav?from_sector={sector_id}").get_data(as_text=True)
    assert f'value="system:{systems[0]["id"]}"' in html

    origin, destination = systems[0]["id"], systems[1]["id"]
    resp = db_client.get(f"/nav?from=system:{origin}&to=system:{destination}")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert "Direct Course" in html and "Same sector" in html and "NAV Map" in html
    assert '<ol class="nav-route">' in html


def test_real_sector_page_unknown_id_is_404(db_client, mysql_config):
    _two_system_sector(mysql_config)
    resp = db_client.get("/sector/999999999")
    assert resp.status_code == 404
    assert "Traceback" not in resp.get_data(as_text=True)


def test_real_sector_page_with_galaxy_placement_renders_neighbor_indicators(db_client, mysql_config):
    from stellarObjects.galaxyGeometry import galactic_radius_pc, sector_position_pc

    edge_pc = 3.526
    shell_index, shell_slot_index = 5, 100
    position = sector_position_pc(shell_index, shell_slot_index, edge_pc)
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            sector_id = _db.insert_sector(conn, SpaceSector(name="Placed Sector"), galaxy_position={
                "center_x_pc": position[0], "center_y_pc": position[1], "center_z_pc": position[2],
                "galactic_radius_pc": galactic_radius_pc(position),
                "shell_index": shell_index, "shell_slot_index": shell_slot_index,
                "vertices_pc": {"inner": [], "outer": []},
            })
    finally:
        conn.close()

    resp = db_client.get(f"/sector/{sector_id}")
    assert resp.status_code == 200
    scene = _scene(resp.get_data(as_text=True))
    assert len(scene["neighbors"]) > 0
    for entry in scene["neighbors"]:
        assert entry["exists"] is False  # nothing else was ever placed
        assert entry["designation"]
        assert "href" not in entry


def test_real_sector_page_lists_and_maps_every_phenomenon_type(db_client, mysql_config):
    from stellarObjects.roguePlanetData import InterstellarComet, RoguePlanet
    from stellarObjects.supernovaRemnantData import SupernovaRemnant

    sector = SpaceSector(name="Phenomena Contents Sector", edge_ly=40.0)
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = False
    sector.add_system(StarSystem(system_config=cfg), system_config=cfg)
    remnant = SupernovaRemnant(SystemConfig())
    remnant.compact_remnant = None
    entries = [
        sector.add_phenomenon(remnant, "supernova-remnant"),
        sector.add_phenomenon(RoguePlanet(SystemConfig()), "rogue-planet"),
        sector.add_phenomenon(InterstellarComet(SystemConfig()), "comet"),
    ]
    sector_id = _db.save_sector(sector, config=mysql_config, galaxy_position={
        "center_x_pc": 500.0, "center_y_pc": 200.0, "center_z_pc": 10.0,
        "galactic_radius_pc": (500.0 ** 2 + 200.0 ** 2 + 10.0 ** 2) ** 0.5,
        "vertices_pc": {"inner": [], "outer": []},
    })

    html = db_client.get(f"/sector/{sector_id}").get_data(as_text=True)
    contents = html.split('id="contents-heading">Contents</h2>', 1)[1].split("</section>", 1)[0]
    conn = _db.get_connection(mysql_config)
    try:
        system_name = conn.execute("SELECT name FROM star_systems WHERE sector_id = ?",
                                   (sector_id,)).fetchone()["name"]
    finally:
        conn.close()
    for name in [system_name] + [entry.phenomenon.name for entry in entries]:
        assert html_lib.escape(name) in contents
    for label in ("Supernova Remnant", "Rogue Planet", "Interstellar Comet", "Star System"):
        assert label in contents
    kinds = {cloud["kind"] for cloud in _scene(html)["clouds"]}
    assert {"supernovaRemnant", "roguePlanet", "interstellarComet"} <= kinds


def test_real_sector_page_escapes_its_name(db_client, mysql_config):
    sector = SpaceSector('"><img src=x onerror=alert(1)>', edge_ly=10.0)
    cfg = _system_config("G2V")
    sector.add_system(StarSystem(system_config=cfg), position=(0.0, 0.0, 0.0), system_config=cfg)
    sector_id = _db.save_sector(sector, config=mysql_config)
    resp = db_client.get(f"/sector/{sector_id}")
    assert resp.status_code == 200
    assert "<img src=x onerror=alert(1)>" not in resp.get_data(as_text=True)


def test_real_admin_action_error_shows_on_the_page(db_client, mysql_config):
    """A real admin session reaches the admin API in-process: generating
    around an unplaced sector fails with the API's own message, shown on
    the page next to the form (no redirect)."""
    sector_id, _systems = _two_system_sector(mysql_config)
    adminAuth.bootstrap_control_schema(mysql_config)
    resp = db_client.post("/api/auth/login", json={
        "username": adminAuth.DEFAULT_ADMIN_USERNAME, "password": adminAuth.DEFAULT_ADMIN_PASSWORD})
    assert resp.status_code == 200
    resp = db_client.post("/api/auth/change-credentials", json={
        "current_password": adminAuth.DEFAULT_ADMIN_PASSWORD,
        "new_username": "sector-admin", "new_password": "a-strong-test-password-123"})
    assert resp.status_code == 200

    html = db_client.get(f"/sector/{sector_id}").get_data(as_text=True)
    assert "never been placed in a galaxy" in html  # no generate form for an unplaced sector
    token_value = _csrf(db_client.application, db_client)
    resp = db_client.post(f"/sector/{sector_id}",
                          data={"action": "generate_neighborhood", csrf.FIELD_NAME: token_value})
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert '<p class="error" role="alert">' in html
    assert "Traceback" not in html
