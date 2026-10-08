# tests/test_web_sector_nav.py

"""
The Flask-served sector page (`/sector/<id>`, `web/sector_page.py`) and
NAV page (`/nav`, `web/nav_page.py`), plus the old `sector.py`/`nav.py`
URLs that now redirect to them.

Most tests fake the data layer by monkeypatching the `apiclient`
functions the views call (see `test_web_pages.py`). The tests at the
bottom use a real throwaway database through the in-process transport and
are skipped without a MySQL test server.
"""

import json
import re

import pytest
from markupsafe import escape

from planetgen.web.app import create_app
from planetgen.api.authz import SESSION_COOKIE_NAME
from planetgen.api.config import Config

from planetgen.web.lib import apiclient  # noqa: E402
from planetgen.db import store  # noqa: E402
from planetgen.admin import auth as adminAuth
from planetgen import tuning
from planetgen.generation.config import SystemConfig  # noqa: E402
from planetgen.galaxy.sector import SpaceSector  # noqa: E402
from planetgen.generation.system import StarSystem  # noqa: E402
from planetgen.web import csrf, generate_page, jobs, sector_page  # noqa: E402
from planetgen.web.helpers import page_url  # noqa: E402
from planetgen.web.nav_page import endpoint, nav_url  # noqa: E402


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
        "ring_index": 5 if placed else None, "layer_index": 1 if placed else None,
        "ring_slot_index": 20 if placed else None,
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
            "direction_pc": (1.0, 0.0, 0.0), "ring_index": 5, "layer_index": 1, "ring_slot_index": 21,
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
        self.nav_galaxy_scope = False
        self.action_error = None
        self.bright_stars = {}

    def get_sector(self, db, sector_id):
        self.calls.append(("get_sector", db, sector_id))
        try:
            return self.sectors[int(sector_id)]
        except (KeyError, ValueError):
            raise apiclient.NotFoundError(f"No such sector: {sector_id}")

    def get_galaxy_shape(self, db):
        self.calls.append(("get_galaxy_shape", db))
        return None

    def get_bright_stars_in_cell(self, db, ring_index, layer_index, ring_slot_index):
        self.calls.append(("get_bright_stars_in_cell", ring_index, layer_index, ring_slot_index))
        return self.bright_stars.get((ring_index, layer_index, ring_slot_index), [])

    def get_sector_facilities(self, db, sector_id):
        self.calls.append(("get_sector_facilities", db, sector_id))
        return []

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
            "scope": "galaxy" if self.nav_galaxy_scope else "sector",
            "direct": {"distance_ly": 3.25, "bearing_deg": 45.2, "mark_deg": 357.5,
                       "elevation_deg": -2.5, "frame": "sector"},
            "warp_times": [{"warp_factor": 1, "velocity_multiple_of_c": 1.0, "formatted": "3 years"}],
            "fold_times": [{"fold_factor": 4, "velocity_multiple_of_c": 256.0, "formatted": "4 days"}],
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

    estimate = {"sectors": 4, "summary": "About 280 KB and 3 s for 4 sectors.", "refused": False, "refusal": None}

    def generate_sector_neighborhood(self, cookie_header, sector_id, radius_ly=None, estimate_only=False):
        self.calls.append(("estimate" if estimate_only else "generate", sector_id))
        if self.action_error:
            raise self.action_error
        if estimate_only:
            return {"generated": 0, "already_existed": 2, "candidates": 6, "estimate": self.estimate}
        return {"generated": 4, "already_existed": 2, "candidates": 6}


_FAKED = ("get_sector", "get_galaxy_shape", "get_sector_facilities", "get_bright_stars_in_cell", "get_sectors", "get_system", "get_phenomenon", "get_nav", "get_wiki_config",
          "auth_me", "upload_sector_to_wiki", "generate_sector_neighborhood")


@pytest.fixture
def fake(monkeypatch):
    data = FakeData()
    for name in _FAKED:
        monkeypatch.setattr(apiclient, name, getattr(data, name))
    # The map's first frame of galaxy tiles (the sector page embeds the Galaxy Map's engine).
    monkeypatch.setattr(sector_page, "fetch_tiles", lambda db, keys, known_stamp=None: {
        "stamp": "0" * 16, "tiles": {}, "density": {}, "edge_pc": 3.066, "has_shape": False})
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
    session = client.get_cookie(SESSION_COOKIE_NAME)
    with app.app_context():
        return csrf._sign(nonce, session.value if session else "")  # bound to the login session


def _page_scene(html):
    """The map's scene data a sector page embeds (`#galaxymap3d-data`)."""
    match = re.search(r'<script type="application/json" id="galaxymap3d-data">(.*?)</script>', html, re.S)
    return json.loads(match.group(1))


def _scene(client, path):
    """The sector's scene JSON (`/sector/<id>/scene`, what the map draws) for a
    sector page's `path`, query and all."""
    base, _, query = path.partition("?")
    response = client.get(f"{base}/scene" + (f"?{query}" if query else ""))
    assert response.status_code == 200
    return response.get_json()


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
    assert "Cube edge 3.07 pc (10 ly)" in html and "2 systems" in html and "1 phenomenon" in html
    assert 'href="/galaxy?quadrant=' in html
    assert re.search(r'<script type="module" src="/static/galaxymap3d.js\?v=[^"]+"></script>', html)
    # The map is the Galaxy Map locked to this sector (MAP.68).
    assert 'id="galaxymap3d-canvas"' in html and "galaxymap3d-crumbs" not in html
    assert json.loads(re.search(r'id="galaxymap3d-data">(.*?)</script>', html, re.S).group(1))["pinned"] == {
        "ring": 5, "layer": 1, "slot": 20}
    # Contents: nearest first, systems and phenomena, plain links.
    contents = html[html.index('id="sector-contents"'):]
    assert contents.index("Other") < contents.index("Alpha") < contents.index("Veil")
    assert '<a href="/system/1001">Alpha</a>' in contents
    assert '<a href="/phenomenon/nebula/3">Veil</a>' in contents
    assert "Emission, 2.5 ly radius" in contents
    # Location's neighbor names link too.
    assert 'nearest: <a href="/system/1002">Other</a> (1.0 ly)' in contents
    # Anonymous visitors get no forms at all.
    assert "<form method=\"post\"" not in html and "data-nav" not in html


def test_sector_map_entries_are_plain_links(client, fake):
    scene = _scene(client, "/sector/5")
    assert {star["href"] for star in scene["stars"]} == {
        "/system/1001", "/system/1002"}
    assert scene["clouds"][0]["href"] == "/phenomenon/nebula/3"
    assert scene["neighbors"][0]["href"] == "/sector/6"
    assert not any("navTarget" in entry for entry in scene["stars"] + scene["clouds"] + scene["neighbors"])


def test_sector_scene_json_is_the_pages_own_scene_for_the_galaxy_map(client, fake):
    """MAP.66: /sector/<id>/scene is the scene the map draws, with where the
    sector sits (centerPc) and its half edge (halfEdgePc, in parsecs), for
    the Galaxy Map to open the sector in place."""
    resp = client.get("/sector/5/scene")
    assert resp.status_code == 200 and resp.mimetype == "application/json"
    assert resp.headers["Cache-Control"] == "no-store"
    scene = resp.get_json()
    assert len(scene["centerPc"]) == 3
    assert scene["halfEdgePc"] == pytest.approx(fake.sectors[5]["edge_mpc"] / 2000)
    assert {star["href"] for star in scene["stars"]} == {"/system/1001", "/system/1002"}


def test_sector_scene_names_each_endpoint_and_leaves_the_nav_buttons_to_the_page(client, fake):
    scene = client.get("/sector/5/scene").get_json()
    assert {star["endpoint"] for star in scene["stars"]} == {"system:1001", "system:1002"}
    assert not any("nav" in entry for entry in scene["stars"] + scene["clouds"])


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
    assert len(_scene(client, "/sector/5?contents_page=2")["stars"]) == 55


def test_sector_contents_list_systems_then_phenomena_then_rogues(client, fake):
    """UX.24: systems first, then other phenomena, then rogue planets, each
    nearest first; the rogue group is one full-width row whose members
    show their octant and location."""
    rogue = {"type": "rogue_planet", "descriptor": "jupiter", "class": "J", "radius_ly": 0.0,
             "offset_x_ly": 0.1, "offset_y_ly": 0.0, "offset_z_ly": 0.0,
             "octant": "+x+y+z", "nearest": [{"id": 1002, "name": "Other", "distance_ly": 0.5}]}
    fake.sectors[5]["phenomena"] += [
        {**rogue, "id": 21, "name": "Drifter", "distance_ly": 0.1},
        {**rogue, "id": 22, "name": "Wanderer", "distance_ly": 0.2, "octant": "-x-y-z"},
    ]
    fake.sectors[5]["phenomena"][0]["distance_ly"] = 9.0
    html = client.get("/sector/5").get_data(as_text=True)
    contents = html[html.index('id="sector-contents"'):]
    assert contents.index(">Other<") < contents.index(">Alpha<") < contents.index(">Veil<") < contents.index(">Drifter<")
    group = re.search(r'<tr class="contents-group-row">\s*<td colspan="6">(.*?)</details>', contents, re.S).group(1)
    assert ">Drifter</a>" in group and ">Wanderer</a>" in group
    assert "+x+y+z" in group and "-x-y-z" in group
    assert 'Nearest: <a href="/system/1002">Other</a>' in group
    # UX.25: "Show on map" is a small map icon beside each name, its
    # words in the aria-label and tooltip.
    assert re.search(r'>Drifter</a> <button type="button" class="icon-btn" data-map-target="rogue_planet:21" '
                     r'title="Show on map" aria-label="Show Drifter on the map" hidden><svg class="icon" '
                     r'aria-hidden="true" focusable="false"><use href="/static/icons\.svg\?v=[^"]+#show-on-map">',
                     group)


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


@pytest.mark.parametrize("wiki_url", ["javascript:alert(1)", "data:text/html,x", "//evil.example/x"])
def test_sector_wiki_link_never_links_a_non_http_url(client, fake, wiki_url):
    # One saved before the API checked it is left out of the page.
    fake.sectors[5]["wiki_url"] = wiki_url
    html = client.get("/sector/5").get_data(as_text=True)
    assert "View on Wiki" not in html
    assert f'href="{wiki_url}"' not in html


def test_admin_sees_forms_with_csrf_token(client, fake):
    _log_in(client, fake)
    html = client.get("/sector/5").get_data(as_text=True)
    every_form = re.findall(r'<form method="post" action="/sector/5".*?</form>', html, re.S)
    for form in every_form:
        assert f'name="{csrf.FIELD_NAME}"' in form
        assert 'name="db"' not in form and 'name="id"' not in form
    # The Edit panel's Regenerate and Delete forms (ADM.8) aside:
    forms = [form for form in every_form if 'name="edit_action"' not in form]
    assert len(forms) == 2 and len(every_form) == 4
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


def test_generate_neighborhood_shows_the_estimate_first(app, client, fake):
    """PERF.3: the first press only works out the size and time; the
    page asks before anything is generated."""
    _log_in(client, fake)
    resp = client.post("/sector/5", data={"action": "generate_neighborhood", csrf.FIELD_NAME: _csrf(app, client)})
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert ("estimate", 5) in fake.calls and ("generate", 5) not in fake.calls
    assert "About 280 KB and 3 s for 4 sectors." in html
    assert 'name="estimate_ok" value="1"' in html and "Generate these 4 sectors" in html


def test_a_refused_neighborhood_offers_generate_anyway(app, client, fake):
    """ADM.33: no room on the disk stops the plain button, but the admin
    can generate anyway."""
    _log_in(client, fake)
    fake.estimate = dict(fake.estimate, refused=True, refusal="Refused: not enough disk.")
    resp = client.post("/sector/5", data={"action": "generate_neighborhood", csrf.FIELD_NAME: _csrf(app, client)})
    html = resp.get_data(as_text=True)
    assert "Refused: not enough disk." in html and "recorded in the activity log" in html
    assert 'name="estimate_ok" value="1"' in html and 'name="generate_anyway" value="1"' in html
    assert "Generate these 4 sectors anyway" in html
    assert ("generate", 5) not in fake.calls


def test_generate_anyway_starts_the_neighborhood_job_and_logs_it(app, client, fake, monkeypatch, tmp_path):
    monkeypatch.setenv("PLANETGEN_JOBS_DIR", str(tmp_path / "jobs"))
    real_start = jobs.start_job
    monkeypatch.setattr(jobs, "start_job", lambda kind, title, steps, **kw: real_start(kind, title, steps,
                                                                                       spawn=False, **kw))
    events = []
    monkeypatch.setattr(sector_page.activity_log, "event", lambda *args, **kwargs: events.append((args, kwargs)))
    _log_in(client, fake)
    resp = client.post("/sector/5", data={"action": "generate_neighborhood", "estimate_ok": "1",
                                          "generate_anyway": "1", csrf.FIELD_NAME: _csrf(app, client)})
    assert resp.status_code == 303
    names = [args[1] for args, _kwargs in events if args[0] == "GEN"]
    assert names == ["job.start", "job.generate_anyway"]


def test_generate_neighborhood_starts_a_job_then_redirects_to_get(app, client, fake, monkeypatch, tmp_path):
    """ADM.11: the confirmed run is a Generate page job (detached from the
    request), not a request that only ends when the generation does."""
    monkeypatch.setenv("PLANETGEN_JOBS_DIR", str(tmp_path / "jobs"))
    started = []
    real_start = jobs.start_job
    monkeypatch.setattr(jobs, "start_job", lambda kind, title, steps, **kw: started.append((kind, steps, kw))
                        or real_start(kind, title, steps, spawn=False, **kw))
    _log_in(client, fake)
    token = _csrf(app, client)
    resp = client.post("/sector/5?contents_page=2",
                       data={"action": "generate_neighborhood", "estimate_ok": "1", csrf.FIELD_NAME: token})
    assert resp.status_code == 303
    assert resp.headers["Location"] == "/sector/5"
    assert ("generate", 5) not in fake.calls
    kind, steps, kw = started[0]
    assert kind == "galaxy" and kw["database"] == DB
    assert steps[0]["label"] == generate_page.MATH_CHECK_LABEL
    assert steps[-1]["argv"][-5:] == ["galaxy", "--center-sector", "5", "--radius-pc",
                                      str(tuning.DEFAULT_GENERATE_RADIUS_PC)]
    html = client.get("/sector/5").get_data(as_text=True)
    assert "keeps running if you close this page" in html
    # The message is shown once.
    assert "keeps running if you close" not in client.get("/sector/5").get_data(as_text=True)
    # A second one while it runs is refused on the page.
    resp = client.post("/sector/5", data={"action": "generate_neighborhood", "estimate_ok": "1",
                                          csrf.FIELD_NAME: _csrf(app, client)})
    assert resp.status_code == 200 and "is still running" in resp.get_data(as_text=True)


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
    assert 'From: <a href="/system/1001">Alpha</a>' in html
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
    assert "3.25 ly" in html and "045 mark 358" in html and "3 years" in html
    assert "Sector Local Frame" in html and "Fold 4" in html and "4 days" in html
    assert "Same sector" in html
    assert 'To: <a href="/system/1002">Other</a>' in html
    route = re.search(r'<ol class="nav-route">.*?</ol>', html, re.S).group(0)
    assert [name for name in re.findall(r">([^<]+)</a>", route)] == ["Alpha", "Waypoint", "Other"]
    assert "3 stops, 1.07 pc (3.5 ly) total." in html
    # The NAV map's points are plain links.
    assert '<a class="navmap-point navmap-hop" href="/system/1500">' in html
    assert "data-nav=" not in html
    assert '<a href="/nav?from=system:1002&amp;to=system:1001">Reverse course</a>' in html


def test_nav_result_links_to_the_galaxy_map(client, fake):
    """The course panel offers "Show on Galaxy Map", which carries both
    endpoints to the map (the drill-down's section 9.4)."""
    html = client.get("/nav?from=system:1001&to=system:1002").get_data(as_text=True)
    assert 'href="/galaxy?course=system:1001,system:1002"' in html
    assert ">Show on Galaxy Map</a>" in html


def test_galaxy_course_is_the_waypoints_in_parsecs(client, fake):
    """`nav_page.galaxy_course` hands the map galaxy-frame parsecs; a
    course inside one sector names that sector instead."""
    from planetgen.physics.units import ly_to_pc
    from planetgen.web.nav_page import galaxy_course

    app = client.application
    with app.test_request_context("/galaxy"):
        same_sector = galaxy_course("system:1001", "system:1002")
        assert same_sector["scope"] == "sector"
        assert same_sector["points"] == []
        assert same_sector["sector"]["id"] == 5
        assert same_sector["navUrl"] == "/nav?from=system:1001&to=system:1002"

        fake.nav_galaxy_scope = True
        course = galaxy_course("system:1001", "system:2001")
        assert course["scope"] == "galaxy"
        assert [p["name"] for p in course["points"]] == ["Alpha", "Waypoint", "Distant"]
        assert [p["role"] for p in course["points"]] == ["origin", "hop", "destination"]
        assert course["points"][-1]["x"] == pytest.approx(ly_to_pc(3.0))
        assert course["points"][0]["url"] == "/system/1001"

        # A pair that can't be navigated together simply isn't drawn.
        fake.nav_error = apiclient.ApiError("planetGen API error (400): different sectors", status_code=400)
        assert galaxy_course("system:1001", "system:2001") is None


def test_nav_phenomenon_origin(client, fake):
    html = client.get("/nav?from=nebula:3&to=system:1002").get_data(as_text=True)
    assert ("get_nav", DB, 3, 1002, "phenomenon", "system", "nebula", None) in fake.calls
    assert 'From: <a href="/phenomenon/nebula/3">Veil</a>' in html
    route = re.search(r'<ol class="nav-route">.*?</ol>', html, re.S).group(0)
    assert 'href="/phenomenon/nebula/3">Veil</a>' in route


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


# --- Old CGI URLs ---------------------------------------------------------------------------

def test_old_sector_url_redirects(client):
    result = client.get("/sector.py?db=x&id=5&contents_page=2")
    assert result.status_code == 301
    assert result.headers["Location"] == "/sector/5?contents_page=2"


@pytest.mark.parametrize("sector_id", ["", "abc", "0"])
def test_old_sector_url_without_valid_id_goes_to_sectors(client, sector_id):
    result = client.get(f"/sector.py?db=x&id={sector_id}")
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
def test_old_nav_url_translates_parameters(client, params, location):
    result = client.get("/nav.py", query_string=params)
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
    sector_id = store.save_sector(sector, config=mysql_config)
    conn = store.get_connection(mysql_config)
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
        assert escape(row["name"]) in html
    assert len(_scene(db_client, f"/sector/{sector_id}")["stars"]) == 2

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
    from planetgen.galaxy.geometry import galactic_radius_pc, sector_position_pc

    edge_pc = 3.526
    address = (5, 1, 20)
    position = sector_position_pc(*address, edge_pc)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            sector_id = store.insert_sector(conn, SpaceSector(name="Placed Sector"), galaxy_position={
                "center_x_pc": position[0], "center_y_pc": position[1], "center_z_pc": position[2],
                "galactic_radius_pc": galactic_radius_pc(position),
                "ring_index": address[0], "layer_index": address[1], "ring_slot_index": address[2],
            })
    finally:
        conn.close()

    resp = db_client.get(f"/sector/{sector_id}")
    assert resp.status_code == 200
    scene = _scene(db_client, f"/sector/{sector_id}")
    assert len(scene["neighbors"]) > 0
    for entry in scene["neighbors"]:
        assert entry["exists"] is False  # nothing else was ever placed
        assert entry["designation"]
        assert "href" not in entry


def test_real_sector_page_lists_and_maps_every_phenomenon_type(db_client, mysql_config):
    from planetgen.generation.phenomena.rogue import InterstellarComet, RoguePlanet
    from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant

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
    sector_id = store.save_sector(sector, config=mysql_config, galaxy_position={
        "center_x_pc": 500.0, "center_y_pc": 200.0, "center_z_pc": 10.0,
        "galactic_radius_pc": (500.0 ** 2 + 200.0 ** 2 + 10.0 ** 2) ** 0.5,
    })

    html = db_client.get(f"/sector/{sector_id}").get_data(as_text=True)
    contents = html.split('id="contents-heading">Contents</h2>', 1)[1].split("</section>", 1)[0]
    conn = store.get_connection(mysql_config)
    try:
        system_name = conn.execute("SELECT name FROM star_systems WHERE sector_id = ?",
                                   (sector_id,)).fetchone()["name"]
    finally:
        conn.close()
    for name in [system_name] + [entry.phenomenon.name for entry in entries]:
        assert escape(name) in contents
    for label in ("Supernova Remnant", "Rogue Planet", "Interstellar Comet", "Star System"):
        assert label in contents
    kinds = {cloud["kind"] for cloud in _scene(db_client, f"/sector/{sector_id}")["clouds"]}
    assert {"supernovaRemnant", "roguePlanet", "interstellarComet"} <= kinds


def test_real_sector_page_escapes_its_name(db_client, mysql_config):
    sector = SpaceSector('"><img src=x onerror=alert(1)>', edge_ly=10.0)
    cfg = _system_config("G2V")
    sector.add_system(StarSystem(system_config=cfg), position=(0.0, 0.0, 0.0), system_config=cfg)
    sector_id = store.save_sector(sector, config=mysql_config)
    resp = db_client.get(f"/sector/{sector_id}")
    assert resp.status_code == 200
    assert "<img src=x onerror=alert(1)>" not in resp.get_data(as_text=True)


def test_real_admin_action_error_shows_on_the_page(db_client, mysql_config):
    """A real admin session reaches the admin API in-process: generating
    around an unplaced sector fails with the API's own message, shown on
    the page next to the form (no redirect)."""
    sector_id, _systems = _two_system_sector(mysql_config)
    _username, first_password = adminAuth.bootstrap_control_schema(mysql_config)
    resp = db_client.post("/api/auth/login", json={
        "username": adminAuth.DEFAULT_ADMIN_USERNAME, "password": first_password})
    assert resp.status_code == 200
    resp = db_client.post("/api/auth/change-credentials", json={
        "current_password": first_password,
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


# --- NAV links and pick mode (MAP.21, design doc sections 9.1-9.2) ----------------

def test_sector_map_page_offers_start_here_and_end_here_without_a_pick(client, fake):
    page = _page_scene(client.get("/sector/5").get_data(as_text=True))
    assert page["pick"] is None and page["navUrl"] == "/nav"
    assert "pick-banner" not in client.get("/sector/5").get_data(as_text=True)


def nav_url_for(app, origin, destination):
    with app.test_request_context():
        return nav_url(origin=origin, destination=destination)


def test_pick_destination_mode(app, client, fake):
    html = client.get("/sector/5?pick=to&from=system:12").get_data(as_text=True)
    banner = re.search(r'<p class="pick-banner".*?</p>', html, re.S).group(0)
    assert "Choosing a destination" in banner and "End Here" in banner
    assert f'href="{escape(nav_url_for(app, "system:12", None))}">Cancel</a>' in banner
    page = _page_scene(html)
    assert page["pick"] == "to" and page["pickOther"] == "system:12"
    assert page["navUrl"] == "/nav"


def test_pick_start_mode_without_the_other_end(app, client, fake):
    html = client.get("/sector/5?pick=from").get_data(as_text=True)
    assert "Choosing a start" in html
    assert f'href="{escape(nav_url_for(app, None, None))}">Cancel</a>' in html
    page = _page_scene(html)
    assert page["pick"] == "from" and page["pickOther"] is None


@pytest.mark.parametrize("query", ["pick=sideways", "pick=to&from=planet:3", "pick=to&from=system:x"])
def test_bad_pick_mode_is_ignored(client, fake, query):
    resp = client.get(f"/sector/5?{query}")
    assert resp.status_code == 200
    html = resp.get_data(as_text=True)
    assert "pick-banner" not in html
    assert _page_scene(html)["pick"] is None


# --- Bright stars waiting in an unfilled neighbor -----------------------------------

def test_unfilled_neighbor_lists_its_waiting_bright_stars(client, fake):
    fake.sectors[5]["neighbors"].append({
        "direction_pc": (0.0, 1.0, 0.0), "ring_index": 5, "layer_index": 1, "ring_slot_index": 22,
        "designation": "ABD", "exists": False,
    })
    fake.bright_stars[(5, 1, 22)] = [
        {"star_type": f"B{i}V", "luminosity_sol": 9000.0 - i} for i in range(7)
    ]
    scene = _scene(client, "/sector/5")
    existing, unfilled = scene["neighbors"]
    assert "brightStars" not in existing
    assert unfilled["brightStarCount"] == 7
    assert unfilled["brightStars"][0] == "B0V, 9,000 L\u2609"
    assert len(unfilled["brightStars"]) == 5
    # Only the unfilled cell is asked about.
    assert [c for c in fake.calls if c[0] == "get_bright_stars_in_cell"] == [
        ("get_bright_stars_in_cell", 5, 1, 22)]


def test_bright_stars_fail_open(client, fake, monkeypatch):
    fake.sectors[5]["neighbors"].append({
        "direction_pc": (0.0, 1.0, 0.0), "ring_index": 5, "layer_index": 1, "ring_slot_index": 22,
        "designation": "ABD", "exists": False,
    })

    def broken(*args):
        raise apiclient.ApiError("down")
    monkeypatch.setattr(apiclient, "get_bright_stars_in_cell", broken)
    assert "brightStars" not in _scene(client, "/sector/5")["neighbors"][1]


# --- Map picks on the NAV page (MAP.22) and "Show on Galaxy Map" (MAP.25) ---------------

def _map_picks(html):
    match = re.search(r'<section class="panel" aria-labelledby="map-picks-heading">.*?</section>', html, re.S)
    return match.group(0) if match else ""


def test_nav_origin_step_offers_map_picks(client, fake):
    picks = _map_picks(client.get("/nav").get_data(as_text=True))
    assert "Or pick a start on a map" in picks
    assert 'href="/galaxy?pick=from">Pick on Galaxy Map</a>' in picks
    assert "Pick in this sector" not in picks  # no other endpoint yet
    picks = _map_picks(client.get("/nav?to=system:2001").get_data(as_text=True))
    assert 'href="/sector/9?pick=from&amp;to=system:2001">Pick in this sector</a>' in picks
    assert 'href="/galaxy?pick=from&amp;to=system:2001">Pick on Galaxy Map</a>' in picks
    picks = _map_picks(client.get("/nav?to=nebula:3").get_data(as_text=True))
    assert "Pick in this sector" not in picks  # a phenomenon has no sector
    assert 'href="/galaxy?pick=from&amp;to=nebula:3"' in picks


def test_nav_destination_step_offers_map_picks(client, fake):
    picks = _map_picks(client.get("/nav?from=system:1001").get_data(as_text=True))
    assert "Or pick a destination on a map" in picks
    assert 'href="/sector/5?pick=to&amp;from=system:1001">Pick in this sector</a>' in picks
    assert 'href="/galaxy?pick=to&amp;from=system:1001">Pick on Galaxy Map</a>' in picks


def test_nav_course_has_no_map_picks(client, fake):
    assert _map_picks(client.get("/nav?from=system:1001&to=system:1002").get_data(as_text=True)) == ""


def test_nav_bad_preset_destination_is_404(client, fake):
    assert client.get("/nav?to=system:999").status_code == 404


def test_sector_page_shows_it_on_the_galaxy_map(client, fake):
    from planetgen.galaxy.geometry import provisional_sector_designation

    html = client.get("/sector/5").get_data(as_text=True)
    designation = provisional_sector_designation(5, 1, 20)
    assert f'<a href="/galaxy?sector={designation}">Show on Galaxy Map</a>' in html
    fake.sectors[5] = _sector_detail(placed=False)
    assert "Show on Galaxy Map" not in client.get("/sector/5").get_data(as_text=True)


def test_galaxy_map_redirects(client, fake):
    from planetgen.galaxy.geometry import provisional_sector_designation

    designation = provisional_sector_designation(5, 1, 20)
    resp = client.get("/sector/5/galaxy")
    assert resp.status_code == 302 and resp.headers["Location"].endswith(f"/galaxy?sector={designation}")
    resp = client.get("/system/1001/galaxy")
    assert resp.status_code == 302 and resp.headers["Location"].endswith(f"/galaxy?sector={designation}")
    assert client.get("/system/3000/galaxy").headers["Location"].endswith("/galaxy")
    assert client.get("/sector/404/galaxy").status_code == 404
    fake.sectors[5] = _sector_detail(placed=False)
    assert client.get("/sector/5/galaxy").headers["Location"].endswith("/galaxy")


# --- Bookmarks (MAP.23) and the NAV page's Bookmarks select (MAP.22) ------------------

def _bookmark_button(html):
    match = re.search(r'<button type="button" class="btn btn-small btn-bookmark" data-bookmark-toggle[^>]*>[^<]*</button>',
                      html, re.S)
    return match.group(0) if match else ""


def test_sector_page_has_a_bookmark_button(client, fake):
    from planetgen.galaxy.geometry import provisional_sector_designation

    html = client.get("/sector/5").get_data(as_text=True)
    assert re.search(r'<script type="module" src="/static/bookmarks.js\?v=[^"]+"></script>', html)
    button = _bookmark_button(html)
    designation = provisional_sector_designation(5, 1, 20)
    for attribute in (f'data-bookmark-db="{DB}"', 'data-bookmark-kind="sector"',
                      f'data-bookmark-value="{designation}"', 'data-bookmark-name="Fake Sector"',
                      'data-bookmark-url="/sector/5"', 'data-bookmark-sector-id="5"', 'aria-pressed="false"'):
        assert attribute in button
    # Shown and wired by the script; bookmarks live in the browser.
    assert " hidden>" in button and "☆ Bookmark" in button
    # No galaxy address: keyed by its page instead.
    fake.sectors[5] = _sector_detail(placed=False)
    assert 'data-bookmark-value="/sector/5"' in _bookmark_button(client.get("/sector/5").get_data(as_text=True))


def _nav_bookmarks(html):
    match = re.search(r'<form class="search-form nav-bookmarks" data-bookmarks-nav[^>]*>.*?</form>', html, re.S)
    return match.group(0) if match else ""


def test_nav_origin_step_has_a_bookmarks_select(client, fake):
    html = client.get("/nav?to=system:2001").get_data(as_text=True)
    assert re.search(r'<script type="module" src="/static/bookmarks.js\?v=[^"]+"></script>', html)
    form = _nav_bookmarks(_map_picks(html))
    for attribute in (f'data-bookmark-db="{DB}"', 'data-nav-url="/nav"', 'data-pick="from"',
                      'data-keep-name="to"', 'data-keep-value="system:2001"'):
        assert attribute in form
    assert " hidden>" in form  # until bookmarks.js finds bookmarks to offer
    assert '<select name="bookmark"></select>' in form
    assert ">Use as start</button>" in form
    assert 'data-keep-value=""' in _nav_bookmarks(client.get("/nav").get_data(as_text=True))


def test_nav_destination_step_has_a_bookmarks_select(client, fake):
    form = _nav_bookmarks(client.get("/nav?from=system:1001").get_data(as_text=True))
    assert 'data-pick="to"' in form and 'data-keep-name="from"' in form
    assert 'data-keep-value="system:1001"' in form
    assert ">Use as destination</button>" in form


def test_nav_course_offers_bookmarks_for_either_end(client, fake):
    """NAV.40: once both ends are set, a bookmark can replace either."""
    html = client.get("/nav?from=system:1001&to=system:1002").get_data(as_text=True)
    assert re.search(r'<script type="module" src="/static/bookmarks.js\?v=[^"]+"></script>', html)
    group = re.search(r'<section class="panel"[^>]*data-bookmarks-nav-group hidden>.*?</section>', html, re.S)
    assert group, "the Change an End section, hidden until bookmarks.js fills it"
    forms = re.findall(r'<form class="search-form nav-bookmarks" data-bookmarks-nav[^>]*>.*?</form>', group.group(0), re.S)
    assert len(forms) == 2
    assert 'data-pick="from"' in forms[0] and 'data-keep-name="to"' in forms[0]
    assert 'data-keep-value="system:1002"' in forms[0] and "New start" in forms[0]
    assert 'data-pick="to"' in forms[1] and 'data-keep-name="from"' in forms[1]
    assert 'data-keep-value="system:1001"' in forms[1] and "New destination" in forms[1]


def test_sector_pick_mode_has_a_bookmarks_menu_that_keeps_the_pick(client, fake):
    """NAV.40: the sector page offers bookmarks while picking, carrying
    the pick for static/bookmarks.js."""
    html = client.get("/sector/5?pick=to&from=system:1001").get_data(as_text=True)
    menu = re.search(r'<details class="bookmarks-menu pick-bookmarks"[^>]*>', html, re.S)
    assert menu
    for attribute in ('data-bookmarks-menu', 'data-pick="to"', 'data-keep-name="from"',
                      'data-keep-value="system:1001"', 'data-nav-url="/nav"', f'data-bookmark-db="{DB}"'):
        assert attribute in menu.group(0)
    assert "pick-bookmarks" not in client.get("/sector/5").get_data(as_text=True)
