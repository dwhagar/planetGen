# tests/test_web_facilities.py

"""
Facilities on the web pages (UX.5, schema v42): the system page's
Facilities panel and body-list rows, the admin form (preview, save,
remove), the System Map's facility markers, a colony making its world
"Inhabited", and the sector page's Contents rows for facilities outside a
system.

Most tests fake the `apiclient` functions the views call, like
`test_web_system_phen.py`; the ones at the bottom run against a real
throwaway database through the in-process transport and are skipped
without a MySQL test server.
"""

import re

import pytest

from api.app import create_app
from api.authz import SESSION_COOKIE_NAME
from api.config import Config

import web  # noqa: F401 -- puts src/html/lib on sys.path
import apiclient  # noqa: E402
import queryDb  # noqa: E402
from stellarObjects import _db, adminAuth, physical_constants  # noqa: E402
from stellarObjects import facilities as facility_rules  # noqa: E402
from stellarObjects.config import SystemConfig  # noqa: E402
from stellarObjects.spaceSector import SpaceSector  # noqa: E402
from stellarObjects.systemData import StarSystem  # noqa: E402
from systemmap import render_system_map_panel  # noqa: E402
from web import csrf  # noqa: E402

DB = "planetgen_web_test"
AU_KM = physical_constants.AU_TO_KM


class _FakeConfig(Config):
    WEB_DATABASE = DB
    SESSION_COOKIE_SECURE = False
    SECRET_KEY = "test-secret"


def _moon():
    return {"id": 30, "name": "Luna", "planet_class": "D", "body_type": "t", "habitable": False,
            "inhabited": False, "zone": "e", "distance_km": 384_400.0, "period_years": 0.075, "gravity_g": 0.17,
            "mass_kg": 7.35e22, "hill_radius_km": 61_500.0,
            "orbital_index": 0, "radius_km": 1737.0, "position_x_km": 384_400.0, "position_y_km": 0.0,
            "life_chemical": None}


def _planet(planet_id, name, body_type, distance_au, moons=()):
    return {"id": planet_id, "name": name, "planet_class": "M" if body_type == "t" else "J", "body_type": body_type,
            "habitable": body_type == "t", "inhabited": False, "zone": "e", "distance_km": distance_au * AU_KM,
            "period_years": distance_au ** 1.5, "gravity_g": 1.0, "orbital_index": planet_id, "star_id": 1,
            "radius_km": 6371.0 if body_type == "t" else 69_911.0, "position_x_km": distance_au * AU_KM,
            "mass_kg": 5.97e24 if body_type == "t" else 1.898e27,
            "hill_radius_km": 1.5e6 if body_type == "t" else 5.3e7,
            "position_y_km": 0.0, "life_chemical": None, "moons": list(moons)}


def _system():
    return {
        "id": 5, "name": "Sol", "sector_id": None, "quadrant": None, "location": None,
        "is_binary": 0, "binary_type": None, "binary_configuration": None,
        "binary_mutual_position_x_km": None, "binary_mutual_position_y_km": None,
        "binary_mutual_position_z_km": None, "wikijs_url": None, "mediawiki_url": None,
        "stars": [{"id": 1, "role": "single", "name": "Sol", "star_type": "G2V", "mass_kg": 1.989e30,
                   "radius_km": 696_000.0, "temperature_k": 5778.0, "luminosity_w": 3.828e26,
                   "heliosphere_radius_km": 120 * AU_KM}],
        "planets": [_planet(10, "Terra", "t", 1.0, [_moon()]), _planet(11, "Jove", "g", 5.2)],
        "belts": [{"id": 20, "star_id": 1, "distance_km": 2.8 * AU_KM, "lower_limit_km": 2.2 * AU_KM,
                   "upper_limit_km": 3.3 * AU_KM, "density": "typical", "composition_summary": "rock",
                   "orbital_index": 3}],
        "comets": [], "sector_siblings": [], "nearest_neighbors": [],
    }


def _facility(facility_id, name, kind, placement, host_type, host_id, host_name, distance_km=None, **extra):
    row = {
        "id": facility_id, "name": name, "kind": kind, "placement": placement, "host_type": host_type,
        "host_id": host_id, "host_name": host_name, "star_system_id": 5, "description": None,
        "orbit_distance_km": distance_km, "orbit_period_years": None, "orbital_speed_kms": None,
        "orbit_phase_deg": None, "center_x_pc": None, "center_y_pc": None, "center_z_pc": None,
    }
    if distance_km is not None:
        row.update(orbit_period_years=0.0027, orbital_speed_kms=3.07, orbit_phase_deg=45.0)
    row.update(extra)
    return row


def _facilities():
    return [
        _facility(1, "High <Yard>", "starbase", "orbital", "planet", 11, "Jove", distance_km=4.2e5),
        _facility(2, "New Hope", "colony", "terrestrial", "planet", 10, "Terra"),
        _facility(3, "Sunwatch", "outpost", "orbital", "star", 1, "Sol", distance_km=0.5 * AU_KM),
        _facility(4, "Rockpile", "mining-colony", "asteroid", "asteroid_belt", 20, None),
        _facility(5, "Tidewatch", "station", "orbital", "moon", 30, "Luna", distance_km=1.0e4),
    ]


class FakeData:
    """Stands in for the `apiclient` functions these pages call."""

    def __init__(self):
        self.system = _system()
        self.facilities = _facilities()
        self.sections = {"overview": "", "stars": {}, "planets": {}, "moons": {}, "belts": {}, "comets": {}}
        self.admin = None
        self.orbit_error = None
        self.create_error = None
        self.calls = []

    def get_system(self, db, system_id):
        if int(system_id) != 5:
            raise apiclient.NotFoundError(f"no such system: {system_id}")
        return self.system

    def get_system_sections(self, db, system_id):
        return self.sections

    def get_system_facilities(self, db, system_id):
        self.calls.append(("get_system_facilities", db, system_id))
        return self.facilities

    def get_wiki_config(self):
        return {"wikijs": False, "mediawiki": False}

    def auth_me(self, cookie_header):
        return self.admin

    def get_facility_orbit(self, db, host_type, host_id, distance_km=None):
        self.calls.append(("orbit", db, host_type, host_id, distance_km))
        if self.orbit_error:
            raise self.orbit_error
        return {"distance_km": distance_km or 3 * 696_000.0, "period_years": 0.25, "orbital_speed_kms": 42.125}

    def create_facility(self, cookie_header, db, body):
        self.calls.append(("create", db, body))
        if self.create_error:
            raise self.create_error
        return 99

    def delete_facility(self, cookie_header, db, facility_id):
        self.calls.append(("delete", db, facility_id))


_FAKED = ("get_system", "get_system_sections", "get_system_facilities", "get_wiki_config", "auth_me",
          "get_facility_orbit", "create_facility", "delete_facility")


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
    fake.admin = {"username": "boss", "must_change_credentials": False}
    client.set_cookie(SESSION_COOKIE_NAME, "token-value")


def _csrf(app, client):
    nonce = "n" * 43
    client.set_cookie(csrf.COOKIE_NAME, nonce)
    session = client.get_cookie(SESSION_COOKIE_NAME)
    with app.app_context():
        return csrf._sign(nonce, session.value if session else "")


def _panel(html):
    return re.search(r'<section class="panel" id="facilities".*?</section>', html, re.S).group(0)


def _post(app, client, **form):
    form.setdefault(csrf.FIELD_NAME, _csrf(app, client))
    return client.post("/system/5", data=form)


_STAR_LIMITS = facility_rules.orbit_limits(696_000.0, 120 * AU_KM)
_STAR_OUTPOST = {"name": "Far Point", "host": "star:1", "placement": "orbital", "kind": "outpost",
                 "orbit_step": str(facility_rules.step_for_distance(0.75 * AU_KM, *_STAR_LIMITS))}


# --- Listing --------------------------------------------------------------------------

def test_system_page_lists_facilities(client, fake):
    html = client.get("/system/5").get_data(as_text=True)
    panel = _panel(html)
    assert "High &lt;Yard&gt;" in panel and "<Yard>" not in html
    for text in ("Starbase", "Jove (planet)", "In orbit", "420,000 km", "3.07 km/s",
                 "Mining colony", "Asteroid belt", "Among the asteroids", "On the surface", "Luna (moon)"):
        assert text in panel, text
    # Visitors get no form and no Remove buttons.
    assert 'name="facility_action"' not in html
    # Each facility also sits in its host's row of the system list.
    terra_row = html[html.index('<span class="body-name">Terra</span>'):]
    terra_row = terra_row[:terra_row.index("</details>")]
    assert "1 facility" in terra_row and "New Hope" in terra_row and "High" not in terra_row
    assert '<ul class="facility-list" aria-label="Facilities">' in html
    assert ("get_system_facilities", DB, 5) in fake.calls


def test_system_page_without_facilities_shows_no_panel(client, fake):
    fake.facilities = []
    html = client.get("/system/5").get_data(as_text=True)
    assert 'id="facilities"' not in html
    assert "facility-list" not in html


# --- Admin form -----------------------------------------------------------------------

def test_admin_sees_the_form_with_every_host(client, fake):
    _as_admin(client, fake)
    panel = _panel(client.get("/system/5").get_data(as_text=True))
    for value in ("star:1", "planet:10", "planet:11", "moon:30", "asteroid_belt:20"):
        assert f'<option value="{value}"' in panel
    # Name, then placement, then host (ADM.9); each host says which
    # placements it takes, for static/facilityform.js to filter on.
    assert panel.index('name="name"') < panel.index('name="placement"') < panel.index('name="host"')
    assert '<option value="planet:10" data-placements="terrestrial orbital"' in panel
    assert '<option value="planet:11" data-placements="orbital"' in panel
    assert '<option value="asteroid_belt:20" data-placements="asteroid"' in panel
    assert "Asteroid belt from 2.2 AU to 3.3 AU" in panel
    # The orbit is a slider, shown for the default "in orbit" placement.
    assert 'type="range" id="facility-orbit-step" name="orbit_step" min="0" max="1000"' in panel
    assert 'name="distance"' not in panel and 'name="distance_unit"' not in panel
    assert re.search(r'<div class="search-field facility-orbit" data-facility-orbit>', panel)
    assert 'value="preview">Preview</button>' in panel
    assert 'value="save">Save</button>' in panel
    assert panel.count('name="facility_action" value="remove"') == 5
    assert "<script" not in panel


def test_preview_shows_the_orbit_without_saving(app, client, fake):
    _as_admin(client, fake)
    resp = _post(app, client, facility_action="preview", **_STAR_OUTPOST)
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert ("orbit", DB, "star", 1, pytest.approx(0.75 * AU_KM, rel=0.01)) in fake.calls
    assert not [call for call in fake.calls if call[0] == "create"]
    panel = _panel(html)
    assert "Allowed by the placement rules." in panel
    assert "42.12 km/s" in panel or "42.13 km/s" in panel
    assert "3 months" in panel or "91 days" in panel
    # The form keeps what was typed.
    assert 'value="Far Point"' in panel
    assert re.search(r'<option value="star:1"[^>]* selected>', panel)
    assert re.search(r'<option value="orbital" selected>', panel)


def test_preview_reports_a_rule_violation(app, client, fake):
    _as_admin(client, fake)
    html = _post(app, client, facility_action="preview", name="Floaters", host="planet:11",
                 placement="terrestrial", kind="colony").get_data(as_text=True)
    assert "Planet: Jove (gas giant) can&#39;t take a facility placed on its surface" in html
    assert not [call for call in fake.calls if call[0] in ("orbit", "create")]


def test_surface_and_belt_placements_hide_the_orbit_slider(app, client, fake):
    _as_admin(client, fake)
    panel = _panel(_post(app, client, facility_action="preview", name="Second Hope", host="moon:30",
                         placement="terrestrial", kind="colony").get_data(as_text=True))
    assert "Allowed by the placement rules." in panel
    assert '<div class="search-field facility-orbit" data-facility-orbit hidden>' in panel
    assert not [call for call in fake.calls if call[0] == "orbit"]

    panel = _panel(_post(app, client, facility_action="preview", name="Rockpile", host="asteroid_belt:20",
                         placement="asteroid", kind="mining-colony").get_data(as_text=True))
    assert "random spot in the belt, between 2.2 AU and 3.3 AU from its star" in panel
    assert '<div class="search-field facility-orbit" data-facility-orbit hidden>' in panel
    resp = _post(app, client, facility_action="save", name="Rockpile", host="asteroid_belt:20",
                 placement="asteroid", kind="mining-colony", orbit_step="900")
    assert resp.status_code == 303
    (create,) = [call for call in fake.calls if call[0] == "create"]
    assert "distance_km" not in create[2]


def test_slider_runs_from_the_surface_to_the_sphere_of_influence(app, client, fake):
    _as_admin(client, fake)
    _post(app, client, facility_action="preview", name="Low", host="planet:10", placement="orbital",
          kind="station", orbit_step="0")
    _post(app, client, facility_action="preview", name="High", host="planet:10", placement="orbital",
          kind="station", orbit_step="1000")
    orbits = [call[4] for call in fake.calls if call[0] == "orbit"]
    assert orbits == [pytest.approx(6371.0 * 1.01), pytest.approx(1.5e6)]


def test_preview_reports_the_apis_orbit_error(app, client, fake):
    _as_admin(client, fake)
    fake.orbit_error = apiclient.ApiError(
        "planetGen API error (400): orbit distance 1000 km is inside the host (radius 69911 km)", status_code=400)
    html = _post(app, client, facility_action="preview", name="Low", host="planet:11", placement="orbital",
                 kind="station", orbit_step="0").get_data(as_text=True)
    assert "Not allowed: orbit distance 1000 km is inside the host (radius 69911 km)." in html
    assert "Allowed by the placement rules." not in html


def test_save_posts_to_the_api_and_redirects(app, client, fake):
    _as_admin(client, fake)
    resp = _post(app, client, facility_action="save", description="A relay", **_STAR_OUTPOST)
    assert resp.status_code == 303
    assert resp.headers["Location"] == "/system/5?facility=added#facilities"
    (create,) = [call for call in fake.calls if call[0] == "create"]
    assert create[1] == DB
    assert create[2] == {"name": "Far Point", "kind": "outpost", "placement": "orbital", "host_type": "star",
                         "host_id": 1, "distance_km": pytest.approx(0.75 * AU_KM, rel=0.01),
                         "description": "A relay"}
    assert "Facility added." in client.get(resp.headers["Location"]).get_data(as_text=True)


def test_save_of_a_colony_sends_no_orbit(app, client, fake):
    _as_admin(client, fake)
    resp = _post(app, client, facility_action="save", name="Second Hope", host="moon:30",
                 placement="terrestrial", kind="colony", orbit_step="500")
    assert resp.status_code == 303
    (create,) = [call for call in fake.calls if call[0] == "create"]
    assert create[2] == {"name": "Second Hope", "kind": "colony", "placement": "terrestrial", "host_type": "moon",
                         "host_id": 30}


def test_save_errors_show_next_to_the_form(app, client, fake):
    _as_admin(client, fake)
    html = _post(app, client, facility_action="save", **dict(_STAR_OUTPOST, name="")).get_data(as_text=True)
    assert "Give the facility a name." in _panel(html)
    html = _post(app, client, facility_action="save",
                 **dict(_STAR_OUTPOST, host="asteroid_belt:20")).get_data(as_text=True)
    assert "can&#39;t take a facility placed in orbit around it; choose another host." in html
    html = _post(app, client, facility_action="save", **dict(_STAR_OUTPOST, host="planet:999")).get_data(as_text=True)
    assert "Choose where the facility goes." in html
    assert not [call for call in fake.calls if call[0] == "create"]

    fake.create_error = apiclient.ApiError("planetGen API error (403): default credentials must be changed",
                                           status_code=403)
    resp = _post(app, client, facility_action="save", **_STAR_OUTPOST)
    assert resp.status_code == 200
    assert "Change the admin username and password the installer set first." in _panel(resp.get_data(as_text=True))


def test_remove_deletes_only_this_systems_facilities(app, client, fake):
    _as_admin(client, fake)
    resp = _post(app, client, facility_action="remove", facility_id="3")
    assert resp.status_code == 303
    assert resp.headers["Location"] == "/system/5?facility=removed#facilities"
    assert ("delete", DB, 3) in fake.calls

    fake.calls.clear()
    html = _post(app, client, facility_action="remove", facility_id="777").get_data(as_text=True)
    assert "That facility is not in this system." in html
    assert not [call for call in fake.calls if call[0] == "delete"]


def test_facility_posts_need_an_admin_and_a_token(app, client, fake):
    token = _csrf(app, client)
    assert client.post("/system/5", data={csrf.FIELD_NAME: token, "facility_action": "save",
                                          **_STAR_OUTPOST}).status_code == 403
    _as_admin(client, fake)
    assert client.post("/system/5", data={"facility_action": "save", **_STAR_OUTPOST}).status_code == 400
    assert not [call for call in fake.calls if call[0] == "create"]


# --- System Map -----------------------------------------------------------------------

def test_system_map_draws_facilities_at_their_hosts():
    system = _system()
    html = render_system_map_panel(system, system["stars"], system["planets"], system["belts"],
                                   facilities=_facilities())
    markers = re.findall(r'<g class="sysmap-facility[^"]*"[^>]*data-kind="facility"[^>]*>', html)
    # The moon's in Terra's moon scene, the rest in the system scene, and
    # Terra's own colony in both (around Terra, and on the moon scene's
    # center).
    assert len(markers) == 6
    assert 'data-name="High &lt;Yard&gt;"' in html
    moon_scene = "".join(re.findall(r'<svg[^>]*data-scene="planet-10"[^>]*>.*?</svg>', html, re.S))
    assert 'data-name="Tidewatch"' in moon_scene and 'data-name="New Hope"' in moon_scene
    assert html.count('class="sysmap-facility-orbit"') == 3
    assert "sysmap-legend-facility" in html

    plain = render_system_map_panel(system, system["stars"], system["planets"], system["belts"])
    assert "sysmap-facility" not in plain and "sysmap-legend-facility" not in plain


# --- Sector page ----------------------------------------------------------------------

def test_sector_contents_list_facilities_outside_systems(app):
    from web.sector_page import _contents

    sector = {"systems": [], "center_x_pc": 0.0, "center_y_pc": 50.0, "center_z_pc": 0.0, "phenomena": [
        {"id": 7, "type": "asteroid_field", "name": "The Shoals", "descriptor": "dense", "radius_ly": 1.0,
         "distance_ly": 2.5, "class": "C3", "octant": None, "nearest": []},
    ]}
    facilities = [
        _facility(8, "Waypoint <1>", "station", "standalone", "space", 3, "Harbor", star_system_id=None,
                  center_x_pc=0.0, center_y_pc=50.0 + 1.0 / 3.26156, center_z_pc=0.0),
        _facility(9, "Shoal Camp", "outpost", "asteroid", "asteroid_field", 7, "The Shoals", star_system_id=None),
    ]
    with app.test_request_context("/sector/3"):
        rows, _map = _contents(sector, facilities)
    by_name = {row["name"]: row for row in rows}
    waypoint = by_name["Waypoint <1>"]
    assert waypoint["type"] == "Facility (Station)"
    assert waypoint["url"] is None
    assert waypoint["details"] == "Stand-alone station, parked in open space"
    assert waypoint["distance_ly"] == pytest.approx(1.0, rel=1e-3)
    camp = by_name["Shoal Camp"]
    assert camp["distance_ly"] == 2.5
    assert str(camp["location"]) == 'On <a href="/phenomenon/asteroid_field/7">The Shoals</a>'
    assert rows[0]["name"] == "Waypoint <1>"  # nearest the center first


# --- Real database, in-process ----------------------------------------------------------

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


def _system_with_terrestrial(mysql_config):
    """A saved single-star system with at least one terrestrial planet
    (retried, since generation is random). Returns `(system_id, ids)`."""
    for _ in range(200):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.PLANETS = True
        cfg.MOONS = False
        cfg.BINARY_SYSTEM = False
        system = StarSystem(system_config=cfg)
        if any(p.body_type == "t" for p in system.planets):
            system_id = _db.save_system(system, cfg, config=mysql_config)
            break
    else:
        pytest.fail("could not generate a system with a terrestrial planet")
    conn = _db.get_connection(mysql_config)
    try:
        ids = {
            "star": conn.execute("SELECT id FROM stars WHERE star_system_id = ?", (system_id,)).fetchone()["id"],
            "terrestrial": conn.execute("SELECT id FROM planets WHERE star_system_id = ? AND body_type = 't'"
                                        " ORDER BY id LIMIT 1", (system_id,)).fetchone()["id"],
        }
    finally:
        conn.close()
    return system_id, ids


def test_real_colony_makes_its_world_inhabited(db_app, mysql_config):
    system_id, ids = _system_with_terrestrial(mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        detail = queryDb.system_detail(conn, system_id)
        before = next(p for p in detail["planets"] if p["id"] == ids["terrestrial"])
        with conn:
            _db.add_facility(conn, "New Hope", "colony", "terrestrial", "planet", ids["terrestrial"])
        detail = queryDb.system_detail(conn, system_id)
    finally:
        conn.close()
    planet = next(p for p in detail["planets"] if p["id"] == ids["terrestrial"])
    assert planet["inhabited"] is True
    assert before["habitable"] == planet["habitable"]  # a colony changes nothing else
    others = [p for p in detail["planets"] if p["id"] != ids["terrestrial"]]
    assert all(p["inhabited"] == (p["life_stage"] == "technological_civilization") for p in others)

    html = db_app.test_client().get(f"/system/{system_id}").get_data(as_text=True)
    row = html[html.index(f'<span class="body-name">{planet["name"]}</span>'):]
    row = row[:row.index("</summary>")]
    assert '<span class="flag flag-yes">Inhabited</span>' in row
    assert "New Hope" in _panel(html)


def _login_fresh_admin(client, mysql_config):
    """Logs in as an admin whose credentials are current (the API's write
    routes refuse the installer's default ones)."""
    _username, first_password = adminAuth.bootstrap_control_schema(mysql_config)
    assert client.post("/api/auth/login", json={"username": "admin", "password": first_password}).status_code == 200
    assert client.post("/api/auth/change-credentials", json={
        "current_password": first_password, "new_username": "boss", "new_password": "a-long-new-password-1",
    }).status_code == 200


def test_real_admin_previews_saves_and_removes(db_app, mysql_config):
    system_id, ids = _system_with_terrestrial(mysql_config)
    client = db_app.test_client()
    _login_fresh_admin(client, mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        star = conn.execute("SELECT radius_km, heliosphere_radius_km FROM stars WHERE id = ?",
                            (ids["star"],)).fetchone()
    finally:
        conn.close()
    limits = facility_rules.orbit_limits(star["radius_km"], star["heliosphere_radius_km"])
    form = {"name": "Sunwatch", "host": f"star:{ids['star']}", "placement": "orbital", "kind": "outpost",
            "orbit_step": str(facility_rules.step_for_distance(0.5 * AU_KM, *limits))}

    resp = client.post(f"/system/{system_id}", data={csrf.FIELD_NAME: _csrf(db_app, client),
                                                     "facility_action": "preview", **form})
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert "Allowed by the placement rules." in html
    # Half an AU around a Sun-like star: about 0.35 years, about 42 km/s.
    speed = float(re.search(r"at ([\d.,]+) km/s", _panel(html)).group(1).replace(",", ""))
    assert 35 < speed < 50
    conn = _db.get_connection(mysql_config)
    try:
        assert queryDb.facilities_for_system(conn, system_id) == []
    finally:
        conn.close()

    resp = client.post(f"/system/{system_id}", data={csrf.FIELD_NAME: _csrf(db_app, client),
                                                     "facility_action": "save", **form})
    assert resp.status_code == 303
    conn = _db.get_connection(mysql_config)
    try:
        (saved,) = queryDb.facilities_for_system(conn, system_id)
    finally:
        conn.close()
    assert saved["orbit_distance_km"] == pytest.approx(0.5 * AU_KM, rel=0.01)
    assert saved["orbital_speed_kms"] == pytest.approx(speed, rel=0.01)
    page = client.get(resp.headers["Location"]).get_data(as_text=True)
    assert "Facility added." in page and "Sunwatch" in _panel(page)
    assert 'data-kind="facility"' in page

    resp = client.post(f"/system/{system_id}", data={csrf.FIELD_NAME: _csrf(db_app, client),
                                                     "facility_action": "remove", "facility_id": str(saved["id"])})
    assert resp.status_code == 303
    assert "Facility removed." in client.get(resp.headers["Location"]).get_data(as_text=True)
    conn = _db.get_connection(mysql_config)
    try:
        assert queryDb.facilities_for_system(conn, system_id) == []
    finally:
        conn.close()


def test_real_sector_page_lists_stand_alone_facilities(db_app, mysql_config):
    sector = SpaceSector("Harbor", edge_ly=13.0)
    sector_id = _db.save_sector(sector, config=mysql_config, galaxy_position={
        "center_x_pc": 0.0, "center_y_pc": 50.0, "center_z_pc": 0.0, "galactic_radius_pc": 50.0,
    })
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            _db.add_facility(conn, "Waypoint", "station", "standalone", "space", sector_id, offset_ly=(1.0, 2.0, 2.0))
    finally:
        conn.close()
    html = db_app.test_client().get(f"/sector/{sector_id}").get_data(as_text=True)
    row = re.search(r"<tr>\s*<td>Waypoint</td>.*?</tr>", html, re.S).group(0)
    assert "Facility (Station)" in row
    assert "Stand-alone station, parked in open space" in row
    assert "3.00 ly" in row or "3 ly" in row, row
