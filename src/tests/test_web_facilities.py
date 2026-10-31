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
from html import escape

import pytest

from planetgen.web.app import create_app
from planetgen.api.authz import SESSION_COOKIE_NAME
from planetgen.api.config import Config

from planetgen.describe import model  # noqa: E402
from planetgen.web.lib import apiclient  # noqa: E402
from planetgen.db import query  # noqa: E402
from planetgen.db import store  # noqa: E402
from planetgen.admin import auth as adminAuth
from planetgen.physics import constants
from planetgen.population import facilities as facility_rules  # noqa: E402
from planetgen.generation.config import SystemConfig  # noqa: E402
from planetgen.galaxy.sector import SpaceSector  # noqa: E402
from planetgen.generation.system import StarSystem  # noqa: E402
from planetgen.web.maps.systemmap import render_system_map_panel  # noqa: E402
from planetgen.web import csrf  # noqa: E402
from tests.publicids import pid

DB = "planetgen_web_test"
AU_KM = constants.AU_TO_KM


class _FakeConfig(Config):
    WEB_DATABASE = DB
    SESSION_COOKIE_SECURE = False
    SECRET_KEY = "test-secret"


def _moon():
    return {"id": 30, "name": "Luna", "planet_class": "D", "body_type": "t", "habitable": False,
            "habitable_class": False, "inhabited": False, "zone": "e", "distance_km": 384_400.0, "period_years": 0.075, "gravity_g": 0.17,
            "mass_kg": 7.35e22, "hill_radius_km": 61_500.0,
            "orbital_index": 0, "radius_km": 1737.0, "position_x_km": 384_400.0, "position_y_km": 0.0,
            "life_class": None}


def _planet(planet_id, name, body_type, distance_au, moons=()):
    return {"id": planet_id, "name": name, "planet_class": "M" if body_type == "t" else "J", "body_type": body_type,
            "habitable": body_type == "t", "habitable_class": body_type == "t", "inhabited": False, "zone": "e", "distance_km": distance_au * AU_KM,
            "period_years": distance_au ** 1.5, "gravity_g": 1.0, "orbital_index": planet_id, "star_id": 1,
            "radius_km": 6371.0 if body_type == "t" else 69_911.0, "position_x_km": distance_au * AU_KM,
            "mass_kg": 5.97e24 if body_type == "t" else 1.898e27,
            "hill_radius_km": 1.5e6 if body_type == "t" else 5.3e7,
            "position_y_km": 0.0, "life_class": None, "moons": list(moons)}


def _system():
    return {
        "id": 5, "name": "Sol", "sector_id": None, "quadrant": None, "location": None,
        "is_binary": 0, "binary_type": None, "binary_configuration": None,
        "binary_mutual_position_x_km": None, "binary_mutual_position_y_km": None,
        "binary_mutual_position_z_km": None,
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
        "type_id": None, "type_name": None, "type_icon": None, "affiliation_polity_id": None,
        "affiliation_name": None, "fields": {},
        "ratings": {name: None for name in ("crime", "housing", "resources", "maintenance", "health")},
        "rating_words": {name: None for name in ("crime", "housing", "resources", "maintenance", "health")},
        "overall": None, "ratings_source": None, "ratings_note": None, "ratings_set_at": None,
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
        self.admin = None
        self.orbit_error = None
        self.create_error = None
        self.types = []
        self.owner = None
        self.suggested = ["crime"]
        self.calls = []

    def get_facility_types(self, db, active_only=True):
        return self.types

    def get_object_children(self, db, ref, limit=None, after=None):
        """The system page's tree (UX.173) reads the children of the system and of each body; this one has none."""
        return {"ref": ref, "kind": ref.partition(":")[0], "total": 0, "items": [], "next": None}

    def get_object_sheet(self, db, ref):
        return {"kind": ref.partition(":")[0], "groups": []}

    def get_system_habitable_species(self, db, system_id):
        return []

    def get_system_text(self, db, system_id, fmt):
        return {"id": system_id, "format": fmt, "content": ""}

    def get_system_owner(self, db, system_id):
        return self.owner

    def update_facility(self, cookie_header, db, facility_id, body):
        self.calls.append(("update", db, facility_id, body))
        return {}

    def suggest_facility_ratings(self, cookie_header, db, facility_id):
        self.calls.append(("suggest", db, facility_id))
        return self.suggested

    def get_system(self, db, system_id):
        if int(system_id) != 5:
            raise apiclient.NotFoundError(f"no such system: {system_id}")
        return self.system

    def get_system_description(self, db, system_id):
        return model.document("system", "system:1", "Sol", [model.section("system", "system:1", "Sol")])

    def get_system_facilities(self, db, system_id):
        self.calls.append(("get_system_facilities", db, system_id))
        return self.facilities

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


_FAKED = ("get_system", "get_system_description", "get_system_facilities", "auth_me", "get_object_children",
          "get_object_sheet", "get_system_habitable_species", "get_system_text",
          "get_facility_orbit", "create_facility", "delete_facility", "get_facility_types", "get_system_owner",
          "update_facility", "suggest_facility_ratings")


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
    """The Facilities table and the Place-a-facility dialog (ADM.34 moved the form into the Admin menu)."""
    section = re.search(r'<section class="panel" id="facilities".*?</section>', html, re.S)
    dialog = re.search(r'<sl-dialog id="admin-facility-place".*?</sl-dialog>', html, re.S)
    return (section.group(0) if section else "") + (dialog.group(0) if dialog else "")


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
    # UX.176: a facility is listed on its host's own page, not in a row of the system page
    # (tests/test_item_pages.py::test_a_planets_page_lists_the_facilities_it_hosts_and_its_chips).
    assert 'class="body-name"' not in html and 'class="facility-list"' not in html
    assert ("get_system_facilities", DB, "5") in fake.calls


def test_the_panel_shows_type_affiliation_and_ratings_worst_first(client, fake):
    fake.facilities = [
        _facility(1, "Calm", "outpost", "orbital", "star", 1, "Sol", distance_km=1e8, overall="green",
                  ratings={"crime": "green", "housing": "green", "resources": None, "maintenance": None,
                           "health": "green"},
                  rating_words={"crime": "Low", "housing": "Ample", "resources": None, "maintenance": None,
                                "health": "Healthy"}),
        _facility(2, "Plain", "outpost", "orbital", "star", 1, "Sol", distance_km=2e8),
        _facility(3, "Troubled", "colony", "terrestrial", "planet", 10, "Terra", type_id=4, type_name="Farm <world>",
                  type_icon="tree", affiliation_name="Free Traders", fields={"hectares": 40.5, "irrigated": True},
                  overall="red",
                  ratings={"crime": "yellow", "housing": "red", "resources": None, "maintenance": None,
                           "health": "green"},
                  rating_words={"crime": "Elevated", "housing": "Overcrowded", "resources": None,
                                "maintenance": None, "health": "Healthy"}),
    ]
    panel = _panel(client.get("/system/5").get_data(as_text=True))
    assert panel.index("Troubled") < panel.index("Calm") < panel.index("Plain")        # worst first, unassessed last
    assert "Farm &lt;world&gt;" in panel and "<world>" not in panel
    assert '<sl-icon name="tree"' in panel and "Free Traders" in panel
    assert "hectares: 40.5, irrigated: yes" in panel
    # Each chip carries a shape, its dimension's name and its word, and a label that does not need colour.
    assert 'aria-label="Housing: Red, overcrowded"' in panel and "x-octagon-fill" in panel
    assert 'aria-label="Crime: Yellow, elevated"' in panel and "exclamation-triangle-fill" in panel
    assert 'aria-label="Health: Green, healthy"' in panel and "check-circle-fill" in panel
    troubled = panel[panel.index("Troubled"):panel.index("Calm")]
    assert troubled.index("Overcrowded") < troubled.index("Elevated") < troubled.index("Healthy")   # worst chip first
    assert "Resources" not in troubled                                          # unassessed ones are left out


def test_system_page_without_facilities_shows_no_panel(client, fake):
    fake.facilities = []
    html = client.get("/system/5").get_data(as_text=True)
    assert 'id="facilities"' not in html
    assert "facility-list" not in html


# --- Admin form -----------------------------------------------------------------------

def test_admin_sees_the_form_with_every_host(client, fake):
    _as_admin(client, fake)
    page = client.get("/system/5").get_data(as_text=True)
    panel = _panel(page)
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
    # Each facility is removed from the Admin menu's "Remove a facility" submenu.
    assert page.count('name="facility_action" value="remove"') == 5
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
    assert "42.1 km/s" in panel
    assert "91.3 days" in panel or "91.2 days" in panel
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
    assert ("delete", DB, "3") in fake.calls

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
    # center). The moon scene is drawn at two body sizes (MAP.308), so its two facilities are there twice.
    assert len(markers) == 8
    assert 'data-name="High &lt;Yard&gt;"' in html
    moon_scene = "".join(re.findall(r'<svg[^>]*data-scene="planet-10"[^>]*>.*?</svg>', html, re.S))
    assert 'data-name="Tidewatch"' in moon_scene and 'data-name="New Hope"' in moon_scene
    assert html.count('class="sysmap-facility-orbit"') == 4
    assert "sysmap-legend-facility" in html

    plain = render_system_map_panel(system, system["stars"], system["planets"], system["belts"])
    assert "sysmap-facility" not in plain and "sysmap-legend-facility" not in plain


# --- Sector page ----------------------------------------------------------------------

def test_sector_contents_list_facilities_outside_systems(app):
    from planetgen.web.sector_page import _contents

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
            system_id = store.save_system(system, cfg, config=mysql_config)
            break
    else:
        pytest.fail("could not generate a system with a terrestrial planet")
    conn = store.get_connection(mysql_config)
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
    conn = store.get_connection(mysql_config)
    try:
        detail = query.system_detail(conn, system_id)
        before = next(p for p in detail["planets"] if p["id"] == ids["terrestrial"])
        with conn:
            store.add_facility(conn, "New Hope", "colony", "terrestrial", "planet", ids["terrestrial"])
        detail = query.system_detail(conn, system_id)
    finally:
        conn.close()
    planet = next(p for p in detail["planets"] if p["id"] == ids["terrestrial"])
    assert planet["inhabited"] is True
    assert before["habitable"] == planet["habitable"]  # a colony changes nothing else
    others = [p for p in detail["planets"] if p["id"] != ids["terrestrial"]]
    assert all(p["inhabited"] == (p["life_stage"] == "technological_civilization") for p in others)

    client = db_app.test_client()
    html = client.get(f"/system/{pid('system', system_id)}").get_data(as_text=True)
    assert "New Hope" in _panel(html)
    # UX.176: the Inhabited chip and the facilities a body hosts are on the body's own page.
    page = client.get(f"/object/planet:{pid('planet', planet['id'])}").get_data(as_text=True)
    assert '<span class="flag flag-yes">Inhabited</span>' in page and "New Hope" in page
    # UX.224: Inhabited says the world is habitable and of a habitable class, so those two chips are not shown beside it.
    assert '<span class="flag flag-yes">Habitable</span>' not in page and '<span class="flag flag-yes">Habitable class</span>' not in page


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
    conn = store.get_connection(mysql_config)
    try:
        star = conn.execute("SELECT radius_km, heliosphere_radius_km FROM stars WHERE id = ?",
                            (ids["star"],)).fetchone()
    finally:
        conn.close()
    limits = facility_rules.orbit_limits(star["radius_km"], star["heliosphere_radius_km"])
    form = {"name": "Sunwatch", "host": f"star:{pid('star', ids['star'])}", "placement": "orbital", "kind": "outpost",
            "orbit_step": str(facility_rules.step_for_distance(0.5 * AU_KM, *limits))}

    resp = client.post(f"/system/{pid('system', system_id)}", data={csrf.FIELD_NAME: _csrf(db_app, client),
                                                     "facility_action": "preview", **form})
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200
    assert "Allowed by the placement rules." in html
    # Half an AU around a Sun-like star: about 0.35 years, about 42 km/s.
    speed = float(re.search(r"at ([\d.,]+) km/s", _panel(html)).group(1).replace(",", ""))
    assert 35 < speed < 50
    conn = store.get_connection(mysql_config)
    try:
        assert query.facilities_for_system(conn, system_id) == []
    finally:
        conn.close()

    resp = client.post(f"/system/{pid('system', system_id)}", data={csrf.FIELD_NAME: _csrf(db_app, client),
                                                     "facility_action": "save", **form})
    assert resp.status_code == 303
    conn = store.get_connection(mysql_config)
    try:
        (saved,) = query.facilities_for_system(conn, system_id)
    finally:
        conn.close()
    assert saved["orbit_distance_km"] == pytest.approx(0.5 * AU_KM, rel=0.01)
    assert saved["orbital_speed_kms"] == pytest.approx(speed, rel=0.01)
    page = client.get(resp.headers["Location"]).get_data(as_text=True)
    assert "Facility added." in page and "Sunwatch" in _panel(page)
    assert 'data-kind="facility"' in page

    resp = client.post(f"/system/{pid('system', system_id)}", data={csrf.FIELD_NAME: _csrf(db_app, client),
                                                     "facility_action": "remove", "facility_id": pid("facility", saved["id"])})
    assert resp.status_code == 303
    assert "Facility removed." in client.get(resp.headers["Location"]).get_data(as_text=True)
    conn = store.get_connection(mysql_config)
    try:
        assert query.facilities_for_system(conn, system_id) == []
    finally:
        conn.close()


def test_real_sector_page_lists_stand_alone_facilities(db_app, mysql_config):
    sector = SpaceSector("Harbor", edge_ly=13.0)
    sector_id = store.save_sector(sector, config=mysql_config, galaxy_position={
        "center_x_pc": 0.0, "center_y_pc": 50.0, "center_z_pc": 0.0, "galactic_radius_pc": 50.0,
    })
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            store.add_facility(conn, "Waypoint", "station", "standalone", "space", sector_id, offset_ly=(1.0, 2.0, 2.0))
    finally:
        conn.close()
    html = db_app.test_client().get(f"/sector/{pid('sector', sector_id)}").get_data(as_text=True)
    row = re.search(r"<tr[^>]*>\s*<td>Waypoint</td>.*?</tr>", html, re.S).group(0)
    assert "Facility (Station)" in row
    assert "Stand-alone station, parked in open space" in row
    assert "3.00 ly" in row or "3 ly" in row, row


# --- Type, affiliation and ratings (POP.10) --------------------------------------------

def _types():
    return [
        {"id": 4, "name": "Research outpost", "kind": "outpost", "icon": "flask", "active": True, "fields": [
            {"id": 1, "key": "staff", "label": "Staff", "datatype": "int", "unit": None, "required": True,
             "min_value": 1, "max_value": 100, "choices": None, "default_text": "12"},
            {"id": 2, "key": "focus", "label": "Focus", "datatype": "choice", "unit": None, "required": False,
             "min_value": None, "max_value": None, "choices": ["Survey", "Relay"], "default_text": None},
            {"id": 3, "key": "armed", "label": "Armed", "datatype": "bool", "unit": None, "required": False,
             "min_value": None, "max_value": None, "choices": None, "default_text": None},
            {"id": 5, "key": "output", "label": "Output", "datatype": "decimal", "unit": "TW", "required": False,
             "min_value": None, "max_value": None, "choices": None, "default_text": None},
        ]},
        {"id": 6, "name": "Farm", "kind": "colony", "icon": None, "active": True, "fields": []},
    ]


def test_the_form_offers_types_with_their_fields_affiliation_and_ratings(client, fake):
    fake.types = _types()
    fake.owner = {"polity_id": 3, "polity_name": "Velarix Concord", "color": "#3366cc", "distance_ly": 2.0}
    _as_admin(client, fake)
    dialog = _panel(client.get("/system/5").get_data(as_text=True))
    assert '<sl-option value="4" data-kind="outpost">Research outpost</sl-option>' in dialog
    assert '<sl-option value="6" data-kind="colony">Farm</sl-option>' in dialog
    assert 'data-facility-type-fields="4" hidden' in dialog and 'data-facility-type-fields="6" hidden' in dialog
    assert 'name="field_staff"' in dialog and 'min="1"' in dialog and 'value="12"' in dialog   # the default
    assert '<sl-option value="Survey">Survey</sl-option>' in dialog and 'name="field_armed"' in dialog and "Output (TW)" in dialog
    assert '<sl-option value="polity:3">Velarix Concord (owns this system)</sl-option>' in dialog
    assert re.search(r'<sl-select name="affiliation"[^>]*value="polity:3"', dialog)   # the default
    for dimension in ("crime", "housing", "resources", "maintenance", "health"):
        assert f'name="rating_{dimension}"' in dialog
    assert '<sl-option value="red">Red</sl-option>' in dialog


def test_save_sends_the_type_fields_affiliation_and_ratings(app, client, fake):
    fake.types = _types()
    _as_admin(client, fake)
    response = _post(app, client, **_STAR_OUTPOST, facility_action="save", type_id="4", field_staff="9",
                     field_focus="Relay", field_armed="no", field_output="1.5", affiliation="polity:3",
                     rating_crime="red", rating_health="green", ratings_note="Survey 2200")
    assert response.status_code == 303 and "facility=added" in response.headers["Location"]
    [(_, _, body)] = [call for call in fake.calls if call[0] == "create"]
    assert body["type_id"] == 4
    assert body["fields"] == {"staff": 9, "focus": "Relay", "armed": False, "output": 1.5}
    assert body["affiliation_polity_id"] == 3 and "affiliation_name" not in body
    assert body["ratings"] == {"crime": "red", "health": "green"} and body["ratings_note"] == "Survey 2200"


def test_free_text_affiliation_wins_and_an_untyped_save_sends_no_fields(app, client, fake):
    fake.types = _types()
    _as_admin(client, fake)
    _post(app, client, **_STAR_OUTPOST, facility_action="save", affiliation="polity:3",
          affiliation_name="Free Traders", field_staff="9")
    [(_, _, body)] = [call for call in fake.calls if call[0] == "create"]
    assert body["affiliation_name"] == "Free Traders" and "affiliation_polity_id" not in body
    assert "type_id" not in body and "fields" not in body and "ratings" not in body


@pytest.mark.parametrize("form, message", [
    ({"type_id": "99"}, "Choose one of the listed facility types."),
    ({"type_id": "4", "field_staff": "many"}, "Staff must be a whole number."),
    ({"type_id": "4", "field_output": "lots"}, "Output must be a number."),
    ({"type_id": "4", "field_armed": "maybe"}, "Armed must be yes or no."),
    ({"rating_crime": "amber"}, "A rating is Green, Yellow or Red."),
])
def test_a_bad_type_field_or_rating_is_reported_before_the_api_is_asked(app, client, fake, form, message):
    fake.types = _types()
    _as_admin(client, fake)
    response = _post(app, client, **_STAR_OUTPOST, facility_action="save", **form)
    assert response.status_code == 200 and message in response.get_data(as_text=True)
    assert not [call for call in fake.calls if call[0] == "create"]


def test_a_facility_can_be_edited_and_have_its_ratings_suggested(app, client, fake):
    fake.facilities[0].update(affiliation_polity_id=3, affiliation_name="Velarix Concord", overall="red",
                              ratings={"crime": "red", "housing": None, "resources": None, "maintenance": None,
                                       "health": None},
                              rating_words={"crime": "High", "housing": None, "resources": None,
                                            "maintenance": None, "health": None})
    fake.owner = {"polity_id": 3, "polity_name": "Velarix Concord", "color": "#3366cc", "distance_ly": 2.0}
    _as_admin(client, fake)
    html = client.get("/system/5").get_data(as_text=True)
    dialog = re.search(r'<sl-dialog id="admin-facility-edit-edit-1".*?</sl-dialog>', html, re.S).group(0)
    assert '<option value="red" selected>Red</option>' in dialog         # its own crime rating
    assert '<option value="polity:3" selected>Velarix Concord (owns this system)</option>' in dialog

    response = _post(app, client, facility_action="edit", facility_id="1", affiliation="", affiliation_name="Nobody",
                     rating_crime="green", rating_health="", ratings_note="")
    assert response.status_code == 303 and "facility=edited" in response.headers["Location"]
    [(_, _, facility_id, body)] = [call for call in fake.calls if call[0] == "update"]
    assert facility_id == "1" and body["affiliation_name"] == "Nobody" and body["affiliation_polity_id"] is None
    assert body["ratings"] == {"crime": "green", "housing": None, "resources": None, "maintenance": None,
                               "health": None} and body["ratings_note"] is None

    response = _post(app, client, facility_action="suggest", facility_id="1")
    assert response.status_code == 303 and "facility=suggested" in response.headers["Location"]
    fake.suggested = []
    assert "facility=nothing" in _post(app, client, facility_action="suggest", facility_id="1").headers["Location"]
    assert "not in this system" in _post(app, client, facility_action="edit",
                                         facility_id="77").get_data(as_text=True)
    assert [call[0] for call in fake.calls].count("update") == 1            # the unknown facility changed nothing
