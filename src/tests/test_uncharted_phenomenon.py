# tests/test_uncharted_phenomenon.py

"""
MAP.211: a phenomenon the scatter placed has a page before its sector is generated (what is known of it, with the
kind of star a hypervelocity star will be), and once the sector is generated the page sends the visitor to the
object that was made. Runs against a real throwaway database through the in-process transport.
"""

import math

import pytest

from planetgen.api.config import Config
from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.compact_remnant import BlackHole
from planetgen.generation.system import StarSystem
from planetgen.web.app import create_app
from planetgen.web.lib import apiclient
from planetgen.web.maps.starmap import uncharted_scene_data
from tests.publicids import pid
from tests.test_api import client  # noqa: F401

MPC_PER_PC = 1000
SPEED_KMS = 640.0
CELL = (7, 2, 3)


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


def _scatter(mysql_config, kind, subtype=None, point_pc=(30.0, 40.0, 0.0), velocity=(None, None, None), seed=1,
             mass=None, cell=CELL, record=True):
    """Inserts one `phenomenon_scatter` row and returns its id."""
    conn = store.get_connection(mysql_config)
    row = (*cell, kind, subtype, *(round(value * MPC_PER_PC) for value in point_pc), *velocity, seed, mass, None)
    store.insert_phenomenon_scatter(conn, [row])
    if record:
        store.record_phenomenon_scatter_classes(conn, {(kind, subtype or ""): 1})
    scatter_id = conn.execute("SELECT MAX(id) AS id FROM phenomenon_scatter").fetchone()["id"]
    conn.commit()
    conn.close()
    return scatter_id


def _build(mysql_config, scatter_id):
    conn = store.get_connection(mysql_config)
    store.mark_phenomena_built(conn, [scatter_id])
    conn.commit()
    conn.close()


def _sector(mysql_config, cell=CELL):
    conn = store.get_connection(mysql_config)
    conn.execute("INSERT INTO sectors (name, edge_mpc, shell_index, band_index, wedge_index) VALUES (?, ?, ?, ?, ?)",
                 ("Cell Seven", 4000.0, *cell))
    sector_id = conn.execute("SELECT MAX(id) AS id FROM sectors").fetchone()["id"]
    conn.commit()
    conn.close()
    return sector_id


def _hypervelocity_system(mysql_config, sector_id, speed=SPEED_KMS):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "B5V"
    cfg.PLANETS = False
    cfg.BINARY_SYSTEM = False
    system = StarSystem(system_config=cfg)
    system.runaway_class = "hypervelocity"
    system.runaway_speed_kms = speed
    system.runaway_direction = (1.0, 0.0, 0.0)
    system_id = store.save_system(system, cfg, config=mysql_config)
    conn = store.get_connection(mysql_config)
    conn.execute("UPDATE star_systems SET sector_id = ? WHERE id = ?", (sector_id, system_id))
    conn.commit()
    conn.close()
    return system_id


def test_an_uncharted_hypervelocity_star_has_a_page_with_its_kind_of_star(db_app, mysql_config):
    scatter_id = _scatter(mysql_config, "hypervelocity-star", velocity=(SPEED_KMS, 0.0, 0.0))
    response = db_app.test_client().get(f"/uncharted-phenomenon/{scatter_id}")
    html = response.get_data(as_text=True)
    assert response.status_code == 200
    assert "Hypervelocity Star (uncharted)" in html
    assert "B-type, a young main-sequence star" in html           # the kind of star it will be
    assert "Ejected from the central black hole" in html
    assert "640" in html and "Speed" in html and "Time in flight" in html
    assert "Uncharted" in html


def test_an_uncharted_black_hole_page_names_its_class(db_app, mysql_config):
    scatter_id = _scatter(mysql_config, "black-hole", "stellar", mass=12.5)
    response = db_app.test_client().get(f"/uncharted-phenomenon/{scatter_id}")
    html = response.get_data(as_text=True)
    # TEST.284: it failed twice in a loaded gate with the class name missing and nothing to say why: name the answer.
    assert "Black Hole (Stellar)" in html and "Mass" in html, (response.status_code, dict(response.headers), html[:1500])
    assert "Expected star" not in html                              # only a hypervelocity star will be a star


def test_an_unknown_scatter_row_is_a_404(db_app):
    assert db_app.test_client().get("/uncharted-phenomenon/999999").status_code == 404


def test_a_built_hypervelocity_star_goes_to_its_star_system(db_app, mysql_config):
    scatter_id = _scatter(mysql_config, "hypervelocity-star", velocity=(SPEED_KMS, 0.0, 0.0))
    sector_id = _sector(mysql_config)
    # A second hypervelocity system in the cell flies at another speed: the row's own speed picks the right one.
    _hypervelocity_system(mysql_config, sector_id, speed=SPEED_KMS + 17.0)
    system_id = _hypervelocity_system(mysql_config, sector_id)
    system_pid = pid("system", system_id)
    _build(mysql_config, scatter_id)
    response = db_app.test_client().get(f"/uncharted-phenomenon/{scatter_id}")
    assert response.status_code == 302
    assert response.headers["Location"].endswith(f"/system/{system_pid}")


def test_a_built_black_hole_goes_to_its_phenomenon_page(db_app, mysql_config):
    point = (30.0, 40.0, 0.0)
    scatter_id = _scatter(mysql_config, "black-hole", "stellar", point_pc=point, mass=12.5)
    sector_id = _sector(mysql_config)
    cfg = SystemConfig()
    near = store.save_phenomenon(BlackHole(cfg), cfg, "black-hole", config=mysql_config)
    far = store.save_phenomenon(BlackHole(cfg), cfg, "black-hole", config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        for row_id, offset in ((near, 0.001), (far, 3.0)):
            conn.execute("UPDATE black_holes SET sector_id = ?, center_x_pc = ?, center_y_pc = ?, center_z_pc = ?,"
                         " galactic_radius_pc = ? WHERE id = ?",
                         (sector_id, point[0] + offset, point[1], point[2], math.dist((0, 0, 0), point), row_id))
        conn.commit()
    finally:
        conn.close()
    near_id = pid("black_hole", near)
    _build(mysql_config, scatter_id)
    response = db_app.test_client().get(f"/uncharted-phenomenon/{scatter_id}")
    assert response.status_code == 302
    assert response.headers["Location"].endswith(f"/phenomenon/black_hole/{near_id}")


def test_a_built_row_that_cannot_be_matched_says_so(db_app, mysql_config):
    scatter_id = _scatter(mysql_config, "neutron-star")
    _sector(mysql_config)
    _build(mysql_config, scatter_id)
    response = db_app.test_client().get(f"/uncharted-phenomenon/{scatter_id}")
    assert response.status_code == 200
    assert "could not be matched" in response.get_data(as_text=True)


def test_the_phenomena_list_links_each_uncharted_row_to_its_page(db_app, mysql_config):
    first = _scatter(mysql_config, "hypervelocity-star", velocity=(SPEED_KMS, 0.0, 0.0), cell=(7, 2, 3), record=False)
    second = _scatter(mysql_config, "neutron-star", cell=(7, 2, 4), record=False)
    conn = store.get_connection(mysql_config)
    store.record_phenomenon_scatter_classes(conn, {("hypervelocity-star", ""): 1, ("neutron-star", ""): 1})
    conn.commit()
    conn.close()
    html = db_app.test_client().get("/phenomena").get_data(as_text=True)
    assert f'href="/uncharted-phenomenon/{first}"' in html
    assert f'href="/uncharted-phenomenon/{second}"' in html


def test_the_api_describes_the_row(client, mysql_config):
    scatter_id = _scatter(mysql_config, "hypervelocity-star", point_pc=(300.0, 400.0, 0.0),
                          velocity=(SPEED_KMS, 0.0, 0.0))
    response = client.get(f"/api/uncharted-phenomena/{scatter_id}")
    assert response.status_code == 200
    detail = response.get_json()
    assert detail["type"] == "hypervelocity_star" and not detail["built"]
    assert detail["speed_kms"] == pytest.approx(SPEED_KMS)
    assert detail["direction"] == pytest.approx([1.0, 0.0, 0.0])
    assert detail["distance_from_core_pc"] == pytest.approx(500.0)
    years = 500.0 * 3.085677581491367e16 / (SPEED_KMS * 1000.0) / 31_557_600.0
    assert detail["flight_years"] == pytest.approx(years)
    assert detail["address"] == list(CELL) and "target" not in detail
    assert client.get("/api/uncharted-phenomena/999999").status_code == 404


def test_the_uncharted_sector_scene_carries_the_link_and_the_speed():
    contents = {"stars": [], "scattered": [
        {"id": 5, "kind": "hypervelocity-star", "subtype": None, "type": "hypervelocity_star",
         "x": 101.0, "y": 50.0, "z": 0.0, "speed_kms": 640.0},
        {"id": 6, "kind": "black-hole", "subtype": "stellar", "type": "black_hole",
         "x": 100.5, "y": 50.0, "z": 0.0, "speed_kms": None}]}
    scene = uncharted_scene_data((3, 0, 1), "ABC123", (100.0, 50.0, 0.0), 4.0, contents,
                                 link_url=lambda name, **params: f"/{name}/{params['scatter_id']}")
    star = scene["stars"][0]
    assert star["typeLabel"] == "Hypervelocity star" and star["href"] == "/uncharted_phenomenon/5"
    assert star["speedText"] and "km" in star["speedText"] or "m/s" in star["speedText"]
    assert scene["clouds"][0]["href"] == "/uncharted_phenomenon/6"
    assert math.isfinite(star["x"])
