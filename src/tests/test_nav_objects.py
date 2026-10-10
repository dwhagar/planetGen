# tests/test_nav_objects.py

"""
NAV.16: `/api/nav` and `queryDb.nav_course` take any object reference; a
course to or from a body has legs inside its system.
"""

import math

import pytest

from planetgen.db import store
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem

# Imported fixtures (see test_bughunt_api_gaps.py for why this works).
from tests.test_api import (  # noqa: F401
    _place_sector, admin_client, client, default_admin_client, first_admin_password,
)
from tests.publicids import pid, pids


def _planet_sector(mysql_config, name, center_pc):
    """Saves a placed sector holding one single-star system with planets."""
    for _ in range(60):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.PLANETS = True
        cfg.MAX_PLANETS = True
        cfg.BINARY_SYSTEM = False
        system = StarSystem(system_config=cfg)
        if any(p.body_type != "a" for p in system.planets):
            break
    sector = SpaceSector(name, edge_ly=10.0)
    sector.add_system(system, position=(0.0, 0.0, 0.0), system_config=cfg)
    sector_id = store.save_sector(sector, config=mysql_config, galaxy_position={
        "center_x_pc": center_pc[0], "center_y_pc": center_pc[1], "center_z_pc": center_pc[2],
        "galactic_radius_pc": math.dist(center_pc, (0.0, 0.0, 0.0)),
    })
    conn = store.get_connection(mysql_config)
    try:
        system_id = conn.execute("SELECT id FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()["id"]
        planets = [r["id"] for r in conn.execute(
            "SELECT id FROM planets WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchall()]
    finally:
        conn.close()
    return system_id, planets


def test_a_course_inside_one_system_is_one_leg(client, mysql_config):
    system_id, planets = _planet_sector(mysql_config, "Home", (10.0, 0.0, 0.0))
    body = client.get(f"/api/nav?from=planet:{pid('planet', planets[0])}&to=system:{pid('system', system_id)}").get_json()
    assert body["scope"] == "system" and body["route"] is None
    assert [leg["kind"] for leg in body["legs"]] == ["within"]
    assert body["direct"]["frame"] == "system" and body["direct"]["distance_ly"] > 0
    assert body["origin"]["ref"] == f"planet:{pid('planet', planets[0])}" and body["destination"]["ref"] == f"system:{pid('system', system_id)}"
    assert "same warp and fold tables" in body["note"]
    assert len(body["warp_times"]) > 0


def test_a_course_between_bodies_in_two_systems_has_three_legs(client, mysql_config):
    a_system, a_planets = _planet_sector(mysql_config, "Near", (10.0, 0.0, 0.0))
    b_system, b_planets = _planet_sector(mysql_config, "Far", (40.0, 0.0, 0.0))
    body = client.get(f"/api/nav?from=planet:{pid('planet', a_planets[0])}&to=planet:{pid('planet', b_planets[0])}").get_json()
    assert body["scope"] == "galaxy"
    assert [leg["kind"] for leg in body["legs"]] == ["out", "between", "into"]
    out, between, into = body["legs"]
    assert out["from"] == f"planet:{pid('planet', a_planets[0])}" and out["to"] == f"system:{pid('system', a_system)}"
    assert between["from"] == f"system:{pid('system', a_system)}" and between["to"] == f"system:{pid('system', b_system)}"
    assert into["to"] == f"planet:{pid('planet', b_planets[0])}"
    assert between["direct"]["distance_ly"] == body["direct"]["distance_ly"]
    assert out["direct"]["frame"] == into["direct"]["frame"] == "system"
    # The in-system legs are a heliopause long -- far shorter than the jump between systems.
    assert 0 < out["direct"]["distance_ly"] < between["direct"]["distance_ly"] / 100
    assert body["route"]["path"][0] == pid("system", a_system) and body["route"]["path"][-1] == pid("system", b_system)


def test_system_to_system_is_unchanged_with_one_leg(client, mysql_config):
    a_system, _ = _planet_sector(mysql_config, "Near", (10.0, 0.0, 0.0))
    b_system, _ = _planet_sector(mysql_config, "Far", (40.0, 0.0, 0.0))
    body = client.get(f"/api/nav?from={pid('system', a_system)}&to=system:{pid('system', b_system)}").get_json()
    assert [leg["kind"] for leg in body["legs"]] == ["between"] and body["note"] is None


def test_nav_refuses_a_sector_and_bad_references(client, mysql_config):
    a_system, _ = _planet_sector(mysql_config, "Near", (10.0, 0.0, 0.0))
    for params, status in (
        (f"from=sector:{pid('sector', 1)}&to=system:{pid('system', a_system)}", 400),
        (f"from=ship:1&to=system:{pid('system', a_system)}", 400),
        (f"from=moon:FFFFFFFFFF-FFFFFFF-FFF&to=system:{pid('system', a_system)}", 404),
        (f"from={pid('system', a_system)}", 400),
    ):
        response = client.get(f"/api/nav?{params}")
        assert response.status_code == status, (params, response.get_json())
