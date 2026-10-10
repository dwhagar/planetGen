# tests/test_api_object_ids.py

"""
API.23, stage 1: sectors, systems and phenomena speak object IDs in the API, never row numbers.
"""

import re

import pytest

from planetgen.api.config import Config
from planetgen.db import store
from planetgen.galaxy import uid as sector_uid
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.web.app import create_app

PRINTED = re.compile(r"^[0-9A-F]+(-[0-9A-F]+){0,2}$")


@pytest.fixture
def api(mysql_config):
    sector = SpaceSector("Id Sector", edge_ly=20.0)
    cfg = SystemConfig()
    cfg.PLANETS = False
    cfg.BINARY_SYSTEM = False
    sector.add_system(StarSystem(system_config=cfg), position=(1.0, 1.0, 1.0), system_config=cfg)
    sector.add_system(StarSystem(system_config=cfg), position=(-2.0, 0.5, 3.0), system_config=cfg)
    sector.add_phenomenon(RoguePlanet(SystemConfig()), "rogue-planet")
    store.save_sector(sector, config=mysql_config)

    class TestConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config

    app = create_app(TestConfig)
    app.testing = True
    return app.test_client()


def _printed(value):
    """An object ID as printed: three hyphen-joined hex groups."""
    return isinstance(value, str) and "-" in value and bool(PRINTED.match(value))


def test_lists_and_details_carry_printed_ids(api):
    sectors = api.get("/api/sectors").get_json()["items"]
    sector_id = sectors[0]["id"]
    assert all(isinstance(row["id"], str) for row in sectors)
    sector = api.get(f"/api/sectors/{sector_id}").get_json()
    assert sector["id"] == sector_id
    systems = sector["systems"]
    assert len(systems) == 2 and all(_printed(row["id"]) for row in systems)
    system = api.get(f"/api/systems/{systems[0]['id']}").get_json()
    assert system["id"] == systems[0]["id"] and _printed(system["id"])
    assert all(_printed(row["id"]) for row in system["nearest_neighbors"])
    phenomena = api.get("/api/phenomena").get_json()["items"]
    assert phenomena and all(_printed(row["id"]) for row in phenomena)
    one = api.get(f"/api/phenomena/{phenomena[0]['type']}/{phenomena[0]['id']}").get_json()
    assert one["id"] == phenomena[0]["id"]


def test_a_row_number_names_nothing(api):
    for path in ("/api/systems/1", "/api/systems/2/scene", "/api/phenomena/rogue_planet/1"):
        assert api.get(path).status_code == 404, path


@pytest.fixture
def bodies(mysql_config):
    """A client over a sector holding a wide binary with planets, moons and a comet, one facility on a planet,
    and one hand-made sector with no grid address and no stored ID."""
    sector = SpaceSector("Body Sector", edge_ly=11.5)
    for _ in range(400):
        cfg = SystemConfig()
        cfg.PLANETS = True
        cfg.MAX_PLANETS = True
        cfg.BINARY_SYSTEM = True
        system = StarSystem(system_config=cfg)
        if len(system.stars) > 1 and system.comets and any(getattr(p, "moons", None) for p in system.planets):
            break
    else:
        pytest.skip("no suitable system generated")
    sector.add_system(system, position=(1.0, 1.0, 1.0), system_config=cfg)
    store.save_sector(sector, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        planet = conn.execute("SELECT id, star_system_id FROM planets ORDER BY id LIMIT 1").fetchone()
        store.add_facility(conn, "Relay", "station", "orbital", "planet", planet["id"], distance_km=None)
        conn.execute("INSERT INTO sectors (name, edge_mpc) VALUES (?, ?)", ("By Hand", 11500))
        conn.commit()
    finally:
        conn.close()

    class TestConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config

    app = create_app(TestConfig)
    app.testing = True
    return app.test_client(), mysql_config


def _row_numbers(node, path=""):
    """Every `id`-like key under `node` holding an integer: a row number that leaked."""
    found = []
    if isinstance(node, dict):
        for key, value in node.items():
            here = f"{path}.{key}"
            if (key == "id" or key.endswith("_id")) and isinstance(value, int) and not isinstance(value, bool):
                found.append(here)
            found += _row_numbers(value, here)
    elif isinstance(node, list):
        for index, value in enumerate(node):
            found += _row_numbers(value, f"{path}[{index}]")
    return found


def test_system_answers_leak_no_row_numbers(bodies):
    api, _config = bodies
    system_id = api.get("/api/systems").get_json()["items"][0]["id"]
    for suffix in ("", "/scene", "/facilities", "/sections"):
        body = api.get(f"/api/systems/{system_id}{suffix}").get_json()
        assert _row_numbers(body) == [], suffix
    sections = api.get(f"/api/systems/{system_id}/sections").get_json()
    for kind in ("stars", "planets", "moons", "belts", "comets"):
        assert all(_printed(key) for key in sections.get(kind, {})), kind
    scene = api.get(f"/api/systems/{system_id}/scene").get_json()
    for planet in scene["planets"]:
        assert _printed(planet["id"]) and planet["ref"] == f"planet:{planet['id']}"
        for moon in planet["moons"]:
            assert moon["parent"] == planet["ref"]


def test_a_body_is_reached_by_its_id_and_a_row_number_names_nothing(bodies):
    api, _config = bodies
    system_id = api.get("/api/systems").get_json()["items"][0]["id"]
    detail = api.get(f"/api/systems/{system_id}").get_json()
    planet = detail["planets"][0]["id"]
    assert api.get(f"/api/objects/planet:{planet}").status_code == 200
    assert api.get("/api/objects/planet:1").status_code in (400, 404)
    facility = api.get(f"/api/systems/{system_id}/facilities").get_json()
    facilities = facility.get("items") or facility.get("facilities") or facility
    assert api.get(f"/api/facilities/{facilities[0]['id']}").get_json()["host_id"] == planet
    assert api.get("/api/facilities/1").status_code == 404


def test_a_hand_made_sector_with_no_stored_id_is_found_by_its_unplaced_id(bodies):
    api, config = bodies
    conn = store.get_connection(config)
    try:
        row = conn.execute("SELECT id FROM sectors WHERE name = 'By Hand'").fetchone()
    finally:
        conn.close()
    printed = format(sector_uid.unplaced_sector_uid(row["id"]), "X")
    found = api.get(f"/api/sectors/{printed}")
    assert found.status_code == 200 and found.get_json()["name"] == "By Hand"
