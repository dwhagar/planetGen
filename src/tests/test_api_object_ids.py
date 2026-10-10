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
from tests.publicids import pid

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


_ROW_LINK = re.compile(r'(?:/(?:system|sector|object)/|/phenomenon/[a-z_]+/|\b(?:system|sector|star|planet|moon|belt|comet|'
                       r'facility|nebula|black_hole|rogue_planet):)(\d+)\b')


def test_pages_link_to_objects_by_id_and_never_by_row_number(bodies):
    web, _config = bodies
    sector = web.get("/api/sectors").get_json()["items"][0]["id"]
    system = web.get("/api/systems").get_json()["items"][0]["id"]
    paths = ["/sectors", "/systems", "/phenomena", "/nav", "/nearby", f"/sector/{sector}", f"/system/{system}",
             f"/system/{system}/scene", f"/nav?from=system:{system}", f"/nearby?place=system:{system}&distance=10",
             f"/galaxy?course=system:{system}"]
    for path in paths:
        response = web.get(path)
        assert response.status_code in (200, 302, 404), path
        found = [m.group(0) for m in _ROW_LINK.finditer(response.get_data(as_text=True)) if len(m.group(1)) < 9]
        assert found == [], (path, found[:5])


def test_the_galaxy_stage_and_tiles_carry_sector_ids_the_map_can_open(mysql_config):
    """The Galaxy Map opens `/sector/<id>/scene` with the ids from the stage payload and the tiles."""
    from planetgen.galaxy.viewport import tiles_intersecting_sphere
    from tests.test_api import _generated_sectors

    _generated_sectors(mysql_config, [0.5, 5.0, 40.0])

    class TestConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config

    app = create_app(TestConfig)
    app.testing = True
    client = app.test_client()
    keys = tiles_intersecting_sphere(9, (60.0, 60.0, 60.0), 130.0)
    tiles = client.get("/api/galaxy/tiles?tiles=" + ",".join(keys)).get_json()["tiles"]
    placed = [row for tile in tiles.values() for row in tile["placed"]]
    assert placed and all(isinstance(row["id"], str) for row in placed)
    assert client.get(f"/sector/{placed[0]['id']}/scene").status_code == 200

    at = None
    for _ in range(6):
        stage = client.get("/api/galaxy/stage" + (f"?at={at}" if at else "")).get_json()
        if stage.get("sectors"):
            break
        child = next((c for c in stage["children"] if c.get("generated")), None)
        if child is None:
            break
        at = f"{stage['child_m']}.{child['ring']}.{child['wedge']}.{child['slab']}"
    assert stage.get("sectors"), "no stage with sectors reached"
    assert all(isinstance(row["id"], str) for row in stage["sectors"])
    assert client.get(f"/sector/{stage['sectors'][0]['id']}/scene").status_code == 200


def test_a_rogue_planet_inside_a_nebula_names_the_nebula_by_id(mysql_config):
    from planetgen.generation.phenomena.nebula import Nebula

    sector = SpaceSector("Cloud Sector", edge_ly=20.0)
    sector.add_phenomenon(RoguePlanet(SystemConfig()), "rogue-planet")
    sector_row = store.save_sector(sector, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            nebula = store.insert_nebula(conn, Nebula(SystemConfig(), name="Test Veil"), sector_id=sector_row,
                                         placement={"center_x_pc": 1.0, "center_y_pc": 2.0, "center_z_pc": 3.0,
                                                    "galactic_radius_pc": 3.7})
            store.assign_uids(conn, phenomenon=("nebulae", nebula))
            conn.execute("UPDATE rogue_planets SET inside_nebula_id = ?", (nebula,))
        nebula_printed = pid("nebula", nebula, conn)
    finally:
        conn.close()

    class TestConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config

    app = create_app(TestConfig)
    app.testing = True
    client = app.test_client()
    rogue = next(row for row in client.get("/api/phenomena").get_json()["items"] if row["type"] == "rogue_planet")
    detail = client.get(f"/api/phenomena/rogue_planet/{rogue['id']}").get_json()
    assert detail["inside_nebula_id"] == nebula_printed
    page = client.get(f"/phenomenon/rogue_planet/{rogue['id']}")
    assert page.status_code == 200
    assert f"/phenomenon/nebula/{nebula_printed}" in page.get_data(as_text=True)
