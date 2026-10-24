# tests/test_facility_pages.py

"""
UX.183: a facility is an object a reference can name: `facility:<ID>` resolves with its place, its property sheet is
served, and its page (`/object/facility:<ID>`) shows where it stands and what its record says.
"""

import math

import pytest

from planetgen.api.config import Config
from planetgen.db import store
from planetgen.galaxy import objectref
from planetgen.galaxy.sector import SpaceSector
from planetgen.galaxy.spherical import sector_position_pc
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.web import item_pages
from planetgen.web.app import create_app


@pytest.fixture
def galaxy(mysql_config):
    """An API client over one sector with a system, a facility on its first planet and a stand-alone one parked in the
    sector, and the two facilities' printed IDs."""
    center = sector_position_pc(10, 0, 3, 4.0)
    sector = SpaceSector("Dock Sector", edge_ly=13.0)
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = False
    sector.add_system(StarSystem(system_config=cfg), position=(1.0, 1.0, 1.0), system_config=cfg)
    store.save_sector(sector, config=mysql_config, galaxy_position={
        "center_x_pc": center[0], "center_y_pc": center[1], "center_z_pc": center[2],
        "galactic_radius_pc": math.hypot(*center), "shell_index": 10, "band_index": 0, "wedge_index": 3})
    conn = store.get_connection(mysql_config)
    try:
        planet = conn.execute("SELECT id FROM planets WHERE body_type = 't' ORDER BY id LIMIT 1").fetchone()
        if planet is None:
            pytest.skip("no rocky planet generated")
        sector_id = conn.execute("SELECT id FROM sectors ORDER BY id LIMIT 1").fetchone()["id"]
        store.add_facility(conn, "New Hope", "colony", "terrestrial", "planet", planet["id"], description="A first colony.")
        store.add_facility(conn, "Lone Light", "station", "standalone", "space", sector_id)
        conn.commit()
    finally:
        conn.close()

    class TestConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config

    app = create_app(TestConfig)
    app.testing = True
    client = app.test_client()
    system_id = client.get("/api/sectors").get_json()["items"][0]["id"]
    system = client.get(f"/api/sectors/{system_id}").get_json()["systems"][0]["id"]
    ids = {f["name"]: f["id"] for f in client.get(f"/api/systems/{system}/facilities").get_json()["items"]}
    ids.update({f["name"]: f["id"] for f in client.get(f"/api/sectors/{system_id}/facilities").get_json()["items"]})
    return client, ids


def test_facility_is_a_kind_a_reference_can_name():
    assert "facility" in objectref.KINDS and objectref.TABLES["facility"] == "facilities"
    assert objectref.parse_public("facility:FE81000A2B-0000005-007") == ("facility", "FE81000A2B-0000005-007")
    assert "facility" in item_pages.page_kinds()


def test_a_facility_on_a_planet_resolves_with_the_chain_down_to_its_host(galaxy):
    client, ids = galaxy
    found = client.get(f"/api/objects/facility:{ids['New Hope']}").get_json()
    assert found["kind"] == "facility" and found["name"] == "New Hope"
    assert [parent["kind"] for parent in found["parents"]] == ["galaxy", "sector", "system", "planet"]
    assert found["positions"]["galaxy_pc"] is not None


def test_a_standalone_facility_stands_in_its_sector(galaxy):
    client, ids = galaxy
    found = client.get(f"/api/objects/facility:{ids['Lone Light']}").get_json()
    assert [parent["kind"] for parent in found["parents"]] == ["galaxy", "sector"]
    assert found["positions"]["galaxy_pc"] is not None


def test_the_sheet_serves_every_stored_property_of_a_facility_and_it_holds_nothing(galaxy):
    client, ids = galaxy
    sheet = client.get(f"/api/objects/facility:{ids['New Hope']}/sheet").get_json()
    rows = {row["key"]: row for group in sheet["groups"] for row in group["rows"]}
    assert sheet["kind"] == "facility" and rows["kind"]["text"] == "Colony" and rows["placement"]["text"] == "Terrestrial"
    assert rows["rating_crime"]["stored"] is False
    children = client.get(f"/api/objects/facility:{ids['New Hope']}/children").get_json()
    assert children["total"] == 0 and children["items"] == []
    assert client.get("/api/objects/facility:FE81000A2B-0000005-999/sheet").status_code == 404
