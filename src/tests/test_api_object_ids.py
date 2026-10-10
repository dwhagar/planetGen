# tests/test_api_object_ids.py

"""
API.23, stage 1: sectors, systems and phenomena speak object IDs in the API, never row numbers.
"""

import re

import pytest

from planetgen.api.config import Config
from planetgen.db import store
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
