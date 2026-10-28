"""
MAP.316: the charted regions the Galaxy Map draws at galactic scale and lists. `GET /api/galaxy/charted-regions` has one entry
per 9-block that holds a generated sector (`region_stats`), each with its node key, extent, sector count and capacity and its
share charted (above 0), by share charted; with a naming key in force each is named with its family, without one nothing is.
"""

import pytest

from planetgen.api.config import Config
from planetgen.db import region_stats, store
from planetgen.galaxy import drill, named_regions
from planetgen.names import naming_key
from planetgen.web.app import create_app
from tests.test_region_stats import _addresses
from tests.test_region_wire_names import KEY, EDGE, galaxy, keyed  # noqa: F401 -- fixtures

pytestmark = pytest.mark.db


@pytest.fixture
def client(galaxy):
    class TestConfig(Config):
        MYSQL_CONFIG = galaxy
        WRITE_MYSQL_CONFIG = galaxy
        CONTROL_MYSQL_CONFIG = galaxy

    app = create_app(TestConfig)
    app.testing = True
    return app.test_client()


def _stored(config):
    conn = store.get_connection(config)
    try:
        return region_stats.regions(conn)
    finally:
        conn.close()


def test_every_charted_region_is_listed_once_by_share_charted(galaxy, client, keyed):
    stored = _stored(galaxy)
    listed = client.get("/api/galaxy/charted-regions").get_json()["regions"]
    assert [r["key"] for r in listed] == [drill.format_drill_key(r["block"]) for r in stored]
    assert len({r["key"] for r in listed}) == len(listed) > 0
    shares = [r["share"] for r in listed]
    assert shares == sorted(shares, reverse=True) and all(0 < share <= 1 for share in shares)
    assert sum(r["sector_count"] for r in listed) == len(_addresses())


def test_a_region_carries_its_cell_extent_counts_and_name(galaxy, client, keyed):
    for stored, entry in zip(_stored(galaxy), client.get("/api/galaxy/charted-regions").get_json()["regions"]):
        cell = stored["block"]
        assert (entry["m"], entry["shell"], entry["band"], entry["wedge"]) == (named_regions.REGION_LEVEL, cell.shell, cell.band, cell.wedge)
        assert entry["centre"] == list(stored["extent_centre_pc"]) and entry["radius"] == stored["extent_radius_pc"]
        assert entry["radius"] >= named_regions.framing_radius_pc(cell, EDGE)
        assert (entry["sector_count"], entry["capacity"]) == (stored["sector_count"], stored["capacity"])
        assert entry["share"] == pytest.approx(stored["sector_count"] / stored["capacity"])
        assert entry["name"] == named_regions.region_name(cell, KEY, EDGE)
        assert entry["parent"] == named_regions.family_name(named_regions.family_of(cell), KEY, EDGE)


def test_without_a_naming_key_the_regions_are_listed_unnamed(galaxy, client):
    assert naming_key.active_key() is None
    listed = client.get("/api/galaxy/charted-regions").get_json()["regions"]
    assert listed and all("name" not in entry and "parent" not in entry for entry in listed)
    assert all(entry["radius"] > 0 for entry in listed)
