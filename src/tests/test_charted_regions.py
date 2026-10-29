"""
MAP.316: the charted regions the Galaxy Map draws at galactic scale and lists. `GET /api/galaxy/charted-regions` has one entry
per 9-block that holds a generated sector (`region_stats`), each with its node key, extent, sector count and capacity and its
share charted (above 0), by share charted; with a naming key in force each is named with its family, without one nothing is.
"""

import pytest

from planetgen.api.config import Config
from planetgen.db import region_stats, store
from planetgen.galaxy import drill, named_regions, spherical
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


def _inhabit(config, rows):
    """Puts an `inhabited` interest row `(system id, n, score, address)` at the middle of each sector (no star system behind it)."""
    conn = store.get_connection(config)
    try:
        with conn:
            conn.execute("SET SESSION FOREIGN_KEY_CHECKS = 0")
            for system_id, n, score, address in rows:
                x, y, z = spherical.sector_position_pc(*address, EDGE)
                conn.execute("INSERT INTO system_interests (star_system_id, reason, n, best_score, galaxy_x_pc, galaxy_y_pc,"
                             " galaxy_z_pc) VALUES (?, ?, ?, ?, ?, ?, ?)", (system_id, "inhabited", n, score, x, y, z))
            conn.execute("SET SESSION FOREIGN_KEY_CHECKS = 1")
    finally:
        conn.close()


def test_a_region_with_inhabited_worlds_says_how_many_and_how_good_the_best_is(galaxy, client, keyed):
    addresses = _addresses()
    _inhabit(galaxy, [(9001, 2, 0.8, addresses[0]), (9002, 1, 1.2, addresses[1]), (9003, 1, None, addresses[-1])])
    listed = {entry["key"]: entry for entry in client.get("/api/galaxy/charted-regions").get_json()["regions"]}
    first = listed[drill.format_drill_key(region_stats.region_cell(addresses[0]))]
    assert first["inhabited"] == {"systems": 2, "worlds": 3, "score": 1.2, "tier": "Blue"}
    far = listed[drill.format_drill_key(region_stats.region_cell(addresses[-1]))]
    assert far["inhabited"] == {"systems": 1, "worlds": 1}
    assert [entry for entry in listed.values() if "inhabited" in entry] == [first, far]


def test_a_region_with_no_inhabited_world_has_no_inhabited_entry(galaxy, client, keyed):
    listed = client.get("/api/galaxy/charted-regions").get_json()["regions"]
    assert listed and all("inhabited" not in entry for entry in listed)
