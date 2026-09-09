# tests/test_api.py

"""
End-to-end tests for the read-only Flask API (`src/api/`) against a real,
throwaway MySQL database (see `conftest.py`'s `mysql_config` fixture) --
every test here is skipped, not failed, when no MySQL test server is
configured/reachable.

Seeds a small, real sector (via `stellarObjects._db.save_sector`, the same
path `sectorGen.py` uses) rather than hand-building rows, so these tests
exercise the exact same read path (`queryDb.py`/`stellarObjects._db.load_sector`/
`load_star_system`) production traffic does.
"""

import pytest

from api.app import create_app
from stellarObjects import _db
from stellarObjects.config import SystemConfig
from stellarObjects.spaceSector import SpaceSector
from stellarObjects.systemData import StarSystem


@pytest.fixture
def seeded_sector(mysql_config):
    """Saves one small sector (two systems, one G2V, one M5V) and returns
    (mysql_config, sector_id, system_ids) for tests to build requests
    against."""
    sector = SpaceSector("Test Sector", edge_ly=10.0)

    cfg_a = SystemConfig()
    cfg_a.STAR_TYPE = "G2V"
    cfg_a.PLANETS = False
    system_a = StarSystem(system_config=cfg_a)
    sector.add_system(system_a, position=(1.0, 1.0, 1.0), system_config=cfg_a)

    cfg_b = SystemConfig()
    cfg_b.STAR_TYPE = "M5V"
    cfg_b.PLANETS = False
    system_b = StarSystem(system_config=cfg_b)
    sector.add_system(system_b, position=(-2.0, 0.5, 3.0), system_config=cfg_b)

    sector_id = _db.save_sector(sector, config=mysql_config)

    conn = _db.get_connection(mysql_config)
    try:
        rows = conn.execute(
            "SELECT id, name FROM star_systems WHERE sector_id = ? ORDER BY id", (sector_id,)
        ).fetchall()
    finally:
        conn.close()

    return mysql_config, sector_id, [row["id"] for row in rows]


@pytest.fixture
def client(mysql_config):
    class TestConfig:
        MYSQL_CONFIG = mysql_config

    app = create_app(TestConfig)
    app.testing = True
    return app.test_client()


def test_health_ok(client):
    response = client.get("/api/health")
    assert response.status_code == 200
    assert response.get_json() == {"status": "ok"}


def test_sectors_lists_seeded_sector(client, seeded_sector):
    _config, sector_id, _system_ids = seeded_sector

    response = client.get("/api/sectors")
    assert response.status_code == 200
    body = response.get_json()
    assert body["total"] >= 1
    assert body["limit"] == 100
    assert body["offset"] == 0
    ids = {entry["id"] for entry in body["items"]}
    assert sector_id in ids


def test_sectors_pagination(client, seeded_sector):
    response = client.get("/api/sectors?limit=1&offset=0")
    assert response.status_code == 200
    body = response.get_json()
    assert len(body["items"]) == 1
    assert body["limit"] == 1
    assert body["offset"] == 0


def test_sectors_rejects_invalid_limit(client):
    response = client.get("/api/sectors?limit=not-a-number")
    assert response.status_code == 400
    assert "error" in response.get_json()


def test_sector_detail_found_and_not_found(client, seeded_sector):
    _config, sector_id, _system_ids = seeded_sector

    response = client.get(f"/api/sectors/{sector_id}")
    assert response.status_code == 200
    body = response.get_json()
    assert body["name"] == "Test Sector"

    response = client.get("/api/sectors/999999999")
    assert response.status_code == 404
    assert "error" in response.get_json()


def test_systems_filters_by_star_type_and_sector(client, seeded_sector):
    _config, sector_id, system_ids = seeded_sector

    response = client.get(f"/api/systems?sector_id={sector_id}")
    assert response.status_code == 200
    body = response.get_json()
    assert body["total"] == len(system_ids)

    response = client.get(f"/api/systems?sector_id={sector_id}&star_type=G")
    assert response.status_code == 200
    body = response.get_json()
    assert body["total"] == 1
    assert body["items"][0]["is_binary"] == 0


def test_systems_rejects_invalid_sector_id(client):
    response = client.get("/api/systems?sector_id=not-an-int")
    assert response.status_code == 400
    assert "error" in response.get_json()


def test_system_detail_found_and_not_found(client, seeded_sector):
    _config, _sector_id, system_ids = seeded_sector

    response = client.get(f"/api/systems/{system_ids[0]}")
    assert response.status_code == 200
    body = response.get_json()
    assert "star" in body

    response = client.get("/api/systems/999999999")
    assert response.status_code == 404


def test_systems_near_requires_radius(client, seeded_sector):
    _config, _sector_id, system_ids = seeded_sector

    response = client.get(f"/api/systems/{system_ids[0]}/near")
    assert response.status_code == 400

    response = client.get(f"/api/systems/{system_ids[0]}/near?radius=-5")
    assert response.status_code == 400

    response = client.get(f"/api/systems/{system_ids[0]}/near?radius=1000")
    assert response.status_code == 200
    body = response.get_json()
    assert isinstance(body, list)
    other_ids = {entry["id"] for entry in body}
    assert system_ids[1] in other_ids


def test_unmatched_route_returns_json_404(client):
    response = client.get("/api/no-such-route")
    assert response.status_code == 404
    assert response.get_json() == {"error": "not found"}
