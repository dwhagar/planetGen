# tests/test_api_retired_ids.py

"""GEN.1139: an ID an object gave up (`designations`, GEN.1138) still answers in the API: a read is a 301 to the object's current ID with the
rest of the URL kept, a write is a 409 that does nothing, the ID of a deleted object is a 410, and an ID nobody ever had is a 404 as before."""

import pytest

from planetgen.db import designation_rows as rows
from planetgen.db import store
from planetgen.galaxy import object_uid as ou
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.generation.system import StarSystem
from tests.test_api import admin_client, client, default_admin_client, first_admin_password  # noqa: F401

pytestmark = pytest.mark.fresh_database


@pytest.fixture
def galaxy(mysql_config, admin_client):  # noqa: F811
    sector = SpaceSector("Retired Sector", edge_ly=20.0)
    config = SystemConfig()
    config.PLANETS = True
    config.BINARY_SYSTEM = False
    sector.add_system(StarSystem(system_config=config), position=(1.0, 1.0, 1.0), system_config=config)
    sector.add_phenomenon(RoguePlanet(SystemConfig()), "rogue-planet")
    store.save_sector(sector, config=mysql_config)
    return mysql_config, admin_client


def _renumber(config, table):
    conn = store.get_connection(config)
    try:
        row = conn.execute(f"SELECT id, uid, name FROM {table} ORDER BY id LIMIT 1").fetchone()
        old = ou.from_bytes(row["uid"])
        new = old + (5000 << 12)
        with conn:
            rows.renumber(conn, table, row["id"], new, "admin renumber")
        return row, ou.format_id(old), ou.format_id(new)
    finally:
        conn.close()


def test_a_read_of_a_retired_system_id_is_a_301_to_the_current_id_with_the_rest_of_the_url(galaxy):
    config, client = galaxy
    row, old, new = _renumber(config, "star_systems")
    assert client.get(f"/api/systems/{new}").status_code == 200
    response = client.get(f"/api/systems/{old}")
    assert response.status_code == 301
    assert response.headers["Location"] == f"/api/systems/{new}"
    assert response.headers["X-Designation-Formerly"] == old
    body = response.get_json()
    assert (body["id"], body["retired_id"]) == (new, old)
    [former] = body["former_designations"]
    assert (former["designation"], former["kind"], former["reason"]) == (old, "uid", "admin renumber")
    deeper = client.get(f"/api/systems/{old.lower()}/scene?view=1")
    assert deeper.status_code == 301 and deeper.headers["Location"] == f"/api/systems/{new}/scene?view=1"


def test_a_write_against_a_retired_id_is_a_409_and_changes_nothing(galaxy):
    config, client = galaxy
    row, old, new = _renumber(config, "star_systems")
    response = client.patch(f"/api/systems/{old}", json={"name": "Never Applied"})
    assert response.status_code == 409
    assert response.get_json()["id"] == new and "nothing was changed" in response.get_json()["error"]
    conn = store.get_connection(config)
    try:
        assert conn.execute("SELECT name FROM star_systems WHERE id = ?", (row["id"],)).fetchone()["name"] == row["name"]
        assert conn.execute("SELECT COUNT(*) AS n FROM designations WHERE designation_kind = 'name'").fetchone()["n"] == 0
    finally:
        conn.close()
    assert client.patch(f"/api/systems/{new}", json={"name": "Applied"}).status_code == 200


def test_a_retired_phenomenon_id_redirects_within_its_type(galaxy):
    config, client = galaxy
    _row, old, new = _renumber(config, "rogue_planets")
    response = client.get(f"/api/phenomena/rogue_planet/{old}")
    assert response.status_code == 301 and response.headers["Location"] == f"/api/phenomena/rogue_planet/{new}"


def test_the_id_of_a_deleted_object_is_a_410_for_a_read_and_a_write(galaxy):
    config, client = galaxy
    row, old, _new = _renumber(config, "star_systems")
    conn = store.get_connection(config)
    try:
        with conn:
            conn.execute("UPDATE designations SET gone_epoch = 7 WHERE object_kind = 'star_systems' AND object_row_id = ?", (row["id"],))
    finally:
        conn.close()
    for response in (client.get(f"/api/systems/{old}"), client.patch(f"/api/systems/{old}", json={"name": "Gone"})):
        assert response.status_code == 410
        assert response.get_json()["gone_epoch"] == 7 and response.get_json()["retired_id"] == old
        assert "Location" not in response.headers


def test_an_id_nobody_had_or_a_retired_id_asked_as_another_kind_is_still_a_404(galaxy):
    config, client = galaxy
    _row, old, _new = _renumber(config, "star_systems")
    never = ou.format_id(ou.parse_id(old) + (9000 << 12))
    assert client.get(f"/api/systems/{never}").status_code == 404
    assert client.get(f"/api/planets/{old}/life-air").status_code == 404
    assert client.get("/api/systems/1").status_code == 404
    assert client.patch(f"/api/systems/{never}", json={"name": "Nope"}).status_code == 404


def test_a_page_asked_for_under_a_retired_id_is_a_301_to_the_page_under_the_current_id(galaxy):
    config, client = galaxy
    _row, old, new = _renumber(config, "star_systems")
    assert client.get(f"/system/{new}").status_code == 200
    response = client.get(f"/system/{old}")
    assert response.status_code == 301 and response.headers["Location"] == f"/system/{new}"
    deeper = client.get(f"/system/{old.lower()}/scene?x=1")
    assert deeper.status_code == 301 and deeper.headers["Location"] == f"/system/{new}/scene?x=1"


def test_a_page_of_a_deleted_object_is_a_410_and_one_nobody_had_is_a_404(galaxy):
    config, client = galaxy
    row, old, _new = _renumber(config, "star_systems")
    never = ou.format_id(ou.parse_id(old) + (9000 << 12))
    conn = store.get_connection(config)
    try:
        with conn:
            conn.execute("UPDATE designations SET gone_epoch = 7 WHERE object_kind = 'star_systems' AND object_row_id = ?", (row["id"],))
    finally:
        conn.close()
    gone = client.get(f"/system/{old}")
    assert gone.status_code == 410 and b"Gone" in gone.data
    assert client.get(f"/system/{never}").status_code == 404
