# tests/test_rogue_donor_pages.py

"""GEN.1111 part 9: a rogue planet a system threw out is served by the same API, page and sheet as any rogue planet, with the planet's columns."""

import pytest

from planetgen.api.config import Config
from planetgen.db import store
from planetgen.web.app import create_app
from tests.publicids import pid
from tests.test_rogue_donor_save import donor, saved  # noqa: F401  (fixtures)

pytestmark = pytest.mark.fresh_database


@pytest.fixture
def served(saved, mysql_config):  # noqa: F811
    conn, sector_id, _address = saved
    uid = conn.execute("SELECT uid FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()["uid"]
    row = conn.execute("SELECT id FROM rogue_planets WHERE donor_uid = ? ORDER BY event_index LIMIT 1", (uid,)).fetchone()

    class TestConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config

    app = create_app(TestConfig)
    app.testing = True
    return app.test_client(), pid("rogue_planet", row["id"]), row["id"]


def test_the_api_serves_a_donor_rogue_with_its_planet_columns_and_where_it_came_from(served):
    client, printed, _row = served
    response = client.get(f"/api/phenomena/rogue_planet/{printed}")
    assert response.status_code == 200
    body = response.get_json()
    for name in ("planet_class", "surface_temperature_k", "rogue_note", "since_ejection_gy", "ejection_speed_kms"):
        assert name in body, name
    assert body["rogue_note"]


def test_the_page_names_where_a_donor_rogue_came_from(served):
    client, printed, _row = served
    response = client.get(f"/phenomenon/rogue_planet/{printed}")
    assert response.status_code == 200
    page = response.get_data(as_text=True)
    assert "Where It Came From" in page and "Speed It Left At" in page


def test_the_sheet_serves_a_donor_rogue(served):
    client, printed, _row = served
    response = client.get(f"/api/objects/rogue_planet:{printed}/sheet")
    assert response.status_code == 200
    keys = [row["key"] for group in response.get_json()["groups"] for row in group["rows"]] if "groups" in response.get_json() else []
    assert "rogue_note" in keys or "rogue_note" in response.get_data(as_text=True)


def test_the_lists_and_the_sector_page_carry_a_donor_rogue(served):
    client, printed, _row = served
    listed = client.get("/api/phenomena").get_json()["items"]
    assert printed in [row["id"] for row in listed if row["type"] == "rogue_planet"]
    found = client.get(f"/api/phenomena/rogue_planet/{printed}").get_json()
    assert found["donor_name"] and "-" in found["donor_uid"]
    assert client.get(f"/sector/{pid('sector', found['sector_id'])}").status_code == 200

