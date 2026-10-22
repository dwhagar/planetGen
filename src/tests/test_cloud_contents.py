# tests/test_cloud_contents.py

"""
UX.144: every item page says where the item is (the places above it and the nebula or remnant it sits in), and a nebula
or remnant lists what sits inside it.
"""

import pytest

from planetgen.api.config import Config
from planetgen.db import cloudcontents, store
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.generation.system import StarSystem
from planetgen.web.app import create_app
from tests.publicids import pid


@pytest.fixture
def galaxy(mysql_config):
    """A sector with two systems (the first inside a nebula) and a rogue planet inside it; a second nebula holds
    nothing."""
    sector = SpaceSector("Cloud Sector", edge_ly=20.0)
    cfg = SystemConfig()
    cfg.BINARY_SYSTEM = False
    sector.add_system(StarSystem(system_config=cfg), position=(1.0, 1.0, 1.0), system_config=cfg)
    sector.add_system(StarSystem(system_config=cfg), position=(-2.0, 0.5, 3.0), system_config=cfg)
    sector.add_phenomenon(RoguePlanet(SystemConfig()), "rogue-planet")
    sector_row = store.save_sector(sector, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            placement = {"center_x_pc": 1.0, "center_y_pc": 2.0, "center_z_pc": 3.0, "galactic_radius_pc": 3.7}
            veil = store.insert_nebula(conn, Nebula(SystemConfig(), name="Test Veil"), sector_id=sector_row,
                                       placement=placement)
            empty = store.insert_nebula(conn, Nebula(SystemConfig(), name="Empty Veil"), sector_id=sector_row,
                                        placement=placement)
            store.assign_uids(conn, phenomenon=("nebulae", veil))
            store.assign_uids(conn, phenomenon=("nebulae", empty))
            inside, outside = [row["id"] for row in conn.execute("SELECT id FROM star_systems ORDER BY id").fetchall()]
            conn.execute("UPDATE star_systems SET inside_nebula_id = ? WHERE id = ?", (veil, inside))
            conn.execute("UPDATE rogue_planets SET inside_nebula_id = ?", (veil,))
        ids = {"veil": veil, "empty": empty, "inside": inside, "outside": outside}
        printed = {"veil": pid("nebula", veil, conn), "empty": pid("nebula", empty, conn),
                   "inside": pid("system", inside, conn), "outside": pid("system", outside, conn)}
    finally:
        conn.close()
    return ids, printed


@pytest.fixture
def client(mysql_config, galaxy):
    class TestConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config

    app = create_app(TestConfig)
    app.testing = True
    return app.test_client()


def test_a_cloud_lists_what_is_inside_it(mysql_config, galaxy):
    ids, _printed = galaxy
    conn = store.get_connection(mysql_config)
    try:
        found = cloudcontents.cloud_contents(conn, "nebula", ids["veil"])
        assert found["counts"] == {"system": 1, "rogue_planet": 1} and found["total"] == 2
        assert [item["kind"] for item in found["items"]] == ["system", "rogue_planet"]
        assert cloudcontents.cloud_contents(conn, "nebula", ids["empty"]) == {"counts": {}, "items": [], "total": 0}
        assert cloudcontents.cloud_contents(conn, "nebula", ids["veil"], limit=1)["items"][0]["kind"] == "system"
        with pytest.raises(ValueError):
            cloudcontents.cloud_contents(conn, "black_hole", 1)
    finally:
        conn.close()


def test_the_contents_route_speaks_object_ids(client, galaxy):
    _ids, printed = galaxy
    found = client.get(f"/api/phenomena/nebula/{printed['veil']}/contents").get_json()
    assert found["total"] == 2
    assert printed["inside"] in [item["id"] for item in found["items"] if item["kind"] == "system"]
    assert client.get(f"/api/phenomena/nebula/{printed['empty']}/contents").get_json()["total"] == 0
    assert client.get(f"/api/phenomena/rogue_planet/{printed['veil']}/contents").status_code == 404
    assert client.get("/api/phenomena/nebula/FFFFFFFFFF-0000000-000/contents").status_code == 404


def test_a_nebula_page_lists_its_contents_and_an_empty_one_says_so(client, galaxy):
    _ids, printed = galaxy
    page = client.get(f"/phenomenon/nebula/{printed['veil']}").get_data(as_text=True)
    assert "What is inside it" in page and f"/system/{printed['inside']}" in page
    assert "Star systems (1)" in page and "Rogue planets (1)" in page
    empty = client.get(f"/phenomenon/nebula/{printed['empty']}").get_data(as_text=True)
    assert "Nothing generated sits inside it." in empty
    # Another kind of phenomenon contains nothing and has no such block.
    rogue = next(row for row in client.get("/api/phenomena").get_json()["items"] if row["type"] == "rogue_planet")
    assert "What is inside it" not in client.get(f"/phenomenon/rogue_planet/{rogue['id']}").get_data(as_text=True)


def test_a_phenomenon_page_has_a_location_block(client, galaxy):
    _ids, printed = galaxy
    rogue = next(row for row in client.get("/api/phenomena").get_json()["items"] if row["type"] == "rogue_planet")
    page = client.get(f"/phenomenon/rogue_planet/{rogue['id']}").get_data(as_text=True)
    assert 'id="location"' in page and "Cloud Sector" in page
    assert f'Inside: <a href="/phenomenon/nebula/{printed["veil"]}">Test Veil</a>' in page


def test_a_body_inside_a_nebula_says_where_it_is(client, galaxy):
    _ids, printed = galaxy
    system = client.get(f"/api/systems/{printed['inside']}").get_json()
    # A system's planets are drawn, so the star is the body that is always there; any planet is checked too.
    bodies = [("star", system["stars"][0])] + [("planet", planet) for planet in system["planets"][:1]]
    for kind, body in bodies:
        page = client.get(f"/object/{kind}:{body['id']}").get_data(as_text=True)
        assert 'id="location"' in page
        assert 'href="/galaxy"' in page and "Cloud Sector" in page and f'href="/system/{printed["inside"]}"' in page
        assert f'Inside: <a href="/phenomenon/nebula/{printed["veil"]}">Test Veil</a>' in page


def test_a_body_in_open_space_says_it_is_in_none(client, galaxy):
    _ids, printed = galaxy
    star = client.get(f"/api/systems/{printed['outside']}").get_json()["stars"][0]
    page = client.get(f"/object/star:{star['id']}").get_data(as_text=True)
    assert "Not inside any nebula or supernova remnant." in page and "Test Veil" not in page
    system_page = client.get(f"/system/{printed['outside']}").get_data(as_text=True)
    assert 'id="location"' in system_page and "Not inside any nebula or supernova remnant." in system_page
