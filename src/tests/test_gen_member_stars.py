# tests/test_gen_member_stars.py

"""
GEN.630: a system holds N stars. The stars beyond the pair are `member_stars`, each a `member` row of `stars` that says
the node and slot it hangs in; a system of three round-trips through the JSON form and through save and load, and every
reader that took a system to have one or two stars handles the third.
"""

import pytest

from planetgen.db import orbit_trees, query, store
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.system import StarSystem
from planetgen.galaxy import system_position
from planetgen.population import facilities
from tests.member_systems import binary, with_member

pytestmark = pytest.mark.fresh_database


def _saved(mysql_config, wide, planets=False):
    system, cfg = binary(wide, planets)
    with_member(system, cfg)
    conn = store.get_connection(mysql_config)
    system_id = store.insert_star_system(conn, system, cfg)
    conn.commit()
    return conn, system, system_id


@pytest.mark.parametrize("wide", [False, True])
def test_a_system_of_three_goes_through_the_json_form(wide):
    system, cfg = binary(wide)
    with_member(system, cfg)
    again = StarSystem.from_dict(system.to_dict())
    assert len(again.stars) == 3 and len(again.member_stars) == 1
    assert again.member_stars[0].mass == pytest.approx(system.member_stars[0].mass)
    assert again.star_slots == system.star_slots and again.orbit_nodes == system.orbit_nodes


def test_a_system_without_members_has_none_in_the_json_form():
    system, _cfg = binary(True)
    data = system.to_dict()
    assert data["member_stars"] == [] and data["star_slots"] is None
    del data["member_stars"], data["star_slots"]
    again = StarSystem.from_dict({**data, "schema_version": 5})
    assert again.member_stars == () and again.star_slots is None and len(again.stars) == 2


@pytest.mark.parametrize("wide", [False, True])
def test_a_system_of_three_is_saved_with_a_member_row_and_every_star_in_its_slot(mysql_config, wide):
    conn, system, system_id = _saved(mysql_config, wide)
    try:
        rows = conn.execute("SELECT role, parent_node_id, slot, offset_x_au, offset_y_au, offset_z_au, mass_kg"
                            " FROM stars WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchall()
        assert [r["role"] for r in rows] == ["primary", "secondary", "member"]
        assert [(r["parent_node_id"], r["slot"]) for r in rows] == [(0, 0), (0, 1), (1, 1)]
        # the member sits on the far side of the pair's centre from the pair, in proportion to the masses
        member = rows[2]
        distance_au = sum(member[f"offset_{axis}_au"] ** 2 for axis in "xyz") ** 0.5
        pair_mass = rows[0]["mass_kg"] + rows[1]["mass_kg"]
        assert 0 < distance_au < 400.0 * 1.3 * (pair_mass / (pair_mass + member["mass_kg"])) + 1.0
        centre = [sum(r[f"offset_{axis}_au"] * r["mass_kg"] for r in rows) for axis in "xyz"]
        assert all(abs(c) < 1e-6 * sum(r["mass_kg"] for r in rows) * 400.0 for c in centre)
        assert orbit_trees.read_nodes(conn, system_id) == system.orbit_nodes
    finally:
        conn.close()


@pytest.mark.parametrize("wide", [False, True])
def test_a_system_of_three_loads_with_its_member(mysql_config, wide):
    conn, system, system_id = _saved(mysql_config, wide)
    try:
        loaded = store.load_star_system(conn, system_id)
        assert len(loaded.stars) == 3 and len(loaded.member_stars) == 1
        assert loaded.member_stars[0].mass == pytest.approx(system.member_stars[0].mass, rel=1e-9)
        assert loaded.member_stars[0].db_id is not None
        assert loaded.star_slots == ((0, 0), (0, 1), (1, 1))
        assert loaded.orbit_nodes == system.orbit_nodes
        assert loaded.binary_type == system.binary_type
    finally:
        conn.close()


def test_a_pair_has_its_slots_and_no_members(mysql_config):
    system, cfg = binary(True)
    conn = store.get_connection(mysql_config)
    try:
        system_id = store.insert_star_system(conn, system, cfg)
        rows = conn.execute("SELECT role, slot FROM stars WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchall()
        assert [(r["role"], r["slot"]) for r in rows] == [("primary", 0), ("secondary", 1)]
        loaded = store.load_star_system(conn, system_id)
        assert loaded.member_stars == () and loaded.star_slots is None and len(loaded.stars) == 2
    finally:
        conn.close()


def test_system_detail_lists_every_star_and_the_scene_keeps_the_member_off_the_primary(mysql_config):
    from planetgen.web.maps import systemscene
    conn, _system, system_id = _saved(mysql_config, True)
    try:
        detail = query.system_detail(conn, system_id)
        assert [star["role"] for star in detail["stars"]] == ["primary", "secondary", "member"]
        scene = systemscene.build_scene(conn, system_id)
        member = scene["stars"][2]
        assert member["role"] == "member" and any(member["position_km"])
        assert [m["slot"] for m in scene["tree"]["members"]] == [0, 1, 1]
        assert [m["parent_node_no"] for m in scene["tree"]["members"]] == [0, 0, 1]
    finally:
        conn.close()


def test_system_detail_carries_the_orbit_tree_and_where_each_star_hangs(mysql_config):
    """MAP.213: the detail lists the `orbit_nodes` in node order and each star's node and slot, so a page can name what an orbit joins."""
    conn, _system, system_id = _saved(mysql_config, True)
    try:
        detail = query.system_detail(conn, system_id)
        assert [(node["node_no"], node["parent_node_no"]) for node in detail["orbits"]] == [(0, 1), (1, None)]
        assert [(star["parent_node_no"], star["slot"]) for star in detail["stars"]] == [(0, 0), (0, 1), (1, 1)]
        assert all(node["separation_km"] > 0 and node["period_years"] > 0 for node in detail["orbits"])
    finally:
        conn.close()


def test_a_single_star_has_no_orbits_in_its_detail(mysql_config):
    from planetgen.generation.config import SystemConfig
    cfg = SystemConfig()
    cfg.BINARY_SYSTEM = False
    system = StarSystem(system_config=cfg)
    conn = store.get_connection(mysql_config)
    try:
        system_id = store.insert_star_system(conn, system, cfg)
        conn.commit()
        detail = query.system_detail(conn, system_id)
        assert detail["orbits"] == [] and detail["stars"][0]["parent_node_no"] is None
    finally:
        conn.close()


def test_a_sector_lists_all_three_stars_in_order(mysql_config):
    system, cfg = binary(False)
    with_member(system, cfg)
    sector = SpaceSector("Member Sector", edge_ly=10.0)
    sector.add_system(system, position=(1.0, 1.0, 1.0), system_config=cfg)
    sector_id = store.save_sector(sector, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        systems = query.sector_detail(conn, sector_id)["systems"]
        assert [s["role"] for s in systems[0]["stars"]] == ["primary", "secondary", "member"]
    finally:
        conn.close()


def test_the_member_is_placed_in_the_galaxy_with_the_pair():
    system, cfg = binary(True)
    with_member(system, cfg)
    offsets = system_position._star_offsets_au(system)
    assert len(offsets) == 3
    assert any(offsets[2][1])
    total_mass = sum(star.mass for star in system.stars)
    for axis in range(3):
        assert sum(star.mass * offset[axis] for star, offset in offsets) == pytest.approx(0.0, abs=1e-6 * total_mass * 400.0)


def test_a_facility_around_a_close_pair_does_not_count_the_member_in_the_pair():
    stars = [{"id": 1, "role": "primary", "mass_kg": 2.0, "radius_km": 10.0, "heliosphere_radius_km": 5.0},
             {"id": 2, "role": "secondary", "mass_kg": 1.0, "radius_km": 7.0, "heliosphere_radius_km": 5.0},
             {"id": 3, "role": "member", "mass_kg": 100.0, "radius_km": 1000.0, "heliosphere_radius_km": 5.0}]
    mass, radius, _sphere = facilities.star_host(stars, 1, "close", 3.0, 9.0)
    assert mass == 3.0 and radius == 13.0
    assert facilities.star_host(stars, 3, "wide")[:2] == (100.0, 1000.0)



def test_a_saved_system_numbers_its_members_in_the_body_slots(mysql_config):
    conn, _system, system_id = _saved(mysql_config, True)
    try:
        slots = store._body_slots(conn, system_id)
        stars = sorted(value for (table, _id), value in slots.items() if table == "stars")
        assert stars == [("star", "member", 0), ("star", "primary", 0), ("star", "secondary", 0)]
    finally:
        conn.close()


@pytest.fixture
def web_app(mysql_config, monkeypatch):
    from planetgen.api.config import Config
    from planetgen.web.app import create_app
    from planetgen.web.lib import apiclient

    class RealConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "test-secret"

    def no_http(*args, **kwargs):
        raise AssertionError("HTTP transport used inside a Flask request")
    monkeypatch.setattr(apiclient, "_http_transport", no_http)
    application = create_app(RealConfig)
    application.testing = True
    return application


@pytest.mark.parametrize("wide", [False, True])
def test_the_system_page_and_its_map_draw_a_system_of_three(web_app, mysql_config, wide):
    from tests.publicids import pid
    system, cfg = binary(wide, planets=True)
    with_member(system, cfg)
    sector = SpaceSector("Member Page Sector", edge_ly=10.0)
    sector.add_system(system, position=(1.0, 1.0, 1.0), system_config=cfg)
    sector_id = store.save_sector(sector, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        system_id = conn.execute("SELECT id FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()["id"]
        member_name = conn.execute("SELECT name FROM stars WHERE star_system_id = ? AND role = 'member'",
                                   (system_id,)).fetchone()["name"]
    finally:
        conn.close()
    response = web_app.test_client().get(f"/system/{pid('system', system_id)}")
    assert response.status_code == 200
    assert member_name in response.get_data(as_text=True)


def test_the_system_page_lists_the_orbits_of_a_system_of_three(web_app, mysql_config):
    """MAP.213: a triple's page has an "Orbits of the stars" panel with a row for the inner pair and one for the orbit about it."""
    import re
    from tests.publicids import pid
    system, cfg = binary(True, planets=False)
    with_member(system, cfg)
    sector = SpaceSector("Orbit Page Sector", edge_ly=10.0)
    sector.add_system(system, position=(1.0, 1.0, 1.0), system_config=cfg)
    sector_id = store.save_sector(sector, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        system_id = conn.execute("SELECT id FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()["id"]
        names = [row["name"] for row in conn.execute("SELECT name FROM stars WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchall()]
    finally:
        conn.close()
    html = web_app.test_client().get(f"/system/{pid('system', system_id)}").get_data(as_text=True)
    panel = re.search(r'<section class="panel" id="system-orbits".*?</section>', html, re.S)
    assert panel, "no orbits panel"
    assert len(re.findall(r"<tr><td>", panel.group(0))) == 2
    assert f"{names[0]} and {names[1]}" in panel.group(0) and f"orbit 1 and {names[2]}" in panel.group(0)


def test_the_system_page_of_a_single_star_has_no_orbits_panel(web_app, mysql_config):
    from planetgen.generation.config import SystemConfig
    from tests.publicids import pid
    cfg = SystemConfig()
    cfg.BINARY_SYSTEM = False
    system = StarSystem(system_config=cfg)
    sector = SpaceSector("Lone Page Sector", edge_ly=10.0)
    sector.add_system(system, position=(1.0, 1.0, 1.0), system_config=cfg)
    sector_id = store.save_sector(sector, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        system_id = conn.execute("SELECT id FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()["id"]
    finally:
        conn.close()
    html = web_app.test_client().get(f"/system/{pid('system', system_id)}").get_data(as_text=True)
    assert 'id="system-orbits"' not in html
