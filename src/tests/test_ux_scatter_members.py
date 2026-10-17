# tests/test_ux_scatter_members.py

"""
UX.114: the pages show the stored members of an unbuilt scatter row and its system mass. A scattered star has a page of
its own (`/systems/uncharted/<id>`) and a Members cell in the Uncharted stars table; a scattered phenomenon lists its
members on its page and in the Phenomena table; the Galaxy Map labels a star with company; and a single object shows
none of it. Runs against a real throwaway database through the in-process transport.
"""

import pytest

from planetgen.db import query, scatter_systems, store
from planetgen.generation import multistar
from planetgen.util import draw
from planetgen.web.maps.starmap import _neighbor_indicator_data, uncharted_scene_data
from tests.test_api import client  # noqa: F401
from tests.test_galactic_motion import SCATTER_RADIUS_PC, _insert_scatter_point, _seed_scatter_grid
from tests.test_gen_system_motion import EPOCH
from tests.test_uncharted_phenomenon import _build, _hypervelocity_system, _scatter, _sector, db_app  # noqa: F401

POINT = (SCATTER_RADIUS_PC, 0.0, 0.0)


def _tree(seed=3, stars=3):
    return multistar.draw_member_tree(multistar.Placed("star", 1.0), 5e9, draw.Stream(seed), n_stars=stars)


def _star(mysql_config, stars=3, seed=3):
    """One `bright_stars` row with a tree of `stars` members (a single object when `stars` is 1); returns its id and tree."""
    _seed_scatter_grid(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            row_id = _insert_scatter_point(conn, "bright_stars", POINT)
            tree = _tree(seed, stars)
            scatter_systems.store_member_tree(conn, "bright_stars", row_id, tree, EPOCH)
    finally:
        conn.close()
    return row_id, tree


def _stored_tree(mysql_config, table, row_id, seed=3, stars=3):
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            tree = _tree(seed, stars)
            scatter_systems.store_member_tree(conn, table, row_id, tree, EPOCH)
    finally:
        conn.close()
    return tree


def _read(mysql_config, function, *args):
    conn = store.get_connection(mysql_config)
    try:
        return function(conn, *args)
    finally:
        conn.close()


def test_the_member_reader_gives_every_member_and_orbit_with_the_separation_in_au(mysql_config):
    row_id, tree = _star(mysql_config)
    found = _read(mysql_config, query.scatter_member_tree, "bright_stars", row_id)
    assert len(found["members"]) == tree.member_count == 3 and len(found["nodes"]) == 2
    assert sum(member["mass_solar"] for member in found["members"]) == pytest.approx(tree.system_mass_solar)
    assert [member["member_no"] for member in found["members"]] == [0, 1, 2]
    stored = {node.node_no: node for node in tree.nodes}
    for node in found["nodes"]:
        assert node["separation_au"] == pytest.approx(stored[node["node_no"]].separation_km / 149_597_870.7)
        assert node["kind"] in {"close", "wide", "compact", "group"}


def test_a_single_object_has_no_members_and_a_zero_count(mysql_config):
    row_id, _tree_ = _star(mysql_config, stars=1)
    assert _read(mysql_config, query.scatter_member_tree, "bright_stars", row_id) == {"members": [], "nodes": []}
    detail = _read(mysql_config, query.uncharted_system_detail, row_id)
    assert detail["member_count"] == 0 and detail["system_mass_solar"] is None and "members" not in detail


def test_every_uncharted_reader_carries_the_count_and_the_system_mass(mysql_config):
    row_id, tree = _star(mysql_config)
    listed = _read(mysql_config, query.list_uncharted_systems)[0]
    assert (listed["member_count"], listed["system_mass_solar"]) == (3, pytest.approx(tree.system_mass_solar))
    detail = _read(mysql_config, query.uncharted_system_detail, row_id)
    assert detail["member_count"] == 3 and len(detail["members"]) == 3 and len(detail["nodes"]) == 2
    assert detail["system_id"] is None
    address = (listed["shell_index"], listed["band_index"], listed["wedge_index"])
    in_cell = _read(mysql_config, query.bright_stars_in_sector, *address)[0]
    assert in_cell["member_count"] == 3 and in_cell["system_mass_solar"] == pytest.approx(tree.system_mass_solar)
    contents = _read(mysql_config, query.uncharted_sector_contents, *address)
    assert contents["stars"][0]["member_count"] == 3


def test_a_scattered_phenomenon_carries_its_members_until_it_is_built(mysql_config):
    scatter_id = _scatter(mysql_config, "black-hole", "stellar", mass=12.5)
    tree = _stored_tree(mysql_config, "phenomenon_scatter", scatter_id)
    detail = _read(mysql_config, query.scattered_phenomenon_detail, scatter_id)
    assert detail["member_count"] == 3 and detail["system_mass_solar"] == pytest.approx(tree.system_mass_solar)
    assert len(detail["members"]) == 3 and len(detail["nodes"]) == 2
    _build(mysql_config, scatter_id)
    built = _read(mysql_config, query.scattered_phenomenon_detail, scatter_id)
    assert built["member_count"] == 3 and "members" not in built       # a built row's members are its system's stars


def test_the_phenomena_list_rows_carry_the_count(mysql_config):
    with_company = _scatter(mysql_config, "neutron-star", cell=(7, 2, 3))
    _scatter(mysql_config, "neutron-star", cell=(7, 2, 4), record=False)
    _stored_tree(mysql_config, "phenomenon_scatter", with_company)
    rows = {row["id"]: row for row in _read(mysql_config, query.list_phenomena) if row["scattered"]}
    assert rows[with_company]["member_count"] == 3
    assert [row["member_count"] for row in rows.values()].count(0) == len(rows) - 1


def test_the_api_serves_one_scattered_star_with_its_members(client, mysql_config):  # noqa: F811
    row_id, tree = _star(mysql_config)
    body = client.get(f"/api/uncharted-systems/{row_id}").get_json()
    assert body["id"] == row_id and body["member_count"] == 3 and len(body["members"]) == 3
    assert body["system_mass_solar"] == pytest.approx(tree.system_mass_solar)
    assert client.get("/api/uncharted-systems/999999").status_code == 404
    listed = client.get("/api/uncharted-systems").get_json()["items"][0]
    assert listed["member_count"] == 3


def test_the_uncharted_star_page_lists_the_members_and_the_system_mass(db_app, mysql_config):  # noqa: F811
    row_id, tree = _star(mysql_config)
    html = db_app.test_client().get(f"/systems/uncharted/{row_id}").get_data(as_text=True)
    assert "Uncharted" in html and "System mass" in html and "Members" in html
    assert html.count("<tr><td>") >= 3 + 2          # three members and two orbits
    assert "Close pair" in html or "Wide pair" in html or "Compact pair" in html


def test_the_page_of_a_single_star_has_no_members_panel(db_app, mysql_config):  # noqa: F811
    row_id, _tree_ = _star(mysql_config, stars=1)
    html = db_app.test_client().get(f"/systems/uncharted/{row_id}").get_data(as_text=True)
    assert "Star Data" in html and "System mass" not in html and 'id="members-heading"' not in html


def test_an_unknown_star_is_a_404_and_a_built_one_goes_to_its_system(db_app, mysql_config):  # noqa: F811
    assert db_app.test_client().get("/systems/uncharted/999999").status_code == 404
    row_id, _tree_ = _star(mysql_config)
    system_id = _hypervelocity_system(mysql_config, _sector(mysql_config))
    conn = store.get_connection(mysql_config)
    conn.execute("UPDATE bright_stars SET star_system_id = ? WHERE id = ?", (system_id, row_id))
    conn.commit()
    conn.close()
    response = db_app.test_client().get(f"/systems/uncharted/{row_id}")
    assert response.status_code == 302 and "/system/" in response.headers["Location"]


def test_the_uncharted_stars_table_has_a_members_cell_and_links_to_the_star_page(db_app, mysql_config):  # noqa: F811
    row_id, _tree_ = _star(mysql_config)
    html = db_app.test_client().get("/systems?uncharted=1").get_data(as_text=True)
    assert "Members" in html and "3 members," in html and f'href="/systems/uncharted/{row_id}"' in html


def test_the_phenomena_table_and_the_uncharted_phenomenon_page_show_the_members(db_app, mysql_config):  # noqa: F811
    scatter_id = _scatter(mysql_config, "black-hole", "stellar", mass=12.5)
    _stored_tree(mysql_config, "phenomenon_scatter", scatter_id)
    client = db_app.test_client()
    assert "3 members," in client.get("/phenomena").get_data(as_text=True)
    page = client.get(f"/uncharted-phenomenon/{scatter_id}").get_data(as_text=True)
    assert 'id="members-heading"' in page and "System mass" in page


def test_the_map_labels_a_star_with_company_and_the_neighbour_panel_counts_its_members():
    star = {"id": 1, "x": 1.0, "y": 2.0, "z": 3.0, "luminosity_sol": 900.0, "temperature_k": 9000.0,
            "radius_sol": 2.0, "star_type": "B2V", "member_count": 3, "system_mass_solar": 4.5}
    alone = dict(star, id=2, member_count=0, system_mass_solar=None)
    scene = uncharted_scene_data((7, 2, 3), "A7", (0.0, 0.0, 0.0), 4.0,
                                 {"stars": [star, alone], "scattered": []})
    assert "3 members," in scene["stars"][0]["typeLabel"] and "typeLabel" not in scene["stars"][1]
    neighbor = {"direction_pc": (1.0, 0.0, 0.0), "shell_index": 7, "band_index": 2, "wedge_index": 4,
                "designation": "A8", "exists": False, "bright_stars": [star, alone]}
    shown = _neighbor_indicator_data(lambda *args, **kwargs: "/", neighbor)["brightStars"]
    assert "3 members," in shown[0] and "member" not in shown[1]
