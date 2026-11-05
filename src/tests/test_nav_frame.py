# tests/test_nav_frame.py

"""
NAV.73: a routed system sits where the database stores it. The route search used to add a system's sector-local offset to
its sector's center without turning it by the sector's axes (`frame.local_to_galaxy_pc`), so a system of a sector away
from the +X axis sat several light-years from its stored `galaxy_x_pc`, `galaxy_y_pc` and `galaxy_z_pc`.
"""

import pytest

from planetgen.db import corridor, query, store
from planetgen.galaxy.frame import galaxy_to_local_pc
from planetgen.physics.units import ly_to_pc, milliparsecs_to_ly, pc_to_ly
from planetgen.util import detmath
from tests.test_nav_route_graph import _all_systems_ly, _cloud_local
from tests.test_route_edge_cases import _save_grid_sector

FAR_WEDGE = 2000
"""int: A wedge at a longitude of about 58 degrees, where the sector's axes are turned far from the galaxy's."""


def _nav(config, west, east):
    conn = store.get_connection(config)
    try:
        return query.nav_between(conn, west, east)
    finally:
        conn.close()


def test_the_corridor_places_a_system_where_the_database_stores_it(mysql_config):
    ids = _save_grid_sector(mysql_config, "Far", FAR_WEDGE, _cloud_local(3, count=6, spread=6.0))
    stored = _all_systems_ly(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        point = stored[ids[0]]
        found = corridor.positions_near_segment(conn, point, point, 200.0)
    finally:
        conn.close()
    assert set(found) == set(ids)
    for system_id in ids:
        assert found[system_id] == pytest.approx(stored[system_id], abs=1e-6)


def test_the_hops_of_a_cross_sector_route_are_the_distances_between_the_stored_positions(mysql_config):
    west = _save_grid_sector(mysql_config, "West", FAR_WEDGE, _cloud_local(1))
    east = _save_grid_sector(mysql_config, "East", FAR_WEDGE + 20, _cloud_local(2))
    stored = _all_systems_ly(mysql_config)
    route = _nav(mysql_config, west[0], east[-1])["route"]
    assert len(route["path"]) >= 2
    for node, position in route["positions"].items():
        assert position == pytest.approx(stored[node], abs=1e-6)
    for hop in route["hops"]:
        assert hop["distance_ly"] == pytest.approx(detmath.dist(stored[hop["from"]], stored[hop["to"]]), abs=1e-6)
    assert route["distance_ly"] == pytest.approx(sum(hop["distance_ly"] for hop in route["hops"]))


def test_a_route_inside_one_sector_that_may_leave_it_comes_back_in_the_sectors_own_axes(mysql_config):
    ids = _save_grid_sector(mysql_config, "Far", FAR_WEDGE, _cloud_local(4, count=12, spread=6.0))
    conn = store.get_connection(mysql_config)
    try:
        local = {row["id"]: tuple(milliparsecs_to_ly(row[axis]) for axis in ("position_x_mpc", "position_y_mpc", "position_z_mpc"))
                 for row in conn.execute("SELECT id, position_x_mpc, position_y_mpc, position_z_mpc FROM star_systems").fetchall()}
        center = conn.execute("SELECT center_x_pc, center_y_pc, center_z_pc FROM sectors").fetchone()
        center = (center["center_x_pc"], center["center_y_pc"], center["center_z_pc"])
        result = query.nav_between(conn, ids[0], ids[-1])
    finally:
        conn.close()
    route = result["route"]
    assert result["scope"] == "sector"
    for node, position in route["positions"].items():
        assert position == pytest.approx(local[node], abs=1e-6)
    # And the galaxy-frame position it was searched at turns into that local one.
    stored = _all_systems_ly(mysql_config)
    for node in route["path"]:
        turned = tuple(pc_to_ly(c) for c in galaxy_to_local_pc(center, tuple(ly_to_pc(v) for v in stored[node])))
        assert turned == pytest.approx(local[node], abs=1e-3)
