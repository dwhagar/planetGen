# tests/test_gen_system_motion.py

"""
GEN.265: the nightly motion of a scattered multi-star system. The barycentre is the scatter point and moves with the
rest of the galaxy; the members are computed on read from the stored node phases, so the update writes no member or
node row (`planetgen/physics/system_tree.py`, `planetgen/db/scatter_systems.py`,
docs/design/systems/multi-star-systems.md section 7).
"""

import math

import pytest

from planetgen.db import scatter_systems, store
from planetgen.generation import multistar
from planetgen.physics import constants
from planetgen.physics.kepler import solve_eccentric_anomaly
from planetgen.physics.system_tree import OrbitNode, TreeMember, member_offsets_au, node_relative_au
from planetgen.util import draw
from tests.test_galactic_motion import SCATTER_RADIUS_PC, _insert_scatter_point, _seed_scatter_grid

EPOCH = 1_700_000_000.0
YEAR = constants.SECONDS_PER_YEAR
AU_KM = constants.AU_TO_KM


def _node(node_no, parent, slot, mass, axis_au, e=0.0, period=10.0, anomaly=0.0, i=0.0, node=0.0, peri=0.0, kind="close"):
    return OrbitNode(node_no, parent, slot, mass, axis_au * AU_KM, e, i, node, peri, anomaly, EPOCH, period, kind)


def _length(vector):
    return math.sqrt(sum(c * c for c in vector))


def _sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def _centre(offsets, masses):
    total = sum(masses.values())
    return tuple(sum(masses[n] * offsets[n][axis] for n in offsets) / total for axis in range(3))


def test_a_circular_pair_swings_round_its_barycentre_and_sits_opposite_after_half_a_period():
    members = [TreeMember(0, 1, 0, 1.0), TreeMember(1, 1, 1, 0.25)]
    nodes = [_node(1, None, 0, 1.25, 10.0)]
    start = member_offsets_au(members, nodes, EPOCH)
    later = member_offsets_au(members, nodes, EPOCH + 5.0 * YEAR)
    assert _length(_sub(start[1], start[0])) == pytest.approx(10.0)
    assert _length(start[0]) == pytest.approx(10.0 * 0.25 / 1.25) and _length(start[1]) == pytest.approx(10.0 * 1.0 / 1.25)
    assert later[1] == pytest.approx(tuple(-c for c in start[1]), abs=1e-9)
    assert _centre(start, {0: 1.0, 1: 0.25}) == pytest.approx((0.0, 0.0, 0.0), abs=1e-12)


def test_an_eccentric_pair_follows_keplers_equation_from_its_stored_phase():
    e, axis, period, m0 = 0.5, 40.0, 100.0, 30.0
    node = _node(1, None, 0, 2.0, axis, e=e, period=period, anomaly=m0)
    for years in (0.0, 12.5, 25.0, 60.0, 99.0):
        anomaly = math.radians((m0 + 360.0 * years / period) % 360.0)
        expected = axis * (1.0 - e * math.cos(solve_eccentric_anomaly(anomaly, e)))
        assert _length(node_relative_au(node, EPOCH + years * YEAR)) == pytest.approx(expected, rel=1e-9)
    assert node_relative_au(node, EPOCH + period * YEAR) == pytest.approx(node_relative_au(node, EPOCH), abs=1e-9)


def test_phase_is_read_at_the_requested_time_before_the_epoch_as_well_as_after():
    node = _node(1, None, 0, 2.0, 5.0, e=0.3, period=7.0, anomaly=100.0)
    assert node_relative_au(node, EPOCH - 7.0 * YEAR) == pytest.approx(node_relative_au(node, EPOCH), abs=1e-9)
    assert node_relative_au(node, EPOCH - 3.0 * YEAR) != pytest.approx(node_relative_au(node, EPOCH), abs=1e-3)


def test_a_hierarchical_triple_keeps_its_barycentre_and_its_inner_pair_orbits_the_outer_centre():
    # members 0 and 1 are the close pair (node 2) that is slot 0 of the wide node 1; member 2 is its slot 1.
    members = [TreeMember(0, 2, 0, 1.0), TreeMember(1, 2, 1, 0.5), TreeMember(2, 1, 1, 0.8)]
    nodes = [_node(1, None, 0, 2.3, 300.0, e=0.4, period=2000.0, anomaly=70.0, i=20.0, node=40.0, peri=10.0, kind="wide"),
             _node(2, 1, 0, 1.5, 2.0, e=0.2, period=2.0, anomaly=10.0, i=5.0)]
    masses = {0: 1.0, 1: 0.5, 2: 0.8}
    for years in (0.0, 0.37, 400.0, 1234.5):
        offsets = member_offsets_au(members, nodes, EPOCH + years * YEAR)
        assert _centre(offsets, masses) == pytest.approx((0.0, 0.0, 0.0), abs=1e-9)
        assert _sub(offsets[1], offsets[0]) == pytest.approx(node_relative_au(nodes[1], EPOCH + years * YEAR), abs=1e-9)
        inner_centre = tuple((1.0 * offsets[0][a] + 0.5 * offsets[1][a]) / 1.5 for a in range(3))
        assert _sub(offsets[2], inner_centre) == pytest.approx(node_relative_au(nodes[0], EPOCH + years * YEAR), abs=1e-9)


def test_a_single_object_sits_at_the_barycentre():
    assert member_offsets_au([TreeMember(0, 0, 0, 1.0)], [], EPOCH) == {0: (0.0, 0.0, 0.0)}


def test_the_members_of_a_group_node_fly_straight_from_the_nodes_epoch():
    members = [TreeMember(0, 1, 0, 1.0, (1.0, 0.0, 0.0), (0.0, 30.0, 0.0)),
               TreeMember(1, 1, 1, 1.0, (-1.0, 0.0, 0.0), (0.0, -30.0, 0.0))]
    nodes = [_node(1, None, 0, 2.0, 1.0, kind="group")]
    offsets = member_offsets_au(members, nodes, EPOCH + 1000.0)
    assert offsets[0] == pytest.approx((1.0, 30.0 * 1000.0 / AU_KM, 0.0)) and offsets[1] == pytest.approx((-1.0, -30.0 * 1000.0 / AU_KM, 0.0))


def test_a_malformed_tree_is_refused():
    with pytest.raises(ValueError, match="slot 1"):
        member_offsets_au([TreeMember(0, 1, 0, 1.0)], [_node(1, None, 0, 1.0, 1.0)], EPOCH)
    with pytest.raises(ValueError, match="hang from no node"):
        member_offsets_au([TreeMember(0, 1, 0, 1.0), TreeMember(1, 1, 1, 1.0), TreeMember(2, 9, 0, 1.0)],
                          [_node(1, None, 0, 2.0, 1.0)], EPOCH)


@pytest.mark.parametrize("seed", range(8))
def test_drawn_trees_place_every_member_with_the_barycentre_at_the_centre_of_mass(seed):
    tree = multistar.draw_member_tree(multistar.Placed("star", 1.0), 5e9, draw.Stream(seed), n_stars=4)
    members = [TreeMember(m.member_no, m.parent_node_no, m.slot, m.mass_solar) for m in tree.members]
    nodes = [OrbitNode(n.node_no, n.parent_node_no, n.slot, n.mass_solar, n.separation_km, n.eccentricity, n.inclination_deg,
                       n.ascending_node_deg, n.arg_periapsis_deg, n.mean_anomaly_deg, EPOCH, n.period_years, n.kind)
             for n in tree.nodes]
    masses = {m.member_no: m.mass_solar for m in tree.members}
    for years in (0.0, 3.3, 777.0):
        offsets = member_offsets_au(members, nodes, EPOCH + years * YEAR)
        assert set(offsets) == set(masses)
        scale = max(_length(o) for o in offsets.values()) or 1.0
        assert _centre(offsets, masses) == pytest.approx((0.0, 0.0, 0.0), abs=1e-9 * scale)


def _store_tree(conn, table, row_id, seed=3):
    tree = multistar.draw_member_tree(multistar.Placed("star", 1.0), 5e9, draw.Stream(seed), n_stars=3)
    scatter_systems.store_member_tree(conn, table, row_id, tree, EPOCH)
    return tree


def _checksums(conn):
    return {name: conn.execute(f"CHECKSUM TABLE {name}").fetchone()["Checksum"] for name in ("scatter_members", "scatter_nodes")}


def test_the_nightly_update_moves_the_barycentre_and_writes_no_member_or_node_row(mysql_config):
    _seed_scatter_grid(mysql_config)
    point = (SCATTER_RADIUS_PC, 0.0, 0.0)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            ids = {table: _insert_scatter_point(conn, table, point, None if table == "bright_stars" else "black-hole")
                   for table in scatter_systems.SOURCES}
            for table, row_id in ids.items():
                _store_tree(conn, table, row_id)
            before_sums = _checksums(conn)
            before = {table: conn.execute(f"SELECT position_x_mpc, position_y_mpc, epoch_unix FROM {table} WHERE id = ?",
                                          (row_id,)).fetchone() for table, row_id in ids.items()}
            motion = store.advance_scatter_points(conn, store.orbit_clock(conn, 1e6))
            after = {table: conn.execute(f"SELECT position_x_mpc, position_y_mpc, epoch_unix FROM {table} WHERE id = ?",
                                         (row_id,)).fetchone() for table, row_id in ids.items()}
            after_sums = _checksums(conn)
            counts = {name: conn.execute(f"SELECT COUNT(*) AS n FROM {name}").fetchone()["n"] for name in ("scatter_members", "scatter_nodes")}
    finally:
        conn.close()
    assert motion["counts"] == {"bright_stars": 1, "phenomenon_scatter": 1}
    for table in ids:
        assert (after[table]["position_x_mpc"], after[table]["position_y_mpc"]) != (before[table]["position_x_mpc"], before[table]["position_y_mpc"])
        assert after[table]["epoch_unix"] != before[table]["epoch_unix"]
    assert after_sums == before_sums
    assert counts == {"scatter_members": 2 * 3, "scatter_nodes": 2 * 2}


def test_a_members_position_is_the_moved_barycentre_plus_its_offset_at_the_stored_phase(mysql_config):
    _seed_scatter_grid(mysql_config)
    point = (SCATTER_RADIUS_PC, 0.0, 0.0)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            row_id = _insert_scatter_point(conn, "bright_stars", point)
            tree = _store_tree(conn, "bright_stars", row_id)
            store.advance_scatter_points(conn, store.orbit_clock(conn, 1e6))
            row = conn.execute("SELECT position_x_mpc, position_y_mpc, position_z_mpc, epoch_unix FROM bright_stars WHERE id = ?",
                               (row_id,)).fetchone()
            when = EPOCH + 123.0 * YEAR
            positions = scatter_systems.member_positions_pc(conn, "bright_stars", row_id, when)
    finally:
        conn.close()
    members = [TreeMember(m.member_no, m.parent_node_no, m.slot, m.mass_solar) for m in tree.members]
    nodes = [OrbitNode(n.node_no, n.parent_node_no, n.slot, n.mass_solar, n.separation_km, n.eccentricity, n.inclination_deg,
                       n.ascending_node_deg, n.arg_periapsis_deg, n.mean_anomaly_deg, EPOCH, n.period_years, n.kind)
             for n in tree.nodes]
    offsets = member_offsets_au(members, nodes, when)
    centre = (row["position_x_mpc"] / 1000.0, row["position_y_mpc"] / 1000.0, row["position_z_mpc"] / 1000.0)
    assert set(positions) == {0, 1, 2}
    for number, offset in offsets.items():
        expected = tuple(centre[a] + offset[a] / constants.AU_PER_PARSEC for a in range(3))
        assert positions[number] == pytest.approx(expected, abs=1e-12)
    # the widest pair of this tree is its root node's orbit, and members never stray further than that from the centre
    widest = max(_length(_sub(positions[a], positions[b])) * constants.AU_PER_PARSEC for a in positions for b in positions)
    assert widest <= tree.nodes[0].separation_au * (1.0 + tree.nodes[0].eccentricity) * (1.0 + 1e-9)


def test_a_stored_tree_and_the_point_columns_round_trip_exactly(mysql_config):
    _seed_scatter_grid(mysql_config)
    point = (SCATTER_RADIUS_PC + 0.123456789012345, -7.0000000000001, 0.5)
    velocity = (12.3456789012345, -200.00000000001, 1.5)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            row_id = _insert_scatter_point(conn, "phenomenon_scatter", point, "black-hole", velocity)
            tree = _store_tree(conn, "phenomenon_scatter", row_id, seed=5)
            row = conn.execute("SELECT position_x_mpc, position_y_mpc, position_z_mpc, velocity_x_kms, velocity_y_kms,"
                               " velocity_z_kms, member_count, system_mass_solar FROM phenomenon_scatter WHERE id = ?",
                               (row_id,)).fetchone()
            members, nodes = scatter_systems.read_member_tree(conn, "phenomenon_scatter", row_id)
    finally:
        conn.close()
    assert (row["position_x_mpc"], row["position_y_mpc"], row["position_z_mpc"]) == tuple(c * 1000.0 for c in point)
    assert (row["velocity_x_kms"], row["velocity_y_kms"], row["velocity_z_kms"]) == velocity
    assert row["member_count"] == tree.member_count and row["system_mass_solar"] == tree.system_mass_solar
    assert [(m.member_no, m.parent_node_no, m.slot, m.mass_solar) for m in members] == \
        [(m.member_no, m.parent_node_no, m.slot, m.mass_solar) for m in tree.members]
    assert [(n.node_no, n.parent_node_no, n.slot, n.mass_solar, n.separation_km, n.eccentricity, n.inclination_deg,
             n.ascending_node_deg, n.arg_periapsis_deg, n.mean_anomaly_deg, n.epoch_unix, n.period_years, n.kind)
            for n in nodes] == \
        [(n.node_no, n.parent_node_no, n.slot, n.mass_solar, n.separation_km, n.eccentricity, n.inclination_deg,
          n.ascending_node_deg, n.arg_periapsis_deg, n.mean_anomaly_deg, EPOCH, n.period_years, n.kind) for n in tree.nodes]


def test_storing_again_replaces_the_tree_and_a_single_object_clears_it(mysql_config):
    _seed_scatter_grid(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            row_id = _insert_scatter_point(conn, "bright_stars", (SCATTER_RADIUS_PC, 0.0, 0.0))
            _store_tree(conn, "bright_stars", row_id, seed=1)
            _store_tree(conn, "bright_stars", row_id, seed=2)
            members, nodes = scatter_systems.read_member_tree(conn, "bright_stars", row_id)
            single = multistar.draw_member_tree(multistar.Placed("star", 1.0), 5e9, draw.Stream(1), n_stars=1)
            scatter_systems.store_member_tree(conn, "bright_stars", row_id, single, EPOCH)
            cleared = scatter_systems.read_member_tree(conn, "bright_stars", row_id)
            row = conn.execute("SELECT member_count, system_mass_solar FROM bright_stars WHERE id = ?", (row_id,)).fetchone()
            alone = scatter_systems.member_positions_pc(conn, "bright_stars", row_id, EPOCH)
    finally:
        conn.close()
    assert (len(members), len(nodes)) == (3, 2)
    assert cleared == ([], []) and row["member_count"] is None and row["system_mass_solar"] is None
    assert alone == {0: (SCATTER_RADIUS_PC, 0.0, 0.0)}


def test_a_point_carries_the_systems_mass_where_it_has_members_for_the_influence_rule(mysql_config):
    _seed_scatter_grid(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            star = _insert_scatter_point(conn, "bright_stars", (SCATTER_RADIUS_PC, 0.0, 0.0))
            pair = _insert_scatter_point(conn, "bright_stars", (SCATTER_RADIUS_PC, 10.0, 0.0))
            hole = _insert_scatter_point(conn, "phenomenon_scatter", (SCATTER_RADIUS_PC, 20.0, 0.0), "black-hole")
            conn.execute("UPDATE phenomenon_scatter SET mass_solar = 9.0 WHERE id = ?", (hole,))
            tree = _store_tree(conn, "bright_stars", pair)
            stars = scatter_systems.point_masses_solar(conn, "bright_stars", [star, pair])
            holes = scatter_systems.point_masses_solar(conn, "phenomenon_scatter", [hole])
            none = scatter_systems.point_masses_solar(conn, "bright_stars", [])
    finally:
        conn.close()
    assert stars[star] == pytest.approx(1e31 / constants.SOLAR_MASS_TO_KG)
    assert stars[pair] == tree.system_mass_solar and holes == {hole: 9.0} and none == {}


def test_a_table_that_is_not_a_scatter_is_refused(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        with pytest.raises(ValueError, match="not a scatter table"):
            scatter_systems.read_member_tree(conn, "stars", 1)
        with pytest.raises(LookupError):
            scatter_systems.member_positions_pc(conn, "bright_stars", 987654321, EPOCH)
    finally:
        conn.close()
