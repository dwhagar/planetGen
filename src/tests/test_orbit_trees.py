# tests/test_orbit_trees.py

"""
GEN.311 and the storage half of GEN.269: a built system's tree of orbits is stored in `orbit_nodes` and read back with
the system. A plain binary is one node made from the pair it was generated with (the old binary columns stay filled and
equal it), the stars say where they hang in the tree, and a triple and a seven-star tree stored under a system read back
equal. The pair node carries the `interaction` its Roche geometry gives.
"""

import types

import pytest

from planetgen.db import orbit_trees, store
from planetgen.generation import multistar
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants
from planetgen.physics.system_tree import OrbitNode
from planetgen.util import draw


def _system(binary, wide=None):
    cfg = SystemConfig()
    cfg.PLANETS = False
    cfg.BINARY_SYSTEM = binary
    cfg.WIDE_BINARY = wide
    return StarSystem(system_config=cfg), cfg


def _saved(mysql_config, binary, wide=None):
    for _ in range(40):
        system, cfg = _system(binary, wide)
        if (system.binary_type is not None) == binary:
            break
    else:
        pytest.fail("no system of the asked kind in 40 draws")
    conn = store.get_connection(mysql_config)
    try:
        system_id = store.insert_star_system(conn, system, cfg)
        conn.commit()
        return conn, system, system_id
    except BaseException:
        conn.close()
        raise


@pytest.mark.parametrize("wide", [False, True])
def test_a_binary_is_one_node_equal_to_its_old_columns(mysql_config, wide):
    conn, system, system_id = _saved(mysql_config, True, wide)
    try:
        nodes = orbit_trees.read_nodes(conn, system_id)
        assert len(nodes) == 1 and nodes == system.orbit_nodes
        node = nodes[0]
        row = conn.execute("SELECT * FROM star_systems WHERE id = ?", (system_id,)).fetchone()
        assert node.node_no == 0 and node.parent_node_no is None and node.kind == ("wide" if wide else "close")
        assert node.separation_km == pytest.approx(row["binary_separation_km"], rel=1e-12)
        assert node.eccentricity == row["binary_eccentricity"]
        assert node.inclination_deg == row["binary_mutual_orbital_inclination_deg"]
        assert node.ascending_node_deg == row["binary_mutual_orbital_ascending_node_deg"]
        assert node.mean_anomaly_deg == row["binary_mutual_orbital_phase_deg"]
        assert node.period_years == row["binary_mutual_orbital_period_years"]
        assert node.mass_solar == pytest.approx(sum(s.mass for s in system.stars) / constants.SOLAR_MASS_TO_KG, rel=1e-12)
        assert node.interaction in multistar.INTERACTIONS
        assert node.epoch_unix > 0
    finally:
        conn.close()


@pytest.mark.parametrize("wide", [False, True])
def test_each_star_hangs_in_the_root_node_at_its_stored_offset(mysql_config, wide):
    conn, system, system_id = _saved(mysql_config, True, wide)
    try:
        rows = conn.execute("SELECT role, parent_node_id, offset_x_au, offset_y_au, offset_z_au FROM stars"
                            " WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchall()
        assert [r["role"] for r in rows] == ["primary", "secondary"]
        assert [r["parent_node_id"] for r in rows] == [0, 0]
        row = conn.execute("SELECT * FROM star_systems WHERE id = ?", (system_id,)).fetchone()
        for r, which in zip(rows, ("primary", "secondary")):
            for axis in "xyz":
                assert r[f"offset_{axis}_au"] == pytest.approx(
                    row[f"binary_{which}_position_{axis}_km"] / constants.AU_TO_KM, rel=1e-12, abs=1e-9)
    finally:
        conn.close()


def test_a_single_star_has_no_nodes_and_sits_at_the_centre(mysql_config):
    conn, system, system_id = _saved(mysql_config, False)
    try:
        assert orbit_trees.read_nodes(conn, system_id) == () and system.orbit_nodes == ()
        row = conn.execute("SELECT parent_node_id, offset_x_au, offset_y_au, offset_z_au FROM stars"
                           " WHERE star_system_id = ?", (system_id,)).fetchone()
        assert row["parent_node_id"] is None
        assert (row["offset_x_au"], row["offset_y_au"], row["offset_z_au"]) == (0.0, 0.0, 0.0)
    finally:
        conn.close()


@pytest.mark.parametrize("wide", [False, True])
def test_a_loaded_binary_has_the_stored_node(mysql_config, wide):
    conn, system, system_id = _saved(mysql_config, True, wide)
    try:
        loaded = store.load_star_system(conn, system_id)
        assert loaded.orbit_nodes == system.orbit_nodes and len(loaded.orbit_nodes) == 1
    finally:
        conn.close()


def test_the_nodes_go_through_the_json_form(mysql_config):
    conn, system, _system_id = _saved(mysql_config, True, False)
    conn.close()
    again = StarSystem.from_dict(system.to_dict())
    assert again.orbit_nodes == system.orbit_nodes and len(again.orbit_nodes) == 1
    assert StarSystem.from_dict({**system.to_dict(), "schema_version": 3, "orbit_nodes": []}).orbit_nodes == ()


@pytest.mark.parametrize("stars", [3, 7])
def test_a_drawn_tree_stored_under_a_system_reads_back_equal(mysql_config, stars):
    draw.set_run_seed(311)
    placed = multistar.Placed("star", 1.0)
    tree = multistar.draw_member_tree(placed, 3.0e9, draw, n_stars=stars)
    assert tree.star_count == stars and len(tree.nodes) == stars - 1
    conn, _system, system_id = _saved(mysql_config, False)
    try:
        nodes = orbit_trees.nodes_from_member_tree(tree, 1_760_000_000)
        orbit_trees.write_nodes(conn, system_id, nodes)
        assert orbit_trees.read_nodes(conn, system_id) == nodes
        assert all(isinstance(node, OrbitNode) for node in nodes)
        assert {node.interaction for node in nodes} <= set(multistar.INTERACTIONS)
        orbit_trees.write_nodes(conn, system_id, nodes[:1])       # a tree replaces the one stored
        assert orbit_trees.read_nodes(conn, system_id) == nodes[:1]
        orbit_trees.write_nodes(conn, system_id, ())
        assert orbit_trees.read_nodes(conn, system_id) == ()
    finally:
        conn.close()


def test_a_tree_the_system_carries_is_stored_as_it_is(mysql_config):
    draw.set_run_seed(312)
    tree = multistar.draw_member_tree(multistar.Placed("star", 1.0), 3.0e9, draw, n_stars=3)
    system, cfg = _system(False)
    system.orbit_nodes = orbit_trees.nodes_from_member_tree(tree, 1_760_000_000)
    conn = store.get_connection(mysql_config)
    try:
        system_id = store.insert_star_system(conn, system, cfg)
        conn.commit()
        assert orbit_trees.read_nodes(conn, system_id) == system.orbit_nodes
        assert store.load_star_system(conn, system_id).orbit_nodes == system.orbit_nodes
    finally:
        conn.close()


def _star(yerkes, radius_au, mass_sol, initial=None):
    return types.SimpleNamespace(yerkes_class=yerkes, radius=radius_au * constants.AU_TO_KM,
                                 mass=mass_sol * constants.SOLAR_MASS_TO_KG, initial_mass_sol=initial)


@pytest.mark.parametrize("pair, separation_au, expected", [
    ((_star("V", 0.01, 1.0), _star("V", 0.01, 1.0)), 1.0, "detached"),
    ((_star("III", 0.4, 1.0), _star("V", 0.005, 1.0)), 1.0, "semi-detached"),
    ((_star("III", 0.4, 1.0), _star("III", 0.4, 1.0)), 1.0, "contact"),
    ((_star("III", 0.6, 1.0), _star("VII", 4e-5, 0.6)), 1.0, "accreting"),
    ((_star("V", 0.005, 1.0), _star("VII", 4e-5, 0.6)), 20.0, "detached"),
    ((_star("V", 0.005, 1.0), _star("VII", 4e-5, 0.6, initial=3.0)), 0.5, "post-common-envelope"),
])
def test_a_pair_has_the_interaction_its_geometry_gives(pair, separation_au, expected):
    assert orbit_trees.pair_interaction(*pair, separation_au, 0.0) == expected


def test_the_closest_approach_of_an_eccentric_pair_decides():
    pair = (_star("III", 0.4, 1.0), _star("V", 0.005, 1.0))
    assert orbit_trees.pair_interaction(*pair, 3.0, 0.0) == "detached"
    assert orbit_trees.pair_interaction(*pair, 3.0, 0.9) == "semi-detached"


# --- Members that are not stars (GEN.1013) ------------------------------------------------------------------------


def _with_a_dark_member(kind="neutron-star", pair_kind="wide"):
    """A star built from a stored pair whose companion is `kind`: a single star with the companion dark in its one node."""
    from planetgen.db import scatter_systems
    from planetgen.generation import stored_trees
    from tests.test_gen_stored_trees import EPOCH, _build, _pair, _stored
    seed, tree, age_yr = _pair(pair_kind)
    stored = _stored(tree)
    plan = stored_trees.plan(
        scatter_systems.StoredTree((stored.members[0], stored.members[1]._replace(kind=kind)), stored.nodes), EPOCH)
    system = _build(seed, tree, age_yr, plan=plan, BINARY_SYSTEM=True)
    return system, SystemConfig()


def test_a_lone_star_with_a_dark_member_stores_the_member_and_the_tree_and_reads_both_back(mysql_config):
    system, cfg = _with_a_dark_member()
    assert system.binary_type is None and len(system.dark_members) == 1
    conn = store.get_connection(mysql_config)
    try:
        system_id = store.insert_star_system(conn, system, cfg)
        conn.commit()
        assert orbit_trees.read_dark_members(conn, system_id) == system.dark_members
        nodes = orbit_trees.read_nodes(conn, system_id)
        assert len(nodes) == 1 and nodes[0].mass_solar == pytest.approx(
            system.star.mass / constants.SOLAR_MASS_TO_KG + system.dark_members[0].mass_solar, rel=1e-9)
        roles = [r["role"] for r in conn.execute("SELECT role FROM stars WHERE star_system_id = ?", (system_id,)).fetchall()]
        assert roles == ["single"]
        loaded = store.load_star_system(conn, system_id)
        assert loaded.dark_members == system.dark_members and loaded.star_slots == system.star_slots == ((0, 0),)
        assert len(loaded.stars) == 1 and loaded.binary_type is None
        orbit_trees.write_dark_members(conn, system_id, ())
        assert orbit_trees.read_dark_members(conn, system_id) == ()
    finally:
        conn.close()


def test_the_star_of_a_system_with_a_dark_member_is_off_the_barycentre_in_proportion_to_the_masses(mysql_config):
    system, _cfg = _with_a_dark_member()
    star_mass = system.star.mass / constants.SOLAR_MASS_TO_KG
    dark = system.dark_members[0]
    node = system.orbit_nodes[0]
    [offset] = orbit_trees.system_tree.star_offsets_au(
        [star_mass], system.star_slots, system.orbit_nodes, node.epoch_unix, [(dark.mass_solar, (dark.parent_node_no, dark.slot))])
    distance = sum(c * c for c in offset) ** 0.5
    assert distance == pytest.approx(node.separation_km / constants.AU_TO_KM * dark.mass_solar / (star_mass + dark.mass_solar),
                                     rel=0.6)    # (an eccentric orbit: the separation at the node's phase, not the axis)
    assert distance > 0.0
