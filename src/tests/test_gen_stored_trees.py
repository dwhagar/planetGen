"""
GEN.315: a system built from the tree its scatter row stored (`planetgen.generation.stored_trees`,
`StarSystem(member_tree=...)`, docs/design/multi-star-systems.md section 6). A scattered star with one companion star is
that pair, on the stored orbit; a single star draws no companion of its own; a tree the pair classes cannot hold builds
the placed object alone and keeps its rows; a built row's members and nodes are deleted in the save.
"""

import dataclasses
import random

import pytest

from planetgen import tuning
from planetgen.db import orbit_trees, scatter_systems, store
from planetgen.generation import bright_stars as brightStars
from planetgen.generation import multistar, run_galaxy, run_plan, stored_trees
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.generation.system import StarSystem
from planetgen.physics import constants, stellar_evolution as evolution
from planetgen.physics.system_tree import OrbitNode
from planetgen.util import draw
from tests.test_bright_star_scatter import EDGE_PC, THRESHOLD, _scatter_args, _seed_galaxy, light_steps  # noqa: F401

EPOCH = 1.7e9
SOLAR = constants.SOLAR_MASS_TO_KG


def _stored(tree, epoch=EPOCH):
    """A drawn tree as the rows `read_stored_trees` returns for it."""
    members = tuple(scatter_systems.StoredMember(
        m.member_no, m.parent_node_no, m.slot, m.kind, m.state, m.flags, m.mass_solar, m.initial_mass_solar, m.seed,
        epoch) for m in tree.members)
    return scatter_systems.StoredTree(members, orbit_trees.nodes_from_member_tree(tree, epoch))


def _tree_of(members, seed, mass=1.0, age_fraction=0.5):
    """The tree seed `seed` draws for a main-sequence star of `mass` half way through its life (or `None`)."""
    age_yr = evolution.main_sequence_lifetime_gy(mass) * 1.0e9 * age_fraction
    tree = multistar.draw_member_tree(multistar.Placed("star", mass, mass), age_yr, draw.Stream(f"{seed}:members"))
    return tree if tree.member_count == members else None, age_yr


def _tree(members, mass=1.0, start=1):
    for seed in range(start, start + 500):
        tree, age_yr = _tree_of(members, seed, mass)
        if tree is not None:
            return seed, tree, age_yr
    pytest.skip(f"no seed drew a tree of {members} members")


def _primary_params(mass, age_yr, seed):
    with draw.bound(seed):
        star = Star(SystemConfig(), name="Placed", initial_mass_sol=mass, age_gy=age_yr / 1.0e9)
    return {"type": star.type, "yerkes_class": star.yerkes_class, "mass_kg": star.mass, "radius_km": star.radius,
            "temperature_k": star.temperature, "luminosity_w": star.luminosity, "age_gy": star.age,
            "lifespan_gy": star.lifespan, "initial_mass_sol": star.initial_mass_sol,
            "phase_end_age_gy": getattr(star, "phase_end_age_gy", None), "metallicity_fe_h": star.metallicity_fe_h}


def _build(seed, tree, age_yr, mass=1.0, plan=None, **config):
    cfg = SystemConfig()
    cfg.POPULATION = "thin-disk"
    for name, value in config.items():
        setattr(cfg, name, value)
    plan = stored_trees.plan(_stored(tree), EPOCH) if plan is None else plan
    with draw.bound(seed):
        return StarSystem(cfg, primary_star_params=_primary_params(mass, age_yr, seed), member_tree=plan)


def _pair(kind, start=1):
    """A seed, tree and age of a pair whose node is `kind`."""
    for seed in range(start, start + 2000):
        tree, age_yr = _tree_of(2, seed)
        if tree is not None and tree.nodes[0].kind == kind and tree.members[1].kind in stored_trees.PAIR_STAR_KINDS:
            return seed, tree, age_yr
    pytest.skip(f"no seed drew a {kind} pair")


# --- The plan -----------------------------------------------------------------------------------------------------


def test_a_row_without_members_is_a_single_star_that_draws_no_companion():
    assert stored_trees.plan(None) is stored_trees.SINGLE
    assert stored_trees.plan(scatter_systems.StoredTree((), ())) is stored_trees.SINGLE
    assert stored_trees.SINGLE.complete and stored_trees.SINGLE.companion is None and stored_trees.SINGLE.node is None


@pytest.mark.parametrize("kind", ["close", "wide"])
def test_two_stars_and_one_node_are_a_pair_plan(kind):
    _seed, tree, _age = _pair(kind)
    plan = stored_trees.plan(_stored(tree), EPOCH)
    assert plan.complete and plan.companion.member_no == 1 and plan.node.kind == kind


def test_three_members_a_remnant_or_a_brown_dwarf_build_the_placed_star_alone():
    _seed, three, _age = _tree(3)
    plan = stored_trees.plan(_stored(three), EPOCH)
    assert not plan.complete and plan.companion is None and plan.node is None
    pair = _stored(_pair("close")[1])
    for kind in ("neutron-star", "black-hole", "brown-dwarf"):
        members = (pair.members[0], pair.members[1]._replace(kind=kind))
        assert not stored_trees.plan(scatter_systems.StoredTree(members, pair.nodes), EPOCH).complete


def test_a_compact_or_group_node_is_not_a_pair_the_classes_hold():
    pair = _stored(_pair("close")[1])
    for kind in ("compact", "group"):
        node = dataclasses.replace(pair.nodes[0], kind=kind)
        assert not stored_trees.plan(scatter_systems.StoredTree(pair.members, (node,)), EPOCH).complete


def test_a_node_is_carried_to_the_orbit_epoch_on_the_same_orbit():
    node = OrbitNode(node_no=1, parent_node_no=None, slot=0, mass_solar=1.5, separation_km=1.5e8, eccentricity=0.2,
                     inclination_deg=30.0, ascending_node_deg=10.0, arg_periapsis_deg=0.0, mean_anomaly_deg=350.0,
                     epoch_unix=100.0, period_years=2.0, kind="close", interaction="detached", flags=0)
    year = constants.SECONDS_PER_YEAR
    assert stored_trees.node_at(node, None) is node and stored_trees.node_at(node, 100.0) is node
    moved = stored_trees.node_at(node, 100.0 + 0.5 * year)
    assert moved.epoch_unix == 100.0 + 0.5 * year and moved.mean_anomaly_deg == pytest.approx(350.0 + 90.0 - 360.0)
    assert dataclasses.replace(moved, mean_anomaly_deg=node.mean_anomaly_deg, epoch_unix=node.epoch_unix) == node
    whole = stored_trees.node_at(node, 100.0 + 2.0 * year)
    assert whole.mean_anomaly_deg == pytest.approx(350.0)


# --- The system ---------------------------------------------------------------------------------------------------


@pytest.mark.parametrize("kind", ["close", "wide"])
def test_a_pair_is_the_stored_stars_on_the_stored_orbit(kind):
    seed, tree, age_yr = _pair(kind)
    system = _build(seed, tree, age_yr)
    stored = _stored(tree)
    assert system.tree_built and system.binary_type == kind and len(system.stars) == 2
    masses = sorted(star.mass / SOLAR for star in system.stars)
    assert masses == pytest.approx(sorted(m.mass_solar for m in tree.members), rel=1e-6)
    (node,) = system.orbit_nodes
    assert node == stored.nodes[0]    # the same epoch, so the very same node
    pair = system.star if kind == "close" else system.wide_binary
    if kind == "close":
        assert pair.binary_separation_au * constants.AU_TO_KM == pytest.approx(node.separation_km)
        assert (pair.binary_mutual_orbital_inclination_deg, pair.binary_mutual_orbital_ascending_node_deg,
                pair.binary_mutual_orbital_phase_deg, pair.binary_mutual_orbital_period_years) == (
            node.inclination_deg, node.ascending_node_deg, node.mean_anomaly_deg, node.period_years)
    else:
        assert pair.separation_au * constants.AU_TO_KM == pytest.approx(node.separation_km)
        assert pair.eccentricity == node.eccentricity
        assert (pair.inclination_deg, pair.ascending_node_deg, pair.phase_deg, pair.period_years) == (
            node.inclination_deg, node.ascending_node_deg, node.mean_anomaly_deg, node.period_years)


def test_a_pair_built_later_holds_the_same_orbit_at_the_later_time():
    seed, tree, age_yr = _pair("wide")
    stored = _stored(tree)
    later = EPOCH + 1.0e7
    system = _build(seed, tree, age_yr, plan=stored_trees.plan(stored, later))
    (node,) = system.orbit_nodes
    assert node.epoch_unix == later and node.separation_km == stored.nodes[0].separation_km
    assert system.wide_binary.phase_deg == node.mean_anomaly_deg


@pytest.mark.parametrize("kind", ["close", "wide"])
def test_the_second_star_is_built_from_its_own_seed_not_the_systems_stream(kind):
    seed, tree, age_yr = _pair(kind)
    first = _build(seed, tree, age_yr)
    other = _build(seed + 1, tree, age_yr)     # another stream for everything but the stored second star
    def second(system):
        return min(system.stars, key=lambda star: star.mass)
    assert (second(first).type, second(first).radius, second(first).temperature, second(first).metallicity_fe_h) == (
        second(other).type, second(other).radius, second(other).temperature, second(other).metallicity_fe_h)


def test_the_same_inputs_build_the_same_system():
    seed, tree, age_yr = _pair("close")
    first, second = _build(seed, tree, age_yr), _build(seed, tree, age_yr)
    assert first.to_dict() == second.to_dict()


@pytest.mark.parametrize("forced", [{"BINARY_SYSTEM": True}, {"BINARY_SYSTEM": True, "WIDE_BINARY": True},
                                    {"HABITABLE_WORLD": True}])
def test_a_single_star_row_draws_no_companion_whatever_the_config_asks(forced):
    _seed, tree, age_yr = _tree(1)
    for seed in range(5):
        system = _build(seed, tree, age_yr, plan=stored_trees.SINGLE, **forced)
        assert system.binary_type is None and len(system.stars) == 1 and system.orbit_nodes == ()


def test_the_stored_kind_stands_when_a_habitable_world_is_asked_of_a_close_pair():
    seed, tree, age_yr = _pair("close")
    system = _build(seed, tree, age_yr, HABITABLE_WORLD=True)
    assert system.binary_type == "close"


def test_a_tree_the_pair_classes_cannot_hold_builds_the_placed_star_alone_and_keeps_its_rows():
    seed, tree, age_yr = _tree(3)
    plan = stored_trees.plan(_stored(tree), EPOCH)
    system = _build(seed, tree, age_yr, plan=plan, BINARY_SYSTEM=True)
    assert not system.tree_built and system.binary_type is None and len(system.stars) == 1
    assert system.orbit_nodes == ()


def test_a_companion_that_is_not_its_stored_member_leaves_the_primary_alone(monkeypatch):
    seed, tree, age_yr = _pair("close")
    monkeypatch.setattr(stored_trees, "MASS_TOLERANCE", -1.0)    # nothing built is the stored mass
    system = _build(seed, tree, age_yr)
    assert not system.tree_built and system.binary_type is None and system.orbit_nodes == ()


def test_a_system_generated_whole_is_not_a_built_tree():
    with draw.bound(5):
        system = StarSystem(SystemConfig())
    assert not system.tree_built


# --- The fill -----------------------------------------------------------------------------------------------------


def test_a_fill_builds_each_row_from_its_stored_tree_and_spends_the_rows_of_the_trees_it_built(mysql_config):
    _seed_galaxy(mysql_config)
    run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass"))
    conn = store.get_connection(mysql_config)
    try:
        target = conn.execute(
            "SELECT shell_index, band_index, wedge_index, COUNT(*) AS n FROM bright_stars WHERE member_count > 1"
            " GROUP BY shell_index, band_index, wedge_index ORDER BY n DESC LIMIT 1").fetchone()
        address = (target["shell_index"], target["band_index"], target["wedge_index"])
        rows = store.bright_stars_for_sector(conn, *address)
        trees = scatter_systems.read_stored_trees(conn, "bright_stars", [r["id"] for r in rows])
        epoch = store.get_orbit_epoch_unix(conn)
    finally:
        conn.close()
    plans = {row["id"]: stored_trees.plan(trees.get(row["id"]), epoch) for row in rows}
    assert any(plan.companion is not None for plan in plans.values()), "the toy galaxy drew no pair to build"

    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 0
    from planetgen.galaxy.spherical import sector_position_pc
    run_galaxy.generate_and_save_sector_at(args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)

    conn = store.get_connection(mysql_config)
    try:
        for row in rows:
            plan, system_id = plans[row["id"]], conn.execute(
                "SELECT star_system_id FROM bright_stars WHERE id = ?", (row["id"],)).fetchone()["star_system_id"]
            assert system_id is not None
            stars = conn.execute("SELECT role, mass_kg FROM stars WHERE star_system_id = ? ORDER BY mass_kg DESC",
                                 (system_id,)).fetchall()
            nodes = conn.execute("SELECT * FROM orbit_nodes WHERE star_system_id = ?", (system_id,)).fetchall()
            left = scatter_systems.read_stored_trees(conn, "bright_stars", [row["id"]])
            if plan.companion is not None:
                assert [s["role"] for s in stars] == ["primary", "secondary"]
                expected = sorted((m.mass_solar for m in trees[row["id"]].members), reverse=True)
                assert [s["mass_kg"] / SOLAR for s in stars] == pytest.approx(expected, rel=1e-6)
                (node,) = nodes
                for field in ("separation_km", "eccentricity", "inclination_deg", "ascending_node_deg", "period_years",
                              "kind", "interaction", "flags"):
                    assert node[field] == pytest.approx(getattr(plan.node, field))
                assert not left
            elif plan.complete:
                assert [s["role"] for s in stars] == ["single"] and not nodes and not left
            else:
                assert [s["role"] for s in stars] == ["single"] and not nodes
                assert left, "the rows of a tree that was not built stay for the build that can"
    finally:
        conn.close()


def test_the_uncharted_route_builds_the_stored_tree_and_spends_it(mysql_config):
    from planetgen.api import routes
    _seed_galaxy(mysql_config)
    run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass"))
    conn = store.get_connection(mysql_config)
    try:
        rows = [r for r in conn.execute("SELECT id FROM bright_stars WHERE member_count = 2").fetchall()]
        trees = scatter_systems.read_stored_trees(conn, "bright_stars", [r["id"] for r in rows])
    finally:
        conn.close()
    chosen = next((i for i, tree in trees.items() if stored_trees.plan(tree).companion is not None), None)
    assert chosen is not None, "the toy galaxy drew no star with one companion star"
    routes.generate_uncharted_system_job(mysql_config, chosen)
    conn = store.get_connection(mysql_config)
    try:
        system_id = conn.execute("SELECT star_system_id FROM bright_stars WHERE id = ?", (chosen,)).fetchone()[
            "star_system_id"]
        assert conn.execute("SELECT COUNT(*) AS n FROM stars WHERE star_system_id = ?", (system_id,)).fetchone()["n"] == 2
        assert conn.execute("SELECT COUNT(*) AS n FROM orbit_nodes WHERE star_system_id = ?",
                            (system_id,)).fetchone()["n"] == 1
        assert not scatter_systems.read_stored_trees(conn, "bright_stars", [chosen])
    finally:
        conn.close()


# --- The reader ---------------------------------------------------------------------------------------------------


def test_stored_trees_are_read_with_every_column_and_deleted_by_row(mysql_config):
    from planetgen.generation import scatter_companions
    seeds = random.Random(3).sample(range(1, 10 ** 9), 40)
    rows = [row for row in (_row for _row in _star_rows(seeds)) if scatter_companions.star_tree(row) is not None]
    assert len(rows) > 5
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            scatter_companions.insert_bright_stars(conn, rows, EPOCH)
        ids = {r["seed"]: r["id"] for r in conn.execute("SELECT id, seed FROM bright_stars").fetchall()}
        read = scatter_systems.read_stored_trees(conn, "bright_stars", ids.values())
        assert set(read) == set(ids.values())
        for row in rows:
            tree = scatter_companions.star_tree(row)
            stored = read[ids[row[scatter_companions._STAR_COLUMN["seed"]]]]
            assert [(m.member_no, m.kind, m.state, m.flags, m.seed) for m in stored.members] == [
                (m.member_no, m.kind, m.state, m.flags, m.seed) for m in tree.members]
            assert [m.initial_mass_solar for m in stored.members] == [m.initial_mass_solar for m in tree.members]
            assert [(n.node_no, n.kind, n.interaction, n.flags) for n in stored.nodes] == [
                (n.node_no, n.kind, n.interaction, n.flags) for n in tree.nodes]
        keep, gone = list(ids.values())[:2], list(ids.values())[2:]
        with conn:
            scatter_systems.delete_member_trees(conn, "bright_stars", gone)
        assert set(scatter_systems.read_stored_trees(conn, "bright_stars", ids.values())) == set(keep)
        assert scatter_systems.read_stored_trees(conn, "bright_stars", []) == {}
    finally:
        conn.close()


def _star_rows(seeds):
    from tests.test_gen_scatter_companions import _star_row
    for i, seed in enumerate(seeds):
        yield _star_row(seed, initial_mass=0.8 + (i % 5) * 0.1, age_gy=1.0 + 0.1 * i)
