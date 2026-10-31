"""
GEN.315: a system built from the tree its scatter row stored (`planetgen.generation.stored_trees`,
`StarSystem(member_tree=...)`, docs/design/systems/multi-star-systems.md section 6). A scattered star with one companion star is
that pair, on the stored orbit; a single star draws no companion of its own; a member that is not a star (a brown dwarf, a
neutron star, a black hole) hangs in the tree dark (GEN.1013); a built row's members and nodes are deleted in the save.
"""

import dataclasses
import random

import pytest

from planetgen.db import orbit_trees, scatter_systems, store
from planetgen.generation import multistar, run_galaxy, run_plan, stored_trees
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.generation.system import StarSystem
from planetgen import tuning
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
        epoch, m.radius_au) for m in tree.members)
    return scatter_systems.StoredTree(members, orbit_trees.nodes_from_member_tree(tree, epoch))


PRIMARY_FE_H = 0.1
"""float: The [Fe/H] the placed star is built at (`_primary_params`) and its tree is drawn at, as a scatter row's is (GEN.651)."""


def _present_mass(mass, age_yr):
    """The mass a star born with `mass` has at `age_yr`, built at `PRIMARY_FE_H`: slightly under the birth mass on the
    main sequence once the star model loses mass (the mass the stored member must hold for the build to be that star)."""
    with draw.bound(0):
        star = Star(SystemConfig(), name="Placed", initial_mass_sol=mass, age_gy=age_yr / 1.0e9, metallicity_fe_h=PRIMARY_FE_H)
    return star.mass / SOLAR


def _tree_of(members, seed, mass=1.0, age_fraction=0.5):
    """The tree seed `seed` draws for a main-sequence star of `mass` half way through its life (or `None`)."""
    age_yr = evolution.main_sequence_lifetime_gy(mass) * 1.0e9 * age_fraction
    placed = multistar.Placed("star", _present_mass(mass, age_yr), mass, fe_h=PRIMARY_FE_H)
    tree = multistar.draw_member_tree(placed, age_yr, draw.Stream(f"{seed}:members"))
    return tree if tree.member_count == members else None, age_yr


def _tree(members, mass=1.0, start=1):
    for seed in range(start, start + 500):
        tree, age_yr = _tree_of(members, seed, mass)
        if tree is not None:
            return seed, tree, age_yr
    pytest.skip(f"no seed drew a tree of {members} members")


def _primary_params(mass, age_yr, seed):
    with draw.bound(seed):
        star = Star(SystemConfig(), name="Placed", initial_mass_sol=mass, age_gy=age_yr / 1.0e9, metallicity_fe_h=PRIMARY_FE_H)
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


def _pair_kind(node):
    """The kind of pair (`close` P-type, `wide` S-type) a node with this axis builds (GEN.631)."""
    return "close" if node.separation_km / constants.AU_TO_KM < tuning.WIDE_BINARY_SEPARATION_MIN_AU else "wide"


def _pair(kind, start=1):
    """A seed, tree and age of a pair that builds as a pair of this `kind`."""
    for seed in range(start, start + 4000):
        tree, age_yr = _tree_of(2, seed)
        if tree is not None and _pair_kind(tree.nodes[0]) == kind and tree.members[1].kind in stored_trees.PAIR_STAR_KINDS:
            return seed, tree, age_yr
    pytest.skip(f"no seed drew a {kind} pair")


# --- The plan -----------------------------------------------------------------------------------------------------


def test_a_row_without_members_is_a_single_star_that_draws_no_companion():
    assert stored_trees.plan(None) is stored_trees.SINGLE
    assert stored_trees.plan(scatter_systems.StoredTree((), ())) is stored_trees.SINGLE
    assert stored_trees.SINGLE.complete and stored_trees.SINGLE.members == () and stored_trees.SINGLE.nodes == ()


@pytest.mark.parametrize("kind", ["close", "wide"])
def test_two_stars_and_one_node_are_a_pair_plan(kind):
    _seed, tree, _age = _pair(kind)
    plan = stored_trees.plan(_stored(tree), EPOCH)
    assert plan.complete and [m.member_no for m in plan.members] == [0, 1] and plan.pair_node_no() == 0
    assert plan.nodes[0].node_no == 0 and _pair_kind(plan.nodes[0]) == kind


def test_three_stars_are_a_plan_numbered_from_zero_and_a_remnant_or_a_brown_dwarf_hang_in_it_dark():
    _seed, three, _age = _tree(3)
    plan = stored_trees.plan(_stored(three), EPOCH)
    assert plan.complete
    assert [n.node_no for n in plan.nodes] == [0, 1] and plan.nodes[0].parent_node_no is None
    assert plan.nodes[1].parent_node_no == 0 and {m.parent_node_no for m in plan.members} <= {0, 1}
    assert [m.member_no for m in plan.dark_members] == [m.member_no for m in plan.members if m.kind in stored_trees.UNBUILT_KINDS]
    assert plan.pair_node_no() in (0, 1, None) and (plan.pair_node_no() is None) == (len(plan.lit_members) < 2)
    pair = _stored(_pair("close")[1])
    for kind in stored_trees.UNBUILT_KINDS:
        members = (pair.members[0], pair.members[1]._replace(kind=kind))
        dark = stored_trees.plan(scatter_systems.StoredTree(members, pair.nodes), EPOCH)
        assert dark.complete and [m.member_no for m in dark.dark_members] == [1]
        assert dark.pair_node_no() is None and dark.pair_members() is None and len(dark.nodes) == 1
    mystery = (pair.members[0], pair.members[1]._replace(kind="mystery"))
    assert not stored_trees.plan(scatter_systems.StoredTree(mystery, pair.nodes), EPOCH).complete


def _planned(member_no, parent, slot, kind, mass):
    return stored_trees.PlannedStar(member_no, parent, slot, kind, mass, mass, member_no, 0.01, None)


def _node(node_no, parent, slot, mass):
    return OrbitNode(node_no=node_no, parent_node_no=parent, slot=slot, mass_solar=mass, separation_km=1.5e9, eccentricity=0.1,
                     inclination_deg=0.0, ascending_node_deg=0.0, arg_periapsis_deg=0.0, mean_anomaly_deg=0.0, epoch_unix=0.0,
                     period_years=10.0, kind="wide", interaction="detached", flags=0)


def test_the_pair_of_a_tree_with_a_dark_member_between_its_stars_is_the_node_above_both():
    # A (1.0) in slot 0 of the root; B (0.5) and a brown dwarf (0.05) in the node in slot 1: no node holds two stars.
    members = (_planned(0, 0, 0, "main-sequence", 1.0), _planned(1, 1, 0, "main-sequence", 0.5),
               _planned(2, 1, 1, "brown-dwarf", 0.05))
    nodes = (_node(0, None, 0, 1.55), _node(1, 0, 1, 0.55))
    plan = stored_trees.StoredPlan(members, nodes, True)
    assert plan.pair_node_no() == 0 and plan.pair_members() == (0, 1)
    assert [m.member_no for m in plan.lit_members] == [0, 1] and [m.member_no for m in plan.dark_members] == [2]
    swapped = (_planned(0, 1, 1, "main-sequence", 1.0), _planned(1, 0, 0, "main-sequence", 0.5),
               _planned(2, 1, 0, "brown-dwarf", 0.05))
    both = stored_trees.StoredPlan(swapped, (_node(0, None, 0, 1.55), _node(1, 0, 1, 1.05)), True)
    assert both.pair_node_no() == 0 and both.pair_members() == (1, 0)


def test_the_node_of_two_stars_is_the_pair_even_when_a_dark_member_hangs_above_it():
    members = (_planned(0, 1, 0, "main-sequence", 1.0), _planned(1, 1, 1, "main-sequence", 0.5),
               _planned(2, 0, 1, "neutron-star", 1.4))
    plan = stored_trees.StoredPlan(members, (_node(0, None, 0, 2.9), _node(1, 0, 0, 1.5)), True)
    assert plan.pair_node_no() == 1 and plan.pair_members() == (0, 1)


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
    # the same epoch, so the very same orbit, numbered from 0 and of the kind the pair builds as
    # (a P-type pair is tracked on a circle)
    # (its angle from the ascending node is the argument of periapsis plus the anomaly, GEN.709)
    first = stored.nodes[0]
    assert node == dataclasses.replace(
        first, node_no=0, kind=kind,
        **({"eccentricity": 0.0, "arg_periapsis_deg": 0.0,
            "mean_anomaly_deg": (first.mean_anomaly_deg + first.arg_periapsis_deg) % 360.0} if kind == "close" else
           {"eccentricity": min(first.eccentricity, tuning.WIDE_BINARY_ECCENTRICITY_MAX * (1.0 - 1.0e-9))}))
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


def test_a_tree_with_a_remnant_builds_the_placed_star_with_the_remnant_dark_in_its_tree():
    seed, tree, age_yr = _pair("wide")
    stored = _stored(tree)
    plan = stored_trees.plan(
        scatter_systems.StoredTree((stored.members[0], stored.members[1]._replace(kind="neutron-star")), stored.nodes), EPOCH)
    system = _build(seed, tree, age_yr, plan=plan, BINARY_SYSTEM=True)
    assert system.tree_built and system.binary_type is None and len(system.stars) == 1
    assert [m.kind for m in system.dark_members] == ["neutron-star"] and system.star_slots == ((0, 0),)
    assert len(system.orbit_nodes) == 1 and system.orbit_nodes[0].mass_solar == pytest.approx(
        sum(m.mass_solar for m in plan.members))
    assert system.star.a_crit_au is not None
    loaded = StarSystem.from_dict(system.to_dict())
    assert loaded.dark_members == system.dark_members and loaded.star_slots == system.star_slots


def test_a_tree_with_a_member_of_a_kind_nothing_builds_builds_the_placed_star_alone():
    seed, tree, age_yr = _pair("wide")
    stored = _stored(tree)
    plan = stored_trees.plan(
        scatter_systems.StoredTree((stored.members[0], stored.members[1]._replace(kind="mystery")), stored.nodes), EPOCH)
    system = _build(seed, tree, age_yr, plan=plan, BINARY_SYSTEM=True)
    assert not system.tree_built and system.binary_type is None and len(system.stars) == 1
    assert system.orbit_nodes == () and system.dark_members == ()


def test_a_companion_that_is_not_its_stored_member_leaves_the_primary_alone(monkeypatch):
    seed, tree, age_yr = _pair("close")
    monkeypatch.setattr(stored_trees, "RECORD_TOLERANCE", -1.0)    # nothing built is the stored member
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
    assert any(plan.members for plan in plans.values()), "the toy galaxy drew no pair to build"

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
            if plan.members:
                # a dark member (GEN.1013: a brown dwarf, a remnant) hangs in the tree and is no star row
                lit = plan.lit_members
                assert sorted(s["role"] for s in stars) == (
                    ["single"] if len(lit) == 1 else sorted(["primary", "secondary"] + ["member"] * (len(lit) - 2)))
                expected = sorted((m.mass_solar for m in lit), reverse=True)
                assert [s["mass_kg"] / SOLAR for s in stars] == pytest.approx(expected, rel=1e-6)
                assert orbit_trees.read_dark_members(conn, system_id) == plan.dark_members
                assert len(nodes) == len(plan.nodes)
                for node, planned in zip(sorted(nodes, key=lambda n: n["node_no"]), plan.nodes):
                    for field in ("separation_km", "eccentricity", "inclination_deg", "ascending_node_deg", "period_years",
                                  "interaction", "flags"):
                        if field == "eccentricity" and node[field] == 0.0 and planned.kind == "close":
                            continue  # a P-type pair is built, and tracked, on a circle
                        if (field == "eccentricity" and len(plan.members) == 2
                                    and planned.eccentricity > tuning.WIDE_BINARY_ECCENTRICITY_MAX):
                            # only the two-star build holds a pair under the wide limit; an N-star build keeps the node
                            assert node[field] == pytest.approx(tuning.WIDE_BINARY_ECCENTRICITY_MAX, rel=1e-6)
                            continue  # a wide pair is held under the eccentricity limit
                        assert node[field] == pytest.approx(getattr(planned, field))
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
    chosen = next((i for i, tree in trees.items() if stored_trees.plan(tree).members), None)
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
