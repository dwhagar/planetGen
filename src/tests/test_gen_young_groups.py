# tests/test_gen_young_groups.py

"""
GEN.327: the young group. A massive star younger than 0.3 Myr is, at the `young-group` share, born in a group of 3 to 5
stars that flies straight (a `group` node, offset and velocity on each member, a decay age of 1e4 to 4e4 years). A system
younger than the decay age is stored as the group; an older one is stored as its result, the lightest star ejected and the
others in a tighter hierarchy. The sector page says how long a group has left, and the audit finds none past its decay.
"""

import dataclasses
import math

import pytest

from planetgen import tuning
from planetgen.db import query, scatter_systems, store
from planetgen.generation import multistar, scatter_companions, stored_trees
from planetgen.generation.phenomena import system_types
from planetgen.physics import constants
from planetgen.physics.system_tree import OrbitNode, TreeMember, member_offsets_au
from planetgen.util import draw
from planetgen.web.scatter_members_view import members_view
from tests.test_galactic_motion import SCATTER_RADIUS_PC, _insert_scatter_point, _seed_scatter_grid
from tests.test_gen_system_motion import EPOCH

MASSIVE = multistar.Placed("star", 20.0, 20.0)


@pytest.fixture
def always_a_group(monkeypatch):
    """Every eligible star is in a group (the share is 1)."""
    entry = dataclasses.replace(system_types.SYSTEM_TYPE_RATES["young-group"], default=1.0)
    monkeypatch.setitem(system_types.SYSTEM_TYPE_RATES, "young-group", entry)


def _tree(seed, age_yr, placed=MASSIVE):
    return multistar.draw_member_tree(placed, age_yr, draw.Stream(f"group:{seed}"), groups=True)


def _is_group(tree):
    return any(node.kind == "group" for node in tree.nodes)


def test_only_a_young_massive_star_can_be_in_a_group_and_the_others_draw_nothing():
    fresh = draw.Stream("untouched").random()
    for placed, age_yr in ((multistar.Placed("star", 20.0, 20.0), 4.0e5), (multistar.Placed("star", 3.0, 3.0), 2.0e4),
                           (multistar.Placed("star", 8.0, 8.0), 2.0e4), (multistar.Placed("neutron-star", 1.4, 20.0), 2.0e4)):
        rng = draw.Stream("untouched")
        assert multistar.draw_young_group(placed, age_yr, rng) is None
        assert rng.random() == fresh


def test_a_tree_is_a_group_only_when_the_caller_asks_for_one(always_a_group):
    assert all(not _is_group(multistar.draw_member_tree(MASSIVE, 2.0e3, draw.Stream(seed))) for seed in range(40))
    assert all(_is_group(_tree(seed, 2.0e3)) for seed in range(40))


def test_the_group_share_is_the_registry_default():
    share = system_types.SYSTEM_TYPE_RATES["young-group"].default
    runs = 3000
    groups = sum(_is_group(_tree(seed, 2.0e3)) for seed in range(runs))      # under 1e4 years, so never decayed
    sigma = math.sqrt(runs * share * (1 - share))
    assert abs(groups - runs * share) < 4 * sigma


def test_a_group_is_one_node_of_flying_members_around_a_barycentre_at_rest(always_a_group):
    sizes = set()
    for seed in range(60):
        tree = _tree(seed, 2.0e3)
        assert [node.kind for node in tree.nodes] == ["group"] and tree.nodes[0].parent_node_no is None
        sizes.add(tree.member_count)
        assert tree.ejected == ()
        low, high = tuning.YOUNG_GROUP_DECAY_AGE_RANGE_YR
        assert tree.member_count in set(tuning.YOUNG_GROUP_SIZE_WEIGHTS)
        assert [m.slot for m in tree.members] == list(range(tree.member_count))
        assert all(m.parent_node_no == tree.nodes[0].node_no for m in tree.members)
        assert len({m.decay_age_yr for m in tree.members}) == 1 and low <= tree.members[0].decay_age_yr <= high
        assert tree.age_yr < tree.members[0].decay_age_yr
        mass = tree.system_mass_solar
        for axis in range(3):
            assert sum(m.mass_solar * m.offset_au[axis] for m in tree.members) / mass == pytest.approx(0.0, abs=1e-6)
            assert sum(m.mass_solar * m.velocity_kms[axis] for m in tree.members) / mass == pytest.approx(0.0, abs=1e-9)
        assert tree.nodes[0].period_years == pytest.approx(
            tree.nodes[0].separation_au * constants.AU_TO_KM / tuning.YOUNG_GROUP_VELOCITY_DISPERSION_KMS
            / constants.SECONDS_PER_YEAR)
    assert sizes == {3, 4, 5}


def test_a_group_member_flies_straight_from_the_epoch_the_tree_was_drawn_at(always_a_group):
    tree = _tree(1, 2.0e3)
    node = OrbitNode(**{f: getattr(tree.nodes[0], f) for f in ("node_no", "parent_node_no", "slot",
                                        "mass_solar", "separation_km", "eccentricity", "inclination_deg",
                                        "ascending_node_deg", "arg_periapsis_deg", "mean_anomaly_deg", "period_years",
                                        "kind")}, epoch_unix=1000.0)
    members = [TreeMember(m.member_no, m.parent_node_no, m.slot, m.mass_solar, m.offset_au,
                                          m.velocity_kms) for m in tree.members]
    year = constants.SECONDS_PER_YEAR
    now, later = member_offsets_au(members, [node], 1000.0), member_offsets_au(members, [node], 1000.0 + 100 * year)
    for member in tree.members:
        moved = [b - a for a, b in zip(now[member.member_no], later[member.member_no])]
        expected = [v * 100 * year / constants.AU_TO_KM for v in member.velocity_kms]
        assert moved == pytest.approx(expected, rel=1e-9, abs=1e-9)
        assert now[member.member_no] == pytest.approx(member.offset_au)


def test_a_group_older_than_its_decay_age_is_never_drawn_but_its_remnant_is(always_a_group):
    decayed = 0
    for seed in range(120):
        tree = _tree(seed, 2.0e5)           # past every decay age
        assert not _is_group(tree)
        assert all(m.offset_au is None and m.decay_age_yr is None for m in tree.members)
        assert len(tree.ejected) == 1
        thrown = tree.ejected[0]
        assert thrown.flags & multistar.MEMBER_RUNAWAY_DYNAMICAL and thrown.parent_node_no is None
        if tree.nodes:
            decayed += 1
            assert tree.members[0].mass_solar == pytest.approx(MASSIVE.mass_solar)
            assert max(node.separation_au for node in tree.nodes) <= tuning.YOUNG_GROUP_REMNANT_MAX_AXIS_AU
            assert tree.member_count + tree.dropped + 1 in set(tuning.YOUNG_GROUP_SIZE_WEIGHTS)
            # the lightest left, never the placed star
            assert thrown.mass_solar <= min(m.mass_solar for m in tree.members[1:])
    assert decayed > 100


def test_a_group_straddling_its_decay_age_is_a_group_or_a_remnant_by_that_age(always_a_group):
    seen = set()
    for seed in range(200):
        tree = _tree(seed, 2.5e4)
        if _is_group(tree):
            assert tree.members[0].decay_age_yr > 2.5e4
        else:
            assert tree.ejected
        seen.add(_is_group(tree))
    assert seen == {True, False}


def test_no_drawn_group_is_past_its_decay_age_at_any_age(always_a_group):
    for age_yr in (1.0e3, 1.0e4, 2.0e4, 4.0e4, 1.0e5, 2.9e5):
        for seed in range(80):
            tree = _tree(seed, age_yr)
            assert all(m.decay_age_yr is None or tree.age_yr < m.decay_age_yr for m in tree.members)


def test_a_bright_star_row_young_enough_for_a_group_draws_one_and_an_older_one_does_not(always_a_group):
    def row(age_gy):
        columns = scatter_companions._STAR_COLUMN
        values = [None] * (max(columns.values()) + 1)
        values[columns["seed"]] = 11
        values[columns["yerkes_class"]] = "V"
        values[columns["mass_kg"]] = 20.0 * constants.SOLAR_MASS_TO_KG
        values[columns["initial_mass_sol"]] = 20.0
        values[columns["metallicity_fe_h"]] = 0.0
        values[columns["age_gy"]] = age_gy
        return values
    young = scatter_companions.star_tree(row(2.0e-6))
    assert young is not None and _is_group(young)
    old = scatter_companions.star_tree(row(1.0e-3))
    assert old is None or not _is_group(old)


def test_a_stored_group_keeps_its_offsets_velocities_and_decay_age_and_a_build_does_not_take_it(always_a_group, mysql_config):
    _seed_scatter_grid(mysql_config)
    tree = _tree(5, 2.0e3)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            row_id = _insert_scatter_point(conn, "bright_stars", (SCATTER_RADIUS_PC, 0.0, 0.0))
            scatter_systems.store_member_tree(conn, "bright_stars", row_id, tree, EPOCH)
        members, nodes = scatter_systems.read_member_tree(conn, "bright_stars", row_id)
        assert [n.kind for n in nodes] == ["group"]
        for stored, drawn in zip(members, tree.members):
            assert stored.offset_au == pytest.approx(drawn.offset_au) and stored.velocity_kms == pytest.approx(drawn.velocity_kms)
        found = conn.execute("SELECT member_no, decay_age_yr, flags FROM scatter_members WHERE source_id = ?",
                             (row_id,)).fetchall()
        assert {r["decay_age_yr"] for r in found} == {tree.members[0].decay_age_yr}
        positions = scatter_systems.member_positions_pc(conn, "bright_stars", row_id, EPOCH)
        assert len(positions) == tree.member_count
        plan = stored_trees.plan(scatter_systems.read_stored_trees(conn, "bright_stars", [row_id])[row_id], EPOCH)
        assert not plan.complete                  # a system build holds pairs, not a group: the placed star is built alone
    finally:
        conn.close()


def test_the_audit_finds_a_group_only_when_its_star_has_outlived_the_decay_age(always_a_group, mysql_config):
    _seed_scatter_grid(mysql_config)
    tree = _tree(5, 2.0e3)
    decay_yr = tree.members[0].decay_age_yr
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            row_id = _insert_scatter_point(conn, "bright_stars", (SCATTER_RADIUS_PC, 0.0, 0.0))
            scatter_systems.store_member_tree(conn, "bright_stars", row_id, tree, EPOCH)
            conn.execute("UPDATE bright_stars SET age_gy = ? WHERE id = ?", (0.5 * decay_yr / 1.0e9, row_id))
        assert scatter_systems.groups_past_decay(conn) == []
        with conn:
            conn.execute("UPDATE bright_stars SET age_gy = ? WHERE id = ?", (1.5 * decay_yr / 1.0e9, row_id))
        assert scatter_systems.groups_past_decay(conn) == [row_id]
    finally:
        conn.close()


def test_the_page_says_how_long_a_group_has_left(always_a_group, mysql_config):
    _seed_scatter_grid(mysql_config)
    tree = _tree(5, 2.0e3)
    decay_yr = tree.members[0].decay_age_yr
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            row_id = _insert_scatter_point(conn, "bright_stars", (SCATTER_RADIUS_PC, 0.0, 0.0))
            scatter_systems.store_member_tree(conn, "bright_stars", row_id, tree, EPOCH)
            conn.execute("UPDATE bright_stars SET age_gy = ? WHERE id = ?", (0.25 * decay_yr / 1.0e9, row_id))
        detail = query.uncharted_system_detail(conn, row_id)
    finally:
        conn.close()
    view = members_view(detail)
    [orbit] = view["orbits"]
    assert orbit[1].startswith("Young group, ") and orbit[1].endswith(" left") and orbit[4] == "–"
    assert len(view["members"]) == tree.member_count
    detail["age_gy"] = 2.0 * decay_yr / 1.0e9
    assert members_view(detail)["orbits"][0][1] == "Young group, past its decay age"
