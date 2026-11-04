# tests/test_gen_tree_systems.py

"""
GEN.631: a system the fill makes from nothing draws its stars as a tree (`multistar.draw_member_tree`) at the odds of the
`star-count` table, and builds them: the pair of the tree is the system's primary and secondary, every other star one of
its `member_stars`. The old one-companion draw and its probability table are gone.
"""

import collections
import dataclasses

import pytest

from planetgen import tuning
from planetgen.db import orbit_trees, store
from planetgen.generation import multistar, place_terms, stored_trees
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena import hosts
from planetgen.generation.system import StarSystem
from planetgen.physics import constants, system_tree
from planetgen.util import draw


def _system(seed, **config):
    cfg = SystemConfig()
    cfg.PLANETS = False
    for name, value in config.items():
        setattr(cfg, name, value)
    with draw.bound(seed):
        return StarSystem(system_config=cfg), cfg


def _multiple(count=3):
    for seed in range(2000):
        system, cfg = _system(seed)
        if len(system.stars) >= count:
            return system, cfg
    pytest.fail(f"no system of {count} stars in 2000 draws")


def test_the_old_companion_draw_and_its_table_are_gone():
    assert not hasattr(StarSystem, "_companion") and not hasattr(StarSystem, "_stored_companion")
    assert not hasattr(tuning, "BINARY_SYSTEM_PROBABILITY_BY_SPECTRAL_CLASS")


def test_the_stars_per_system_follow_the_star_count_table(monkeypatch):
    seen = {}
    drawn_plan = StarSystem._drawn_plan

    def watch(self):
        seen["class"] = multistar.companion_class(place_terms.born_mass_solar(self.primary_star))
        return drawn_plan(self)

    monkeypatch.setattr(StarSystem, "_drawn_plan", watch)
    counts = collections.defaultdict(collections.Counter)
    for seed in range(4000):
        system, _cfg = _system(seed)
        counts[seen["class"]][len(system.stars)] += 1
    checked = 0
    for letter, tally in counts.items():
        total = sum(tally.values())
        if total < 400 or letter == "M":    # a brown dwarf is no host class, so M dwarfs lose companions (GEN.631)
            continue
        weights = hosts.STAR_COUNT_PROBABILITIES[letter]
        single = weights[0] / sum(weights)
        assert tally[1] / total == pytest.approx(single, abs=4 * (single * (1 - single) / total) ** 0.5), letter
        checked += 1
    assert checked >= 1


def test_every_star_of_the_tree_is_built_and_hangs_in_it():
    system, _cfg = _multiple(3)
    assert len(system.stars) == 2 + len(system.member_stars) and system.member_stars
    assert system.stars[:2] == [system.primary_star, system.secondary_star]
    assert system.primary_star.mass >= system.secondary_star.mass
    assert len(system.star_slots) == len(system.stars) and len(system.orbit_nodes) == len(system.stars) - 1
    assert system.binary_type in ("close", "wide")
    children = collections.Counter(slot for slot in system.star_slots)
    children.update((node.parent_node_no, node.slot) for node in system.orbit_nodes if node.parent_node_no is not None)
    assert set(children.values()) == {1}    # every slot of every node is held by exactly one star or node
    assert sorted({node.parent_node_no for node in system.orbit_nodes if node.parent_node_no is not None}) \
        <= [node.node_no for node in system.orbit_nodes]
    assert len({star.name for star in system.stars}) == len(system.stars)
    total = sum(star.mass for star in system.stars) / constants.SOLAR_MASS_TO_KG
    root = next(node for node in system.orbit_nodes if node.parent_node_no is None)
    assert root.mass_solar == pytest.approx(total, rel=1e-9)
    for node in system.orbit_nodes:    # the period of each axis about the mass the stars came out with
        axis_au = node.separation_km / constants.AU_TO_KM
        assert node.period_years == pytest.approx((axis_au ** 3 / node.mass_solar) ** 0.5, rel=1e-6)


def test_the_nodes_hold_the_generated_epoch_until_the_system_is_saved_and_the_stars_are_placed_by_them():
    system, _cfg = _multiple(3)
    assert all(node.epoch_unix == stored_trees.GENERATED_EPOCH_UNIX for node in system.orbit_nodes)
    offsets = system_tree.star_offsets_au([star.mass / constants.SOLAR_MASS_TO_KG for star in system.stars],
                                          system.star_slots, system.orbit_nodes, system.orbit_nodes[0].epoch_unix)
    assert len(offsets) == len(system.stars) and all(all(abs(v) < 1e9 for v in o) for o in offsets)


def test_a_system_of_more_than_two_stars_is_saved_and_read_back_with_its_tree(mysql_config):
    system, cfg = _multiple(3)
    conn = store.get_connection(mysql_config)
    try:
        system_id = store.insert_star_system(conn, system, cfg)
        conn.commit()
        loaded = store.load_star_system(conn, system_id)
        nodes = orbit_trees.read_nodes(conn, system_id)
        roles = [r["role"] for r in conn.execute(
            "SELECT role FROM stars WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchall()]
    finally:
        conn.close()
    assert roles == ["primary", "secondary"] + ["member"] * len(system.member_stars)
    assert len(loaded.member_stars) == len(system.member_stars) and loaded.star_slots == system.star_slots
    assert all(node.epoch_unix > 0 for node in nodes) and len({node.epoch_unix for node in nodes}) == 1
    assert [dataclasses.replace(n, epoch_unix=0.0) for n in nodes] == [
        dataclasses.replace(n, epoch_unix=0.0) for n in system.orbit_nodes]


def test_members_a_build_cannot_hold_stay_in_the_plan_as_dark_members_and_the_node_masses_are_the_drawn_ones():
    for seed in range(1500):
        tree = multistar.draw_member_tree(multistar.Placed("star", 0.2, 0.2), 4.0e9, draw.Stream(f"{seed}:members"))
        if tree.member_count >= 3 and any(m.kind in stored_trees.UNBUILT_KINDS for m in tree.members):
            break
    else:
        pytest.skip("no tree with a brown dwarf or remnant among three or more members")
    plan = stored_trees.drawn(tree)
    assert len(plan.members) == tree.member_count and len(plan.nodes) == len(tree.nodes)
    assert {m.member_no for m in plan.dark_members} == {m.member_no for m in tree.members if m.kind in stored_trees.UNBUILT_KINDS}
    root = next(n for n in plan.nodes if n.parent_node_no is None)
    assert root.mass_solar == pytest.approx(sum(m.mass_solar for m in plan.members), rel=1.0e-9)
    for node in plan.nodes:
        held = [m.slot for m in plan.members if m.parent_node_no == node.node_no] + [
            n.slot for n in plan.nodes if n.parent_node_no == node.node_no]
        assert sorted(held) == [0, 1]


def test_binary_system_off_is_one_star_and_on_is_a_pair_or_more():
    for seed in range(20):
        assert len(_system(seed, BINARY_SYSTEM=False)[0].stars) == 1
    pairs = [_system(seed, BINARY_SYSTEM=True)[0] for seed in range(20)]
    assert all(system.binary_type is not None or any(wanted.startswith("a companion star") for wanted in system.unmet_requirements)
               for system in pairs)


@pytest.mark.parametrize("wide", [False, True])
def test_wide_binary_forces_the_kind_of_the_pair(wide):
    kinds = {_system(seed, BINARY_SYSTEM=True, WIDE_BINARY=wide)[0].binary_type for seed in range(15)}
    assert kinds <= {"wide" if wide else "close", None}


def test_the_same_seed_draws_the_same_stars():
    first, _ = _system(7)
    second, _ = _system(7)
    assert first.to_dict() == second.to_dict()


@pytest.mark.parametrize("count", [1, 2, 3, 5])
def test_a_forced_star_count_is_the_number_of_stars(count):
    systems = [_system(seed, STARS=count)[0] for seed in range(12)]
    built = [system for system in systems if "stars" not in " ".join(map(str, system.unmet_requirements))
             and not any(wanted.startswith("a companion star") for wanted in system.unmet_requirements)]
    assert built and all(len(system.stars) == count for system in built)
    assert all(len(system.member_stars) == max(count - 2, 0) for system in built)


def test_the_star_cap_and_the_multi_star_prevalence_reach_the_draw():
    capped = [len(_system(seed, MAX_STARS=2)[0].stars) for seed in range(300)]
    assert max(capped) <= 2
    few = [len(_system(seed, PREVALENCE={"multi_star": -100})[0].stars) for seed in range(600)]
    many = [len(_system(seed, PREVALENCE={"multi_star": 100})[0].stars) for seed in range(600)]
    # the scale thins and thickens the table's counts, it does not forbid them (a white dwarf of a massive progenitor can still be one of five)
    assert sum(few) < sum(many) and sum(count >= 4 for count in few) < sum(count >= 4 for count in many)


def test_a_forced_pair_around_a_primary_under_the_hydrogen_burning_limit_is_refused_with_the_reason():
    """GEN.706: a late M primary under 0.08 solar masses has no lighter star to be a companion, and a brown dwarf is no star of
    the count (GEN.899); every other forced pair around it is built."""
    refused = 0
    for seed in range(120):
        system, _ = _system(seed, STAR_TYPE="M9V", BINARY_SYSTEM=True)
        unmet = [wanted for wanted in system.unmet_requirements if wanted.startswith("a companion star")]
        if unmet:
            refused += 1
            assert place_terms.born_mass_solar(system.primary_star) < multistar.companion_floor_solar(
                system.primary_star.metallicity_fe_h) and "hydrogen-burning limit" in unmet[0]
        else:
            assert len(system.stars) == 2
    assert refused > 0
