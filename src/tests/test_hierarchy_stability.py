# tests/test_hierarchy_stability.py

"""GEN.1016: the stability of a tree of stars level by level (`physics/hierarchy_stability.py`): a triple at each level, a 3+1 quadruple judged by
the network of GEN.1029 with a margin on the systems it calls stable, and the fallbacks outside the network's range."""

import math

import pytest

from planetgen.physics import constants, hierarchy_stability as hs, quadruple_stability as quad, stability
from planetgen.physics.system_tree import OrbitNode, TreeMember

AU_KM = constants.AU_TO_KM


def node(node_no, parent, slot, mass, axis_au, e=0.0, inclination=0.0, kind="wide"):
    return OrbitNode(node_no=node_no, parent_node_no=parent, slot=slot, mass_solar=mass, separation_km=axis_au * AU_KM, eccentricity=e,
                     inclination_deg=inclination, ascending_node_deg=0.0, arg_periapsis_deg=0.0, mean_anomaly_deg=0.0, epoch_unix=0.0,
                     period_years=1.0, kind=kind, interaction="detached", flags=0)


def member(member_no, parent, slot, mass):
    return TreeMember(member_no=member_no, parent_node_no=parent, slot=slot, mass_solar=mass)


def three_plus_one(a_in, a_mid, a_out, masses=(1.0, 0.8, 0.6, 0.5), e=(0.1, 0.1, 0.1)):
    """Node 0 the inner pair, node 1 the middle star about it, node 2 the outer star about those."""
    m1, m2, m3, m4 = masses
    nodes = [node(0, 1, 0, m1 + m2, a_in, e[0]), node(1, 2, 0, m1 + m2 + m3, a_mid, e[1]), node(2, None, None, sum(masses), a_out, e[2])]
    members = [member(0, 0, 0, m1), member(1, 0, 1, m2), member(2, 1, 1, m3), member(3, 2, 1, m4)]
    return nodes, members


def two_plus_two(a_left, a_right, a_out, masses=(1.0, 0.8, 0.9, 0.7)):
    m1, m2, m3, m4 = masses
    nodes = [node(0, 2, 0, m1 + m2, a_left), node(1, 2, 1, m3 + m4, a_right), node(2, None, None, sum(masses), a_out)]
    members = [member(0, 0, 0, m1), member(1, 0, 1, m2), member(2, 1, 0, m3), member(3, 1, 1, m4)]
    return nodes, members


def test_one_or_two_stars_have_no_level_to_judge():
    assert hs.verdicts([], [member(0, 0, 0, 1.0)]) == []
    assert hs.verdicts([node(0, None, None, 1.8, 3.0)], [member(0, 0, 0, 1.0), member(1, 0, 1, 0.8)]) == []
    assert hs.tree_is_stable([], [])


def test_a_triple_level_clears_its_criterion_or_does_not():
    nodes = [node(0, 1, 0, 1.8, 1.0, 0.1), node(1, None, None, 2.5, 100.0, 0.2)]
    members = [member(0, 0, 0, 1.0), member(1, 0, 1, 0.8), member(2, 1, 1, 0.7)]
    [verdict] = hs.verdicts(nodes, members)
    assert verdict.shape == "triple" and verdict.stable and verdict.criterion == "eq4" and verdict.confident
    nodes[1] = node(1, None, None, 2.5, 2.0, 0.2)
    [verdict] = hs.verdicts(nodes, members)
    assert not verdict.stable


def test_a_well_separated_3p1_is_stable_by_the_network_and_a_crowded_one_is_not():
    [*_, wide] = hs.verdicts(*three_plus_one(0.05, 2.0, 80.0))
    assert wide.shape == "3+1" and wide.stable and wide.criterion == "network" and wide.confident
    [*_, tight] = hs.verdicts(*three_plus_one(1.0, 3.0, 8.0))
    assert tight.shape == "3+1" and not tight.stable


def test_the_axis_margin_rejects_a_system_the_plain_network_accepts():
    """Along the axis ratio the network turns unstable somewhere; just inside that edge the plain network says stable and the margin does not."""
    q_in, q_mid, q_out = 0.8, 0.5, 0.25
    edge = next(a for a in (x / 1000 for x in range(20, 460)) if not quad.is_stable(q_in, q_mid, q_out, a, a, 0.1, 0.1, 0.1))
    alpha = edge - 0.004
    assert quad.is_stable(q_in, q_mid, q_out, alpha, alpha, 0.1, 0.1, 0.1)
    masses = (1.0, 0.8, 0.9, 0.675)         # q_mid = 0.9/1.8 = 0.5, q_out = 0.675/2.7 = 0.25
    axes = (alpha * alpha * 400.0, alpha * 400.0, 400.0)
    plain = hs.quadruple_verdict(masses, axes, (0.1, 0.1, 0.1), (0.0, 0.0, 0.0), axis_margin=1.0)
    margin = hs.quadruple_verdict(masses, axes, (0.1, 0.1, 0.1), (0.0, 0.0, 0.0))
    assert plain[0] and plain[1] == "network"
    assert not margin[0] and margin[1] == "network"


def test_a_3p1_outside_the_networks_range_falls_back_to_the_nested_margin_and_says_so():
    # the outer star outweighs the others many times: q_out 8 is past the authors' 2
    masses = (1.0, 0.8, 0.6, 19.2)
    stable, criterion, confident, detail = hs.quadruple_verdict(masses, (0.05, 2.0, 2000.0), (0.1, 0.1, 0.1), (0.0, 0.0, 0.0))
    assert criterion == "nested" and not confident and "q_out" in detail
    assert stable
    tight = hs.quadruple_verdict(masses, (0.05, 2.0, 20.0), (0.1, 0.1, 0.1), (0.0, 0.0, 0.0))
    assert not tight[0] and not tight[2], "20 AU passes the plain nested test but not at the fallback margin of 1.5"


def test_a_node_failing_its_own_triple_makes_the_quadruple_unstable_without_asking_the_network():
    stable, criterion, confident, detail = hs.quadruple_verdict((1.0, 0.8, 0.6, 0.5), (1.0, 1.5, 50.0), (0.1, 0.1, 0.1), (0.0, 0.0, 0.0))
    assert not stable and criterion == "nested" and confident


def test_axes_that_do_not_increase_outward_are_not_a_hierarchy():
    assert not hs.quadruple_verdict((1.0, 0.8, 0.6, 0.5), (2.0, 1.0, 50.0), (0.1, 0.1, 0.1), (0.0, 0.0, 0.0))[0]


def test_a_stored_3p1_reports_its_two_triples_and_the_quadruple():
    found = hs.verdicts(*three_plus_one(0.05, 2.0, 80.0))
    assert sorted((v.node_no, v.shape) for v in found) == [(1, "triple"), (2, "3+1"), (2, "triple")]
    assert hs.tree_is_stable(*three_plus_one(0.05, 2.0, 80.0))
    assert not hs.tree_is_stable(*three_plus_one(1.0, 3.0, 8.0))


def test_a_wide_2p2_is_stable_by_the_network_and_a_crowded_one_is_not():
    [wide] = [v for v in hs.verdicts(*two_plus_two(0.05, 0.07, 100.0)) if v.shape == "2+2"]
    assert wide.stable and wide.criterion == "network" and wide.confident
    found = hs.verdicts(*two_plus_two(0.5, 0.7, 2.0))
    assert [v.shape for v in found].count("triple") == 2
    [crowded] = [v for v in found if v.shape == "2+2"]
    assert not crowded.stable and crowded.criterion == "nested"


def test_the_axis_margin_rejects_a_2p2_the_plain_network_accepts():
    q_in1, q_in2, q_out = 0.8, 0.78, 0.8
    edge = next(a for a in (x / 1000 for x in range(20, 460)) if not quad.is_stable_2p2(q_in1, q_in2, q_out, a, a, 0.1, 0.1, 0.1))
    alpha = edge - 0.004
    assert quad.is_stable_2p2(q_in1, q_in2, q_out, alpha, alpha, 0.1, 0.1, 0.1)
    masses = (1.0, 0.8, 0.9, 0.702)         # q_in2 = 0.78, q_out = 1.602 / 1.8 is below 1
    axes = (alpha * 100.0, alpha * 100.0, 100.0)
    plain = hs.quadruple_verdict_2p2(masses, axes, (0.1, 0.1, 0.1), (0.0, 0.0, 0.0), axis_margin=1.0)
    margin = hs.quadruple_verdict_2p2(masses, axes, (0.1, 0.1, 0.1), (0.0, 0.0, 0.0))
    assert plain[0] and plain[1] == "network"
    assert not margin[0] and margin[1] == "network"


def test_a_2p2_outside_the_networks_range_falls_back_and_says_so():
    """A faint companion in the second binary (q_in2 0.05, below the 0.1 the network was trained on) is judged by the nested triples at 1.5."""
    stable, criterion, confident, detail = hs.quadruple_verdict_2p2((1.0, 0.8, 0.9, 0.045), (0.05, 0.05, 100.0), (0.1, 0.1, 0.1), (0.0, 0.0, 0.0))
    assert stable and criterion == "nested" and not confident and "q_in2" in detail
    wide = hs.quadruple_verdict_2p2((1.0, 0.8, 0.9, 0.7), (0.0005, 0.0005, 30.0), (0.1, 0.1, 0.1), (0.0, 0.0, 0.0))
    assert wide[1] == "network" and wide[2]


def test_a_young_group_node_is_not_a_hierarchy():
    nodes = [node(0, None, None, 3.0, 5.0, kind="group")]
    assert hs.verdicts(nodes, [member(0, 0, 0, 1.0), member(1, 0, 1, 1.0), member(2, 0, 1, 1.0)]) == []


def test_mutual_inclination_of_nodes_reaches_the_network():
    nodes, members = three_plus_one(0.05, 2.0, 80.0)
    tilted = [nodes[0], nodes[1], OrbitNode(**{**nodes[2].__dict__, "inclination_deg": 120.0})]
    assert math.isclose(hs._mutual(tilted[1], tilted[2]), math.radians(120.0), rel_tol=1e-9)


# --- what today's generator draws -----------------------------------------------------------------------------------------------------------

@pytest.mark.slow
def _four_star_system(seed):
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw
    cfg = SystemConfig()
    cfg.STARS = 4
    cfg.PLANETS = False
    with draw.bound(seed):
        system = StarSystem(system_config=cfg)
    members = [TreeMember(number, parent, slot, star.mass / constants.SOLAR_MASS_TO_KG)
               for number, (star, (parent, slot)) in enumerate(zip(system.stars, system.star_slots or ()))]
    return system, members


def test_every_quadruple_the_generator_draws_clears_the_network_at_both_shapes():
    """GEN.1016: the tree draw asks `hierarchy_stability.tree_is_stable` of a tree of four stars and draws it again when a level refuses (`_Builder.hierarchy_stable`),
    so every triple level and every 3+1 and 2+2 verdict of a drawn four-star tree passes. Before the wiring 44 percent of the 3+1 trees and 22 percent of the 2+2 trees
    drawn on main were refused by the network."""
    shapes = set()
    for seed in range(1, 25):
        system, members = _four_star_system(seed)
        found = hs.verdicts(system.orbit_nodes, members)
        assert all(verdict.stable for verdict in found), (seed, found)
        shapes |= {verdict.shape for verdict in found}
    assert {"3+1", "2+2"} <= shapes


def test_a_close_pair_is_circular_in_the_tree_that_is_judged():
    """GEN.1615: a pair of two single stars closer than `tuning.WIDE_BINARY_SEPARATION_MIN_AU` is circular at the tree draw (eccentricity 0, no periapsis), so the tree
    the network judges is the tree that is stored. Seed 112 drew an eccentric close 2+2 pair, passed, and was stored circular by `StarSystem`, which the network then refused."""
    from planetgen import tuning
    closes = 0
    for seed in [112, *range(1, 25)]:
        system, members = _four_star_system(seed)
        for node in system.orbit_nodes:
            under = [member for member in members if member.parent_node_no == node.node_no]
            if len(under) == 2 and node.separation_km / constants.AU_TO_KM < tuning.WIDE_BINARY_SEPARATION_MIN_AU:
                closes += 1
                assert node.eccentricity == 0.0 and node.arg_periapsis_deg == 0.0, (seed, node)
        found = hs.verdicts(system.orbit_nodes, members)
        assert all(verdict.stable for verdict in found), (seed, found)
    assert closes


def test_a_tree_the_network_always_refuses_loses_its_lightest_loose_member(monkeypatch):
    """The restart is bounded: when no draw of a four-star shape passes, the tree drops a member as it does for a shape no overlap screen passes (`tree.member_dropped`)
    and the three stars left have no quadruple level to refuse."""
    monkeypatch.setattr(hs, "tree_is_stable", lambda nodes, members, **margins: False)
    system, _ = _four_star_system(3)
    assert len(system.stars) <= 3


def test_a_tree_of_three_stars_is_not_asked_of_the_quadruple_network(monkeypatch):
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw
    calls = []
    monkeypatch.setattr(hs, "tree_is_stable", lambda nodes, members, **margins: calls.append(1) or False)
    cfg = SystemConfig()
    cfg.STARS = 3
    cfg.PLANETS = False
    with draw.bound(5):
        system = StarSystem(system_config=cfg)
    assert len(system.stars) == 3 and not calls
