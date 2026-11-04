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
def test_every_triple_the_generator_draws_in_a_quadruple_passes_and_the_quadruple_verdicts_are_asked_of_both_shapes():
    """The generator tests each node as a nested triple, so every triple level of its four-star trees passes (the baseline GEN.1016's wiring
    keeps); the 3+1 and 2+2 verdicts exist for every tree. How many 3+1 trees the network refuses is the number the wiring brings to zero."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw
    shapes = set()
    for seed in range(1, 25):
        cfg = SystemConfig()
        cfg.STARS = 4
        cfg.PLANETS = False
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        members = [TreeMember(number, parent, slot, star.mass / constants.SOLAR_MASS_TO_KG)
                   for number, (star, (parent, slot)) in enumerate(zip(system.stars, system.star_slots))]
        found = hs.verdicts(system.orbit_nodes, members)
        assert all(verdict.stable for verdict in found if verdict.shape == "triple"), (seed, found)
        shapes |= {verdict.shape for verdict in found}
    assert {"3+1", "2+2"} <= shapes
