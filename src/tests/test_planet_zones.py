# tests/test_planet_zones.py

"""GEN.314: the planet zones of a hierarchy of stars and the light a planet gets from every member."""

import math

import pytest

from planetgen import tuning
from planetgen.generation import multistar
from planetgen.physics import planet_zones as pz
from planetgen.physics.orbits import (
    hill_radius_in_pair_au, holman_wiegert_circumbinary_a_crit_au, holman_wiegert_critical_semimajor_axis)
from planetgen.physics.constants import AU_TO_KM as AU_KM
from planetgen.physics.system_tree import OrbitNode, TreeMember
from planetgen.util import draw


def node(no, parent, slot, mass, axis_au, e=0.0, kind="close"):
    return OrbitNode(no, parent, slot, mass, axis_au * AU_KM, e, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, kind)


def member(no, parent, slot, mass):
    return TreeMember(no, parent, slot, mass)


def zone(zones, kind, host):
    return next(z for z in zones if z.kind == kind and z.host_no == host)


def test_a_lone_star_has_one_window_from_its_floor_to_the_disk_cap():
    zones = pz.planet_zones([], [member(0, None, 0, 1.0)], floors_au={0: 0.4})
    assert zones == [pz.Zone("S", 0, 0.4, tuning.PLANET_ZONE_DISK_CAP_AU)]
    assert pz.planet_zones([], [member(0, None, 0, 1.0)], disk_cap_au=30.0)[0].outer_au == 30.0


def test_a_wide_pair_has_todays_limit_round_each_star():
    members = [member(0, 1, 0, 1.0), member(1, 1, 1, 0.4)]
    zones = pz.planet_zones([node(1, None, None, 1.4, 120.0, 0.3, "wide")], members)
    for number, own, other in ((0, 1.0, 0.4), (1, 0.4, 1.0)):
        expected = min(holman_wiegert_critical_semimajor_axis(120.0, other / 1.4, 0.3),
                       hill_radius_in_pair_au(120.0, own / 1.4, 0.3))
        assert zone(zones, "S", number).outer_au == pytest.approx(expected, rel=1e-12)
        assert zone(zones, "S", number).inner_au == 0.0


def test_a_close_pair_has_todays_circumbinary_floor_and_no_ceiling_but_the_disk():
    members = [member(0, 1, 0, 1.0), member(1, 1, 1, 0.6)]
    zones = pz.planet_zones([node(1, None, None, 1.6, 0.2, 0.1)], members)
    p = zone(zones, "P", 1)
    assert p.inner_au == pytest.approx(holman_wiegert_circumbinary_a_crit_au(0.2, 0.6 / 1.6, 0.1), rel=1e-12)
    assert p.outer_au == tuning.PLANET_ZONE_DISK_CAP_AU
    floored = zone(pz.planet_zones([node(1, None, None, 1.6, 0.2, 0.1)], members, floors_au={1: 3.0}), "P", 1)
    assert floored.inner_au == 3.0  # a giant's engulfment floor beats the fit


def _triple(outer_axis_au):
    """A close pair (nodes 2, members 0 and 1) round which a third star (member 2) goes at `outer_axis_au`."""
    nodes = [node(1, None, None, 2.4, outer_axis_au, 0.2, "wide"), node(2, 1, 0, 1.6, 0.3, 0.1)]
    members = [member(0, 2, 0, 1.0), member(1, 2, 1, 0.6), member(2, 1, 1, 0.8)]
    return nodes, members


def test_a_pair_inside_a_larger_orbit_is_capped_by_the_pairs_own_s_type_limit():
    nodes, members = _triple(20.0)
    zones = pz.planet_zones(nodes, members)
    pair = zone(zones, "P", 2)
    limit = min(holman_wiegert_critical_semimajor_axis(20.0, 0.8 / 2.4, 0.2), hill_radius_in_pair_au(20.0, 1.6 / 2.4, 0.2))
    assert pair.outer_au == pytest.approx(limit, rel=1e-12)
    assert pair.exists
    # The outer star sees the pair as one point of 1.6 solar masses.
    outer = zone(zones, "S", 2)
    assert outer.outer_au == pytest.approx(min(
        holman_wiegert_critical_semimajor_axis(20.0, 1.6 / 2.4, 0.2), hill_radius_in_pair_au(20.0, 0.8 / 2.4, 0.2)))
    # The root pair's P-type zone is the pair of (the pair, the third star): no ceiling but the disk.
    assert zone(zones, "P", 1).outer_au == tuning.PLANET_ZONE_DISK_CAP_AU


def test_a_pair_with_its_sibling_too_near_has_no_circumbinary_window():
    nodes, members = _triple(1.0)
    pair = zone(pz.planet_zones(nodes, members), "P", 2)
    assert not pair.exists and pair.width_au == 0.0


def test_a_group_node_gives_its_members_no_window():
    nodes = [node(1, None, None, 2.0, 5.0, kind="group")]
    members = [member(0, 1, 0, 1.0), member(1, 1, 1, 1.0)]
    assert pz.planet_zones(nodes, members) == []


def test_a_node_missing_a_child_is_refused():
    with pytest.raises(ValueError):
        pz.planet_zones([node(1, None, None, 1.0, 5.0)], [member(0, 1, 0, 1.0)])


def test_a_planets_flux_sums_over_every_member():
    nodes, members = _triple(20.0)
    luminosity = {0: 1.0, 1: 0.2, 2: 0.5}
    s_zone = zone(pz.planet_zones(nodes, members), "S", 2)
    far = 3.0
    # The third star's light, plus the close pair's two stars at the root orbit's mean inverse square.
    mean_inverse_square = 1.0 / (20.0 ** 2 * math.sqrt(1.0 - 0.2 ** 2))
    assert pz.flux_at(s_zone, far, luminosity, nodes, members) == pytest.approx(
        0.5 / far ** 2 + (1.0 + 0.2) * mean_inverse_square, rel=1e-12)
    p_zone = zone(pz.planet_zones(nodes, members), "P", 2)
    assert pz.flux_at(p_zone, 2.0, luminosity, nodes, members) == pytest.approx(
        (1.0 + 0.2) / 2.0 ** 2 + 0.5 * mean_inverse_square, rel=1e-12)
    # A star of the pair meets its partner at the pair's own orbit and the third star at the root's.
    s_in_pair = zone(pz.planet_zones(nodes, members), "S", 0)
    expected = 1.0 / 0.05 ** 2 + 0.2 / (0.3 ** 2 * math.sqrt(1 - 0.1 ** 2)) + 0.5 * mean_inverse_square
    assert pz.flux_at(s_in_pair, 0.05, luminosity, nodes, members) == pytest.approx(expected, rel=1e-12)


def test_a_lone_stars_flux_is_the_inverse_square_law():
    only = [member(0, None, 0, 1.0)]
    assert pz.flux_at(pz.planet_zones([], only)[0], 2.0, {0: 3.0}, [], only) == pytest.approx(0.75)
    with pytest.raises(ValueError):
        pz.flux_at(pz.planet_zones([], only)[0], 0.0, {0: 3.0}, [], only)


def test_a_planet_meeting_the_host_at_a_group_node_has_no_flux_sum():
    nodes = [node(1, None, None, 2.0, 5.0, kind="group")]
    members = [member(0, 1, 0, 1.0), member(1, 1, 1, 1.0)]
    host = pz.Zone("S", 0, 0.0, 10.0)
    with pytest.raises(ValueError):
        pz.flux_at(host, 1.0, {0: 1.0, 1: 1.0}, nodes, members)


# The table of docs/design/systems/multistar-and-compact-systems.md 5.6, on a seeded draw of the generator's own trees.
TABLE_S_TYPE = {2: 0.77, 3: 0.40, 4: 0.28, 5: 0.20, 6: 0.15, 7: 0.13}
"""dict: Stars with an S-type window above 0.3 AU, by number of stars, as the design page prints them."""

DRAWS = 300
"""int: Trees drawn for each number of stars."""


@pytest.mark.parametrize("stars", sorted(TABLE_S_TYPE))
def test_the_share_of_stars_with_an_s_type_window_falls_as_the_design_table_says(stars):
    wide = total = 0
    for seed in range(DRAWS):
        tree = multistar.draw_member_tree(multistar.Placed("star", 1.0), 5e9, draw.Stream(seed), n_stars=stars)
        if tree.member_count != stars:
            continue
        for found in pz.planet_zones(tree.nodes, tree.members):
            if found.kind == pz.S_TYPE:
                total += 1
                wide += found.width_au > 0.3
    assert total >= 0.8 * stars * DRAWS
    assert wide / total == pytest.approx(TABLE_S_TYPE[stars], abs=0.07)


# GEN.543: "usable" is a width above 0.3 AU, as for the S-type row; the design table's own row (34% falling to 18%) was
# an earlier estimate that neither of two candidate definitions reproduced (docs/design/systems/multistar-and-compact-systems.md 5.6).
MEASURED_P_TYPE = {2: 0.45, 3: 0.48, 4: 0.43, 5: 0.40, 6: 0.35, 7: 0.28}
"""dict: Pair nodes with a P-type window above 0.3 AU, by number of stars, on the same seeded draw (5 and 7 were re-measured
after GEN.270 weighted the tight tail of the period law, and 2 after GEN.272 moved the star draws; 2 again after GEN.1009 part 2 (the tree and the
build read one member record) moved the seeded draw, 0.49 to 0.45 on 300 pairs, which is inside the 0.03 the 300 draws allow, and 7 after GEN.1011 drew each
pair's period and eccentricity once, cut to the room it has, 0.33 to 0.28 on 1800 pairs, a draw change of 3.7 standard errors. Neither is
the dark members of GEN.1013, which leave the tree draw as it was)."""


@pytest.mark.parametrize("stars", sorted(MEASURED_P_TYPE))
def test_the_share_of_pair_nodes_with_a_p_type_window_is_the_measured_one(stars):
    wide = total = 0
    for seed in range(DRAWS):
        tree = multistar.draw_member_tree(multistar.Placed("star", 1.0), 5e9, draw.Stream(seed), n_stars=stars)
        if tree.member_count != stars:
            continue
        for found in pz.planet_zones(tree.nodes, tree.members):
            if found.kind == pz.P_TYPE:
                total += 1
                wide += found.width_au > 0.3
    assert total >= 0.8 * (stars - 1) * DRAWS
    assert wide / total == pytest.approx(MEASURED_P_TYPE[stars], abs=0.03)
