# tests/test_brown_dwarf_companions.py

"""
GEN.331: the substellar members of a drawn tree. A brown dwarf companion (13 to 80 Jupiter masses) has 3.2% of solar-type stars at 28 to 1,590 AU and
under 0.6% at 3 to 10 AU (the desert, [R3048]); 15% of them are a pair of two brown dwarfs; and a planetary-mass companion (5 to 13 Jupiter masses,
`brown_dwarf.is_planetary_mass`) is one more member on its own chance. The orbit of each about the rest of its tree is the law of `multistar.period_days`.
"""

import math
from collections import Counter

import pytest

from planetgen import tuning
from planetgen.generation import multistar
from planetgen.physics import brown_dwarf, stability
from planetgen.util import draw

TREES = 20000
AGE_YR = 4.5e9
SOLAR_MASS = 1.0


def _members_under(tree, node_no, slot):
    """Every member that hangs below the `slot` child of node `node_no`."""
    nodes = {node.node_no: node for node in tree.nodes}
    found = [member for member in tree.members if (member.parent_node_no, member.slot) == (node_no, slot)]
    for node in tree.nodes:
        if (node.parent_node_no, node.slot) == (node_no, slot):
            found += [member for member in tree.members if _hangs_below(nodes, member.parent_node_no, node.node_no)]
    return found


def _hangs_below(nodes, node_no, ancestor):
    while node_no is not None:
        if node_no == ancestor:
            return True
        node_no = nodes[node_no].parent_node_no
    return False


def _companions(tree):
    """`(axis_au, members)` of every node of `tree` that joins an all-substellar side to a side that is not (a brown dwarf, a pair of them or a planetary-mass
    companion about the rest of the tree), and `axis_au` of every node that joins two all-substellar sides (a pair of brown dwarfs)."""
    companions, pairs = [], []
    for node in tree.nodes:
        sides = [_members_under(tree, node.node_no, slot) for slot in (0, 1)]
        substellar = [bool(side) and all(member.kind == "brown-dwarf" for member in side) for side in sides]
        axis_au = node.separation_km / 1.495978707e8
        if all(substellar):
            pairs.append(axis_au)
        elif any(substellar):
            companions.append((axis_au, sides[substellar.index(True)]))
    return companions, pairs


@pytest.fixture(scope="module")
def census():
    window = desert = hosts = pairs_of_dwarfs = planetary = 0
    planetary_masses, lone_axes = [], []
    for seed in range(TREES):
        tree = multistar.draw_member_tree(multistar.Placed("star", SOLAR_MASS, SOLAR_MASS, fe_h=0.0), AGE_YR, draw.Stream(seed))
        companions, pair_axes = _companions(tree)
        for axis_au, side in companions:
            mass = max(member.mass_solar for member in side)
            if brown_dwarf.is_planetary_mass(mass):
                planetary += 1
                planetary_masses.append(mass)
                continue
            hosts += 1
            pairs_of_dwarfs += len(side) == 2
            window += 28.0 <= axis_au < 1590.0
            desert += 3.0 <= axis_au < 10.0
            if len(tree.members) == 2:
                lone_axes.append(axis_au)     # a star and one brown dwarf: the pair has the whole room (the root pair)
    return Counter(window=window, desert=desert, hosts=hosts, pairs=pairs_of_dwarfs, planetary=planetary), planetary_masses, lone_axes


def _within_four_sigma(count, expected_share, draws=TREES):
    expected = expected_share * draws
    return abs(count - expected) <= 4.0 * math.sqrt(draws * expected_share * (1.0 - expected_share))


def test_brown_dwarf_companions_at_28_to_1590_au_are_3_2_percent_of_stars(census):
    counts, _masses, _axes = census
    assert _within_four_sigma(counts["window"], 0.032), counts


def test_the_desert_holds_under_0_6_percent_of_stars_and_the_rest_of_the_run_is_flat_in_log_axis(census):
    counts, _masses, axes = census
    assert counts["desert"] / TREES < 0.006, counts
    inner, middle = sum(10.0 <= a < 100.0 for a in axes), sum(100.0 <= a < 1000.0 for a in axes)
    assert abs(inner - middle) <= 4.0 * math.sqrt(inner + middle), (inner, middle)


def test_fifteen_percent_of_the_brown_dwarf_companions_are_a_pair_of_two(census):
    counts, _masses, _axes = census
    share = tuning.BROWN_DWARF_PAIR_SHARE
    assert abs(counts["pairs"] - share * counts["hosts"]) <= 4.0 * math.sqrt(counts["hosts"] * share * (1.0 - share)), counts


def test_a_planetary_mass_companion_is_under_the_deuterium_limit_and_at_its_own_chance(census):
    counts, masses, _axes = census
    low, high = brown_dwarf.planetary_mass_range_solar()
    assert masses and all(low <= mass < high for mass in masses)
    assert all(brown_dwarf.is_planetary_mass(mass) for mass in masses)
    assert _within_four_sigma(counts["planetary"], tuning.PLANETARY_MASS_COMPANION_CHANCE), counts


def test_a_drawn_brown_dwarf_is_over_the_deuterium_limit_and_under_the_hydrogen_limit():
    low, high = brown_dwarf.mass_range_solar()
    masses = [brown_dwarf.sample_mass_solar(draw.Stream(seed)) for seed in range(200)]
    assert all(low <= mass <= high for mass in masses) and not any(brown_dwarf.is_planetary_mass(mass) for mass in masses)


def test_the_companion_law_is_flat_in_axis_outside_the_desert_and_thin_inside_it():
    mass = 1.05
    axes = [stability.axis_au(multistar.period_days("brown-companion", mass, 0.02, math.inf, False, 1.0, (k + 0.5) / 20000) / 365.25, mass)
            for k in range(20000)]
    edge, widest = tuning.BROWN_DWARF_AXIS_RANGE_AU
    assert max(axes) <= widest * 1.0001
    inside = sum(a < edge for a in axes) / len(axes)
    expected = tuning.BROWN_DWARF_DESERT_KEEP * math.log10(edge / 0.02)
    expected /= expected + math.log10(widest / edge)
    assert inside == pytest.approx(expected, abs=0.002)
    decade = [sum(lo <= a < lo * 10.0 for a in axes) for lo in (10.0, 100.0)]
    assert decade[0] == pytest.approx(decade[1], rel=0.02)


def test_the_pair_law_peaks_near_four_au_and_is_cut_to_the_room():
    mass = 0.06
    axes = sorted(stability.axis_au(multistar.period_days("brown-pair", mass, 0.01, math.inf, False, 1.0, (k + 0.5) / 4000) / 365.25, mass)
                  for k in range(4000))
    assert axes[2000] == pytest.approx(10.0 ** tuning.BROWN_DWARF_PAIR_AXIS_LOG10_AU[0], rel=0.05)
    cut = [stability.axis_au(multistar.period_days("brown-pair", mass, 0.5, 3.0, False, 1.0, (k + 0.5) / 400) / 365.25, mass) for k in range(400)]
    assert 0.5 <= min(cut) and max(cut) <= 3.0 * 1.0001
    assert multistar.period_days("brown-pair", mass, 5.0, 3.0, False, 1.0, 0.5) is None


def test_a_forced_count_adds_no_substellar_member():
    for seed in range(400):
        tree = multistar.draw_member_tree(multistar.Placed("star", SOLAR_MASS, SOLAR_MASS, fe_h=0.0), AGE_YR, draw.Stream(seed), n_stars=2)
        assert tree.member_count == 2 and all(member.kind != "brown-dwarf" for member in tree.members)


def test_the_cooling_track_continues_under_the_lightest_mass_the_table_holds():
    lightest = 0.01
    at_edge = brown_dwarf.state(lightest, 3.0)
    below = [brown_dwarf.state(mass, 3.0) for mass in (0.008, 0.006, brown_dwarf.planetary_mass_range_solar()[0])]
    assert all(state["temperature_k"] < at_edge["temperature_k"] for state in below)
    assert [state["temperature_k"] for state in below] == sorted((state["temperature_k"] for state in below), reverse=True)
    assert all(0.08 < state["radius_sol"] < 0.1 for state in below)
    assert all(state["spectral_type"].startswith("Y") for state in below)
