# tests/test_helium_core_white_dwarf.py

"""GEN.996: the helium-core white dwarf of a giant that fills its lobe on the red giant branch (Rappaport et al. 1995's
core mass and radius relation; the helium flash at a core of about 0.45 solar masses, Tauris and van den Heuvel 2006)."""

import dataclasses
import math
from statistics import NormalDist

import pytest

from planetgen import tuning
from planetgen.generation import multistar
from planetgen.physics import constants
from planetgen.physics import helium_core_white_dwarf as hecore
from planetgen.physics import roche


def test_the_giant_radius_is_rappaports_fit_for_population_one():
    # R = 5500 Mc^4.5 / (1 + 4 Mc^4) + 0.5 solar radii (Z 0.02): the numbers by hand
    for core, radius in ((0.15, 1.57), (0.3, 24.07), (0.4, 81.3), (0.45, 130.5)):
        assert hecore.giant_radius_rsun(core) == pytest.approx(radius, rel=0.01)


def test_a_metal_poorer_giant_is_smaller_at_the_same_core_between_the_two_published_constants():
    at_core = [hecore.giant_radius_rsun(0.4, fe_h) for fe_h in (0.5, 0.0, -0.5, -1.0, -1.5, -2.0)]
    assert at_core == sorted(at_core, reverse=True) and at_core[0] == at_core[1]       # held above Z 0.02
    assert hecore.giant_radius_rsun(0.4, -1.3) == pytest.approx(3300.0 * 0.4 ** 4.5 / (1.0 + 4.0 * 0.4 ** 4) + 0.5, rel=0.01)
    assert at_core[-1] == at_core[-2]                                                   # held under Z 0.001


def test_the_core_at_a_radius_inverts_the_radius_and_is_none_outside_the_relation():
    for core in (0.19, 0.25, 0.33, 0.44):
        assert hecore.core_mass_at_radius_sol(hecore.giant_radius_rsun(core)) == pytest.approx(core, abs=1.0e-9)
    assert hecore.core_mass_at_radius_sol(2.0) is None          # a subgiant, before the relation holds
    assert hecore.core_mass_at_radius_sol(200.0) is None        # past the radius the helium flash is reached at


def test_a_giant_in_a_close_pair_leaves_the_core_it_had_when_it_filled_its_lobe():
    giant, companion = 1.2, 0.66
    for core in (0.2, 0.3, 0.4):
        axis = hecore.overflow_axis_au(core, giant, companion)
        assert hecore.overflow_core_mass_sol(giant, companion, axis) == pytest.approx(core, abs=1.0e-6)
        # at periastron: an eccentric orbit of the same periastron leaves the same core
        wider = axis / (1.0 - 0.3)
        assert hecore.overflow_core_mass_sol(giant, companion, wider, eccentricity=0.3) == pytest.approx(core, abs=1.0e-6)


def test_the_axis_of_a_helium_white_dwarf_runs_from_a_few_solar_radii_to_under_two_au():
    low = hecore.overflow_axis_au(tuning.HELIUM_CORE_MIN_MASS_SOL, 1.2, 0.66)
    high = hecore.overflow_axis_au(tuning.HELIUM_FLASH_CORE_MASS_SOL, 1.2, 0.66)
    assert 0.02 < low < 0.1 and 1.0 < high < 2.0


def test_a_wide_pair_a_tight_pair_or_a_heavy_progenitor_leaves_no_helium_core_white_dwarf():
    assert hecore.overflow_core_mass_sol(1.2, 0.66, 3.0) is None        # the lobe is not reached before the flash
    assert hecore.overflow_core_mass_sol(1.2, 0.66, 0.02) is None       # filled while a subgiant: contact, not this
    assert hecore.overflow_core_mass_sol(2.5, 0.66, hecore.overflow_axis_au(0.3, 2.5, 0.66)) is None   # not a degenerate core


def test_the_lobe_is_eggletons():
    axis = 0.8
    lobe = axis * roche.roche_lobe_fraction(1.2, 0.66) * 215.03       # solar radii in an AU
    assert hecore.core_mass_at_radius_sol(lobe) == pytest.approx(hecore.overflow_core_mass_sol(1.2, 0.66, axis), abs=1.0e-4)


def test_the_generators_own_period_law_gives_about_the_20_pc_share_of_helium_core_white_dwarfs():
    """No rate is set for the channel: the share of a solar-type star's companions that sit where the lobe is reached with a
    core of 0.18 to 0.45 solar masses, times its chance of a companion, is the share of its white dwarfs that are helium-core.
    The 20 pc sample has 7.4 percent under 0.45 solar masses (Giammichele et al. 2012), a tenth of the solar neighbourhood's
    white dwarfs are low-mass (Liebert et al. 2005)."""
    mu, sigma = multistar.PERIOD_LAWS_LOG10_DAYS["solar"]
    for progenitor in (1.0, 1.2, 1.5):
        companion = 0.55 * progenitor
        low = hecore.overflow_axis_au(tuning.HELIUM_CORE_MIN_MASS_SOL, progenitor, companion)
        high = hecore.overflow_axis_au(tuning.HELIUM_FLASH_CORE_MASS_SOL, progenitor, companion)
        total = progenitor + companion
        log_period = [math.log10(365.25 * math.sqrt(axis ** 3 / total)) for axis in (low, high)]
        inside = NormalDist(mu, sigma).cdf(log_period[1]) - NormalDist(mu, sigma).cdf(log_period[0])
        weights = multistar.star_count_weights(progenitor)
        chance = 1.0 - weights[0] / sum(weights)
        assert 0.10 < inside < 0.18, progenitor
        assert 0.04 < chance * inside < 0.10, progenitor


# ---- the tree and the stored member (GEN.996, wired) ----

AGE_YR = 9.0e9


def _pair_tree(index, extra_star=False):
    """The tree of a star and a white dwarf (and a second star when `extra_star`) at nine gigayears, under stream `index`."""
    from planetgen.util import draw
    rng = draw.Stream(f"gen996:{index}")
    leaves = [multistar._star_leaf(0.8 + 0.6 * rng.random(), AGE_YR, rng, fe_h=0.0)]
    dwarf = multistar._remnant("white-dwarf", AGE_YR, 0.0, rng, fe_h=0.0)
    if dwarf is None:
        return None
    leaves.append(dwarf)
    if extra_star:
        leaves.append(multistar._star_leaf(0.5 + 0.5 * rng.random(), AGE_YR, rng, fe_h=0.0))
    return multistar._grow_tree(leaves, list(range(len(leaves))), AGE_YR, rng, None, math.inf, False, 1.0, 0.0)


def _helium_trees(extra_star, wanted, limit):
    found = []
    for index in range(limit):
        tree = _pair_tree(index, extra_star)
        if tree is not None and any(m.kind == "white-dwarf" and m.mass_solar < tuning.HELIUM_FLASH_CORE_MASS_SOL
                                    for m in tree.members):
            found.append(tree)
            if len(found) == wanted:
                break
    return found


def test_a_white_dwarf_the_tree_gives_under_the_flash_core_is_a_helium_core_in_a_post_common_envelope_pair():
    trees = _helium_trees(False, 3, 80)
    assert trees, "no helium-core white dwarf in 80 pairs"
    for tree in trees:
        dwarf = _the_helium_core(tree)
        node = tree.nodes[0]
        assert node.interaction == "post-common-envelope"
        assert tuning.HELIUM_CORE_MIN_MASS_SOL <= dwarf.mass_solar < tuning.HELIUM_FLASH_CORE_MASS_SOL
        # the axis is the one at which its giant, with this core, filled its lobe at periastron
        companion = next(m for m in tree.members if m is not dwarf)
        partner = companion.initial_mass_solar if companion.initial_mass_solar <= dwarf.initial_mass_solar else companion.mass_solar
        assert hecore.overflow_core_mass_sol(dwarf.initial_mass_solar, partner, node.separation_au, node.eccentricity) \
            == pytest.approx(dwarf.mass_solar, rel=1.0e-6)
        # the pair's mass and period follow the lighter white dwarf, and the progenitor's largest radius is the lobe it filled
        assert node.mass_solar == pytest.approx(sum(m.mass_solar for m in tree.members))
        assert node.period_years == pytest.approx(math.sqrt(node.separation_au ** 3 / node.mass_solar))
        assert dwarf.radius_max_au == pytest.approx(hecore.giant_radius_rsun(dwarf.mass_solar) * constants.SOLAR_RADIUS_M / constants.AU_M, rel=1.0e-3)


def test_the_nodes_above_a_helium_core_white_dwarf_lose_its_mass_and_keep_their_axes():
    trees = _helium_trees(True, 2, 160)
    assert trees, "no helium-core white dwarf in 160 triples"
    for tree in trees:
        for node in tree.nodes:
            mass_below = sum(m.mass_solar for m in tree.members
                             if _hangs_below(tree, m.parent_node_no, node.node_no))
            assert node.mass_solar == pytest.approx(mass_below, rel=1.0e-9)
            assert node.period_years == pytest.approx(math.sqrt(node.separation_au ** 3 / node.mass_solar), rel=1.0e-9)


def _the_helium_core(tree):
    return next(m for m in tree.members if m.kind == "white-dwarf" and m.mass_solar < tuning.HELIUM_FLASH_CORE_MASS_SOL)


def _hangs_below(tree, parent_no, node_no):
    nodes = {n.node_no: n for n in tree.nodes}
    while parent_no is not None:
        if parent_no == node_no:
            return True
        parent_no = nodes[parent_no].parent_node_no
    return False


def test_the_rebuilt_white_dwarf_of_a_stored_helium_core_is_the_member_the_tree_drew():
    from planetgen.generation import member_record, stored_trees
    tree = _helium_trees(False, 1, 80)[0]
    dwarf = _the_helium_core(tree)
    planned = stored_trees.PlannedStar(dwarf.member_no, 0, dwarf.slot, dwarf.kind, dwarf.mass_solar, dwarf.initial_mass_solar,
                                       dwarf.seed, dwarf.radius_au, AGE_YR)
    assert stored_trees.helium_core_mass_sol(planned) == dwarf.mass_solar
    record = member_record.make(dwarf.seed, dwarf.initial_mass_solar, AGE_YR, 0.0,
                                white_dwarf_mass_sol=stored_trees.helium_core_mass_sol(planned))
    assert record.mass_solar == pytest.approx(dwarf.mass_solar) and record.radius_au == pytest.approx(dwarf.radius_au, rel=1.0e-9)
    assert record.radius_max_au == pytest.approx(dwarf.radius_max_au, rel=1.0e-9)
    ordinary = member_record.make(dwarf.seed, dwarf.initial_mass_solar, AGE_YR, 0.0)
    assert ordinary.mass_solar >= tuning.WD_MASS_RANGE_SOL[0]
    assert stored_trees.helium_core_mass_sol(dataclasses.replace(planned, mass_solar=ordinary.mass_solar)) is None


def test_a_star_that_is_not_a_white_dwarf_ignores_the_core_mass():
    from planetgen.generation import member_record
    alive = member_record.make(7, 1.0, 1.0e9, 0.0, white_dwarf_mass_sol=0.3)
    plain = member_record.make(7, 1.0, 1.0e9, 0.0)
    assert alive == plain
