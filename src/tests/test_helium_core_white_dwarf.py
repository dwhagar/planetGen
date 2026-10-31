# tests/test_helium_core_white_dwarf.py

"""GEN.996: the helium-core white dwarf of a giant that fills its lobe on the red giant branch (Rappaport et al. 1995's
core mass and radius relation; the helium flash at a core of about 0.45 solar masses, Tauris and van den Heuvel 2006)."""

import math
from statistics import NormalDist

import pytest

from planetgen import tuning
from planetgen.generation import multistar
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
