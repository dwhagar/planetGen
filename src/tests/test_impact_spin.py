# tests/test_impact_spin.py

"""
GEN.978: the birth spin of a rocky planet from the mass ratio of its last impacts (`planetgen.physics.impact_spin`, Kokubo and Genda 2010 equation 4
with a stand-in distribution of the ratio by planet mass).
"""

import math
import random

import pytest

from planetgen.physics import constants, impact_spin, spin

EARTH_MASS = constants.EARTH_MASS_TO_KG


def _radius_km(mass_earth):
    return 6371.0 * mass_earth ** 0.27


def _period_hours(mass_earth, rng):
    """One draw of a planet's birth day: gamma from the stand-in, three Gaussian components, capped as `spin.accretion_spin` caps it."""
    mass, radius = mass_earth * EARTH_MASS, _radius_km(mass_earth)
    critical = spin.critical_spin_rad_s(mass, radius)
    gamma = impact_spin.sample_gamma(mass_earth, rng.random(), rng.random())
    sigma = impact_spin.component_sigma_rad_s(critical, gamma)
    rate = math.sqrt(sum(rng.gauss(0.0, sigma) ** 2 for _ in range(3)))
    return 2.0 * math.pi / min(rate, spin.BREAKUP_FRACTION * critical) / 3600.0


def test_g_is_symmetric_and_has_the_values_of_the_paper():
    assert impact_spin.g_function(0.5) == pytest.approx(0.315, abs=0.001)
    assert impact_spin.g_function(0.2) == pytest.approx(impact_spin.g_function(0.8))
    assert impact_spin.g_function(0.0) == 0.0 == impact_spin.g_function(1.0)
    assert impact_spin.g_function(0.01) == pytest.approx(0.01, rel=0.15)      # a small impactor adds spin in proportion to its mass
    with pytest.raises(ValueError):
        impact_spin.g_function(1.2)


def test_the_equation_reproduces_the_paper_s_own_rms_at_its_mass_ratio():
    """Kokubo and Genda's realistic model: rms 2.35 rad/h at a breakup of 3.3 rad/h, 0.71. Equation 4 gives it for last impacts of mass ratio 0.28, a
    protoplanet colliding with a body about a third its size."""
    assert impact_spin.rms_breakup_fraction(0.28) == pytest.approx(2.35 / 3.3, rel=0.05)
    assert impact_spin.component_breakup_fraction(0.28) == pytest.approx((1.32 + 1.26 + 1.49) / 3.0 / 3.3, rel=0.05)   # the paper's three components


def test_an_equal_mass_impact_spins_a_planet_near_breakup_and_a_small_one_barely():
    assert impact_spin.rms_breakup_fraction(0.5) > 0.8
    assert impact_spin.rms_breakup_fraction(0.01) < 0.03


def test_the_chance_of_a_giant_impact_rises_with_the_planet_s_mass():
    chances = [impact_spin.giant_impact_probability(m) for m in (0.01, 0.05, 0.1, 0.3, 1.0, 5.0)]
    assert chances == sorted(chances)
    assert impact_spin.giant_impact_probability(impact_spin.ISOLATION_MASS_EARTH) == pytest.approx(0.5)
    assert chances[0] < 0.05 and chances[-1] > 0.99
    with pytest.raises(ValueError):
        impact_spin.giant_impact_probability(0.0)


def test_sample_gamma_stays_in_the_range_of_its_branch():
    low_giant, high_giant = impact_spin.GIANT_GAMMA_RANGE
    low_small, high_small = impact_spin.PLANETESIMAL_GAMMA_RANGE
    assert impact_spin.sample_gamma(1.0, 0.0, 0.0) == pytest.approx(low_giant)
    assert impact_spin.sample_gamma(1.0, 0.0, 1.0) == pytest.approx(high_giant)
    assert impact_spin.sample_gamma(0.01, 0.999, 0.0) == pytest.approx(low_small)
    assert impact_spin.sample_gamma(0.01, 0.999, 1.0) == pytest.approx(high_small)


def test_mean_rms_falls_with_smaller_planets_and_lies_between_the_two_branches():
    large, mars, small = (impact_spin.mean_rms_breakup_fraction(m) for m in (3.0, 0.107, 0.02))
    assert large > mars > small
    assert impact_spin.rms_breakup_fraction(impact_spin.PLANETESIMAL_GAMMA_RANGE[1]) < mars < impact_spin.rms_breakup_fraction(impact_spin.GIANT_GAMMA_RANGE[1])


def test_mars_is_reachable_and_an_earth_mass_planet_seldom_turns_that_slowly():
    """Mars turns in 24.6 hours (0.07 of breakup): as likely as not for a planet of its mass, rare for an Earth."""
    rng = random.Random(11)
    mars = sum(_period_hours(0.107, rng) >= 24.0 for _ in range(3000)) / 3000
    earth = sum(_period_hours(1.0, rng) >= 24.0 for _ in range(3000)) / 3000
    assert 0.25 < mars < 0.75
    assert 0.01 < earth < 0.2


def test_the_median_free_day_is_in_the_range_of_vesta_ceres_and_the_earth_moon_system():
    """Vesta 5.3 h, Ceres 9.1 h, the Earth after the Moon-forming impact about 5 h: planets of 0.1 to 5 Earth masses have a median day of 5 to 10 hours
    (the birth spin of the code today has a median of 2.9 h)."""
    rng = random.Random(3)
    days = sorted(_period_hours(10 ** rng.uniform(-1.0, 0.7), rng) for _ in range(6000))
    assert 5.0 < days[len(days) // 2] < 10.0


def test_the_range_is_the_widest_of_the_two_sources_in_question():
    """Kokubo and Genda: 2 to 3 hours; Miguel and Brunini 2010: 10 to 10,000 hours for most planets. The draws cover both ends."""
    rng = random.Random(4)
    days = sorted(_period_hours(10 ** rng.uniform(-1.3, 0.7), rng) for _ in range(8000))
    assert days[int(0.02 * len(days))] < 2.6
    assert days[int(0.98 * len(days))] > 100.0


def test_the_component_sigma_scales_with_the_breakup_rate_of_the_body():
    mass, radius = EARTH_MASS, 6371.0
    critical = spin.critical_spin_rad_s(mass, radius)
    sigma = impact_spin.component_sigma_rad_s(critical, 0.28)
    assert sigma == pytest.approx(impact_spin.component_breakup_fraction(0.28) * critical)
    twin_critical = spin.critical_spin_rad_s(mass, radius * 2.0 ** (1.0 / 3.0) * 0.5 ** (1.0 / 3.0))
    assert impact_spin.component_sigma_rad_s(twin_critical, 0.28) == pytest.approx(sigma)


def test_the_branch_is_the_one_the_gamma_comes_from():
    assert impact_spin.had_giant_impact(1.0, 0.5) and not impact_spin.had_giant_impact(0.01, 0.5)
    for mass_earth in (0.02, 0.1, 1.0):
        for branch_u in (0.05, 0.5, 0.95):
            giant = impact_spin.had_giant_impact(mass_earth, branch_u)
            low, high = impact_spin.GIANT_GAMMA_RANGE if giant else impact_spin.PLANETESIMAL_GAMMA_RANGE
            assert low <= impact_spin.sample_gamma(mass_earth, branch_u, 0.5) <= high


# --- Wired into the spin (GEN.979) ---------------------------------------------------------------------------------------

def test_the_birth_spin_of_a_small_planet_is_prograde_and_of_a_large_one_isotropic():
    """The planetesimal branch spins with the shear (tilts under 90 degrees); after a giant impact the axis is isotropic (the paper)."""
    from planetgen.util import draw
    draw.set_run_seed(979)
    small = [spin.accretion_spin(0.01 * EARTH_MASS, _radius_km(0.01))[1] for _ in range(3000)]
    large = [spin.accretion_spin(3.0 * EARTH_MASS, _radius_km(3.0))[1] for _ in range(3000)]
    giant_chance = impact_spin.giant_impact_probability(0.01)          # about one in a hundred of the small ones had a giant impact
    assert sum(t > 90.0 for t in small) / len(small) < 2.0 * giant_chance
    assert 0.42 < sum(t > 90.0 for t in large) / len(large) < 0.58


def test_the_wired_spin_is_the_stand_in_distribution_and_the_cap():
    from planetgen.util import draw
    draw.set_run_seed(980)
    mass, radius = EARTH_MASS, _radius_km(1.0)
    critical = spin.critical_spin_rad_s(mass, radius)
    rates = [spin.accretion_spin(mass, radius)[0] for _ in range(4000)]
    assert max(rates) <= spin.BREAKUP_FRACTION * critical + 1e-12
    days = sorted(2.0 * math.pi / rate / 3600.0 for rate in rates)
    assert 3.0 < days[len(days) // 2] < 12.0
    assert sum(day >= 24.0 for day in days) / len(days) < 0.2
