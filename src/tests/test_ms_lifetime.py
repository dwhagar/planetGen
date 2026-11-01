# tests/test_ms_lifetime.py

"""
GEN.467: the main-sequence lifetime of `physics/ms_lifetime.py`, the fit of [R2055] in place of the
`10 Gy * M^-2.5` law: the published scale of its Table 1, the metallicity dependence, the inverse the companion
draws use, and the places that read it.
"""

import math

import pytest

from planetgen import tuning
from planetgen.generation import multistar
from planetgen.physics import ms_lifetime, stellar_evolution

MASSES = (0.1, 0.2, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 10.0, 20.0, 50.0, 100.0, 150.0)


@pytest.mark.parametrize("mass, order_myr", [(1.0, 1.0e4), (2.0, 1.0e3), (5.0, 1.0e2)])
def test_the_time_to_the_giant_branch_has_the_scale_of_the_papers_table_1(mass, order_myr):
    """Table 1 lists tBGB as about 1e4, 1e3 and 1e2 Myr at 1, 2 and 5 solar masses (Z = 0.02)."""
    reference = math.log10(0.02 / tuning.MS_LIFETIME_REFERENCE_Z)
    time_myr = ms_lifetime.base_of_giant_branch_gy(mass, fe_h=reference) * 1.0e3
    assert order_myr / 3.0 < time_myr < order_myr * 3.0


def test_a_solar_mass_of_the_reference_composition_lasts_about_ten_gigayears():
    assert ms_lifetime.lifetime_gy(1.0) == pytest.approx(10.2, abs=0.1)


def test_the_lifetime_falls_with_mass_over_the_whole_range():
    times = [ms_lifetime.lifetime_gy(m) for m in MASSES]
    assert times == sorted(times, reverse=True)
    fine = [ms_lifetime.lifetime_gy(0.05 * 6000.0 ** (i / 2000)) for i in range(2001)]
    assert fine == sorted(fine, reverse=True)


def test_the_massive_stars_flatten_to_a_few_million_years():
    assert 2.5e-3 < ms_lifetime.lifetime_gy(150.0) < 4.0e-3
    assert ms_lifetime.lifetime_gy(20.0) == pytest.approx(9.0e-3, rel=0.1)


def test_a_poorer_star_lives_shorter_and_a_richer_one_longer():
    assert (ms_lifetime.lifetime_gy(1.0, fe_h=-2.0) < ms_lifetime.lifetime_gy(1.0, fe_h=-1.0)
            < ms_lifetime.lifetime_gy(1.0) < ms_lifetime.lifetime_gy(1.0, fe_h=0.3))


@pytest.mark.parametrize("fe_h", [-1.5, 0.0, 0.2])
@pytest.mark.parametrize("mass", MASSES)
def test_the_inverse_finds_the_mass_back(mass, fe_h):
    lifetime = ms_lifetime.lifetime_gy(mass, fe_h)
    assert ms_lifetime.mass_for_lifetime_gy(lifetime, fe_h) == pytest.approx(mass, rel=1e-6)


def test_the_inverse_clamps_outside_the_range():
    assert ms_lifetime.mass_for_lifetime_gy(1.0e9) == ms_lifetime.MIN_MASS_SOL
    assert ms_lifetime.mass_for_lifetime_gy(1.0e-6) == ms_lifetime.MAX_MASS_SOL


def test_the_stellar_evolution_function_is_the_fit():
    for mass in MASSES:
        assert stellar_evolution.main_sequence_lifetime_gy(mass) == ms_lifetime.lifetime_gy(mass)


def test_the_oldest_sun_like_star_the_universe_holds_is_a_bit_under_a_solar_mass():
    mass = ms_lifetime.mass_for_lifetime_gy(tuning.UNIVERSE_AGE_GY)
    assert 0.9 < mass < 0.96


def test_the_companion_draw_caps_its_mass_by_the_age():
    """`_alive_star_mass` draws below the heaviest mass still on the main sequence at the system's age."""
    import random
    rng = random.Random(3)
    age_yr = 1.0e8
    cap = ms_lifetime.mass_for_lifetime_gy(age_yr / 1.0e9)
    assert 4.0 < cap < 6.0
    assert all(multistar._alive_star_mass(0.1, 50.0, age_yr, rng) <= cap * 1.0001 for _ in range(200))


@pytest.mark.parametrize("fe_h", [-2.0, -0.4, 0.0, 0.25])
@pytest.mark.parametrize("lifetime", [1.0e-4, 0.003, 0.05, 1.7, 9.0, 40.0, 1.0e3])
def test_the_mass_for_a_lifetime_is_what_the_bisection_over_lifetime_gy_gives_bit_for_bit(fe_h, lifetime):
    """PERF.163: the search reads the metallicity's terms once; the mass is what the plain bisection over
    `lifetime_gy` found."""
    low, high = math.log(ms_lifetime.MIN_MASS_SOL), math.log(ms_lifetime.MAX_MASS_SOL)
    if lifetime >= ms_lifetime.lifetime_gy(ms_lifetime.MIN_MASS_SOL, fe_h):
        wanted = ms_lifetime.MIN_MASS_SOL
    elif lifetime <= ms_lifetime.lifetime_gy(ms_lifetime.MAX_MASS_SOL, fe_h):
        wanted = ms_lifetime.MAX_MASS_SOL
    else:
        for _ in range(ms_lifetime._BISECTIONS):
            middle = 0.5 * (low + high)
            if ms_lifetime.lifetime_gy(math.exp(middle), fe_h) > lifetime:
                low = middle
            else:
                high = middle
        wanted = math.exp(0.5 * (low + high))
    assert ms_lifetime.mass_for_lifetime_gy(lifetime, fe_h) == wanted
