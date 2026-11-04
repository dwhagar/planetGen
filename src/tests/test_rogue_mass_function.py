# tests/test_rogue_mass_function.py

"""GEN.1112 RP2: the masses of the bodies a system throws out (`planetgen.physics.rogue_mass_function`): the cascade's slope, the cut at the scatterer's mass and the draw."""

import random

import pytest

from planetgen.physics import constants, rogue_mass_function as mf

EARTH = constants.EARTH_MASS_TO_KG
JUPITER = constants.JUPITER_MASS_TO_KG
SATURN = 95.16 * EARTH
SUMI_SLOPE = 0.96
SUMI_RANGE = (0.96 - 0.27, 0.96 + 0.47)           # [R4254, R4255], [R4254]: alpha = 0.96 +0.47 -0.27


def test_a_strength_cascade_gives_dohnanyis_mass_exponent():
    assert mf.supply_slope(3.5) == pytest.approx(5.0 / 6.0)          # n(m) ~ m^-11/6, so dN/dlogM ~ m^-5/6
    assert mf.supply_slope(4.0) == pytest.approx(1.0)                # the velocity-evolved cascade of [R3421]
    with pytest.raises(mf.RogueMassFunctionError):
        mf.supply_slope(1.0)


def test_the_ejected_slope_sits_inside_the_microlensing_measurement_without_being_fitted_to_it():
    alpha = mf.log_slope(0.33 * EARTH, 30.0 * EARTH, JUPITER)
    assert alpha == pytest.approx(5.0 / 6.0)
    assert SUMI_RANGE[0] < alpha < SUMI_RANGE[1]
    assert abs(alpha - SUMI_SLOPE) < 0.15


def test_nothing_heavier_than_the_scatterer_is_thrown_out():
    assert mf.mass_function(0.9 * JUPITER, JUPITER) > 0.0
    assert mf.mass_function(1.1 * JUPITER, JUPITER) == 0.0
    assert mf.log_slope(1.0 * EARTH, 2.0 * JUPITER, JUPITER) == 0.0
    assert mf.heaviest_ejected_kg(SATURN) == pytest.approx(SATURN)


def test_a_lighter_scatterer_throws_out_a_smaller_range_with_the_same_slope():
    assert mf.number_share(1.0 * JUPITER, 13.0 * JUPITER, 0.1 * EARTH, SATURN) == 0.0
    assert mf.log_slope(1.0 * EARTH, 10.0 * EARTH, SATURN) == pytest.approx(mf.log_slope(1.0 * EARTH, 10.0 * EARTH, JUPITER))


def test_the_shares_of_the_ejected_bodies_add_up_and_light_bodies_dominate_the_count():
    edges = [0.1 * EARTH, 2.0 * EARTH, 20.0 * EARTH, SATURN, JUPITER]
    shares = [mf.number_share(a, b, 0.1 * EARTH, JUPITER) for a, b in zip(edges, edges[1:])]
    assert sum(shares) == pytest.approx(1.0)
    assert shares == sorted(shares, reverse=True)
    assert shares[0] > 0.9                                            # terrestrial bodies are most of what a giant throws out [R0258]
    assert mf.number_share(0.0, 1e40, 0.1 * EARTH, JUPITER) == pytest.approx(1.0)


def test_the_draw_inverts_the_cumulative_number():
    lo, hi = 0.1 * EARTH, JUPITER
    assert mf.sample_mass_kg(0.0, lo, JUPITER) == pytest.approx(lo)
    assert mf.sample_mass_kg(1.0, lo, JUPITER) == pytest.approx(hi)
    middle = mf.sample_mass_kg(0.5, lo, JUPITER)
    assert mf.number_share(lo, middle, lo, JUPITER) == pytest.approx(0.5)
    rng = random.Random(4)
    masses = [mf.sample_mass_kg(rng.random(), lo, JUPITER) for _ in range(20000)]
    assert all(lo <= m <= hi * (1 + 1e-9) for m in masses)
    assert sum(m < 2.0 * EARTH for m in masses) / len(masses) == pytest.approx(mf.number_share(lo, 2.0 * EARTH, lo, JUPITER), abs=0.01)


def test_the_safronov_cut_comes_from_the_kick_module():
    sun = constants.SOLAR_MASS_TO_KG
    assert mf.can_eject(JUPITER, 7.1492e7, sun, 5.2)
    assert not mf.can_eject(EARTH, 6.371e6, sun, 1.0)


def test_bad_inputs_are_refused():
    with pytest.raises(mf.RogueMassFunctionError):
        mf.mass_function(0.0, JUPITER)
    with pytest.raises(mf.RogueMassFunctionError):
        mf.number_share(1.0, 2.0, 2.0 * JUPITER, JUPITER)
    with pytest.raises(mf.RogueMassFunctionError):
        mf.sample_mass_kg(1.5, EARTH, JUPITER)
    with pytest.raises(mf.RogueMassFunctionError):
        mf.heaviest_ejected_kg(0.0)
