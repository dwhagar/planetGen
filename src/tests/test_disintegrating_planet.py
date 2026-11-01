# src/tests/test_disintegrating_planet.py

"""GEN.1041: the wind of a disintegrating rocky planet against the numbers of [R3487] for KIC 12557548 b (Kepler-1520 b)."""

import math

import pytest

from planetgen.physics import constants, disintegrating_planet as dp

# eq. 6 of the paper: T_star 4450 K, R_star 0.65 solar, a = 0.013 AU
KEPLER_1520_B_K = dp.dayside_temperature_k(4450.0, 0.65 * constants.SOLAR_RADIUS_M, 0.013 * constants.AU_TO_M)


def test_the_dayside_of_kepler_1520_b_is_2150_k():
    assert KEPLER_1520_B_K == pytest.approx(2150.0, rel=0.01)


def test_the_vapour_pressure_of_olivine_at_2145_k_is_a_few_pascal():
    pascal = dp.vapour_pressure_pa(2145.0)
    assert pascal == pytest.approx(3.3, rel=0.1)                    # about 33 dyne/cm^2 from eq. 7 with the olivine constants
    assert dp.vapour_pressure_pa(2300.0) > 4.0 * pascal             # the exponential in the temperature


def test_a_planet_of_two_hundredths_of_an_earth_loses_a_tenth_to_one_earth_mass_per_gyr_and_lasts_tens_to_hundreds_of_myr():
    rate = dp.mass_loss_earth_per_gy(0.02, KEPLER_1520_B_K)
    assert 0.1 < rate < 1.0                                         # the observed lower limit and the fiducial (3.4)
    assert 0.04 < dp.lifetime_gy(0.02, KEPLER_1520_B_K) < 0.4       # 40 to 400 Myr (3.4)


def test_the_planet_that_survives_ten_gyr_is_near_a_tenth_of_an_earth_in_the_isothermal_wind_and_lighter_in_the_full_one():
    isothermal_low, isothermal_high = (dp.lifetime_gy(m, KEPLER_1520_B_K, factor=1.0) for m in (0.09, 0.13))
    assert isothermal_low < 10.0 < isothermal_high                  # the dividing mass is about 0.11 (3.3)
    full_low, full_high = (dp.lifetime_gy(m, KEPLER_1520_B_K) for m in (0.05, 0.11))
    assert full_low < 10.0 < full_high
    assert dp.lifetime_gy(0.3, KEPLER_1520_B_K) > 1.0e6             # a Mercury-size planet is not disintegrating


def test_the_rate_falls_steeply_with_the_mass_and_rises_with_the_temperature():
    rates = [dp.mass_loss_earth_per_gy(m, KEPLER_1520_B_K) for m in (0.01, 0.02, 0.04, 0.08)]
    assert rates == sorted(rates, reverse=True) and rates[0] > 300.0 * rates[-1]
    assert dp.mass_loss_earth_per_gy(0.05, 2400.0) > 100.0 * dp.mass_loss_earth_per_gy(0.05, 2000.0)


def test_the_mass_lost_in_an_age_is_under_the_mass_of_the_body():
    for mass in (0.02, 0.05, 0.1, 0.3):
        for age in (0.1, 1.0, 10.0):
            left = dp.remaining_mass_earth(mass, KEPLER_1520_B_K, age)
            assert 0.0 <= left <= mass
    # a body that outlives the age loses less than it has, by the rate at its start (the rate only rises as it thins, so this is a lower bound on the loss)
    for mass in (0.11, 0.2, 0.3):
        assert dp.DUTY_CYCLE * dp.mass_loss_earth_per_gy(mass, KEPLER_1520_B_K) * 10.0 < mass
        assert dp.remaining_mass_earth(mass, KEPLER_1520_B_K, 10.0) > 0.9 * mass


def test_a_planet_is_gone_after_its_lifetime_and_whole_before_a_tenth_of_it():
    life = dp.lifetime_gy(0.03, KEPLER_1520_B_K)
    assert dp.remaining_mass_earth(0.03, KEPLER_1520_B_K, 2.0 * life) == 0.0
    assert dp.remaining_mass_earth(0.03, KEPLER_1520_B_K, 0.05 * life) > 0.9 * 0.03


def test_the_gate_wants_a_light_hot_planet_on_an_orbit_under_a_day():
    assert dp.can_disintegrate(0.02, 2150.0, 0.65)                  # Kepler-1520 b: 15.7 h
    assert not dp.can_disintegrate(1.0, 2150.0, 0.65)               # the Earth is too heavy
    assert not dp.can_disintegrate(0.02, 1500.0, 0.65)              # too cool to vaporise
    assert not dp.can_disintegrate(0.02, 2150.0, 3.0)               # the period of a day or more


def test_the_loss_is_finite_for_a_body_of_a_ten_thousandth_of_an_earth_and_a_very_hot_one():
    for mass, kelvin in ((1.0e-4, 2150.0), (0.02, 3500.0), (0.5, 2150.0)):
        rate = dp.mass_loss_rate_kg_s(mass, kelvin)
        assert math.isfinite(rate) and rate >= 0.0
