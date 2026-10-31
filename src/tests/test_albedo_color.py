# tests/test_albedo_color.py

"""GEN.912: the Bond albedo range of a bright world follows the colour of its star (physics/albedo_color.py)."""

import pytest

from planetgen import tuning
from planetgen.physics import albedo_color, climate_draws


def test_the_sun_is_the_reference():
    for regime in ("frozen", "wet"):
        assert albedo_color.spectral_factor(regime, 5770.0) == pytest.approx(1.0)
    assert albedo_color.shifted_range(tuning.ALBEDO_FROZEN_RANGE, "frozen", 5770.0) == pytest.approx(tuning.ALBEDO_FROZEN_RANGE)


def test_the_factors_are_kept_at_the_four_stars():
    for regime, factors in albedo_color.SPECTRAL_FACTOR.items():
        share = albedo_color.WET_SHARE if regime == "wet" else 1.0
        for temperature, value in zip(albedo_color.SPECTRAL_TEMPERATURES_K, factors):
            assert albedo_color.spectral_factor(regime, temperature) == pytest.approx(1.0 + share * (value - 1.0))


def test_the_frozen_factors_agree_with_the_freezing_lines_of_shields_2013():
    """19 percent less flux to freeze an M star's planet and 6 percent more for an F star's, against the G star's, at a G albedo of 0.55."""
    for factor, flux_ratio in ((albedo_color.SPECTRAL_FACTOR["frozen"][0], 0.74 / 0.93), (albedo_color.SPECTRAL_FACTOR["frozen"][3], 0.99 / 0.93)):
        absorbed_ratio = (1.0 - 0.55 * factor) * flux_ratio / (1.0 - 0.55)
        assert absorbed_ratio == pytest.approx(1.0, abs=0.06)


def test_a_red_star_darkens_and_a_blue_one_brightens_in_order():
    for regime in ("frozen", "wet"):
        factors = [albedo_color.spectral_factor(regime, t) for t in (2500.0, 3410.0, 4000.0, 5040.0, 5500.0, 5770.0, 6900.0, 9000.0)]
        assert factors == sorted(factors)
        assert factors[0] == factors[1] < 1.0 < factors[-1] == factors[-2]
    assert albedo_color.spectral_factor("frozen", 3410.0) < albedo_color.spectral_factor("wet", 3410.0)      # snow suffers most


def test_dry_rock_and_no_temperature_are_untouched():
    assert albedo_color.spectral_factor("dry", 3400.0) == 1.0
    assert albedo_color.spectral_factor("frozen", None) == 1.0
    assert albedo_color.spectral_factor("wet", 0.0) == 1.0
    assert climate_draws.albedo_range("t", 1e-6, 1361.0, 3400.0) == climate_draws.DEFAULT_ALBEDO_RANGE
    assert climate_draws.albedo_range("g", 0.3, 1361.0, 3400.0) == climate_draws.GIANT_BOND_ALBEDO_RANGE      # GEN.1030: a giant has its own range, which the star's colour does not move


def test_the_generator_range_follows_the_star():
    sun = 1361.0
    assert climate_draws.albedo_range("t", 1e-3, sun) == tuning.ALBEDO_WET_RANGE
    assert climate_draws.albedo_range("t", 1e-3, sun, 5770.0) == pytest.approx(tuning.ALBEDO_WET_RANGE)
    m_low, m_high = climate_draws.albedo_range("t", 1e-3, sun, 3400.0)
    assert m_low < tuning.ALBEDO_WET_RANGE[0] and m_high < tuning.ALBEDO_WET_RANGE[1]
    f_low, f_high = climate_draws.albedo_range("t", 1e-3, 0.3 * sun, 3400.0)
    assert f_high < tuning.ALBEDO_FROZEN_RANGE[1] and f_low < tuning.ALBEDO_FROZEN_RANGE[0]


def test_a_cold_star_keeps_the_ranges_below_one():
    low, high = albedo_color.shifted_range((0.9, 0.99), "frozen", 9000.0)
    assert high <= 0.95 and low <= high
