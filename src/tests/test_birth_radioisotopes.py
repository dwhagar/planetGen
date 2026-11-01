# tests/test_birth_radioisotopes.py

"""GEN.476: aluminium-26 and iron-60 at a star's birth from the gamma-ray masses, the star formation of the region and the granularity of the
sources (`physics/birth_radioisotopes.py`)."""

import statistics

import pytest

from planetgen.physics import birth_radioisotopes as birth
from planetgen.util import draw

SOLAR_SYSTEM_AL26 = 5.23e-5      # [R2128], via [R2768]: a test point only
SOLAR_SYSTEM_FE60 = 1.01e-8      # [R4307], via [R2768]
MILKY_WAY_DEPLETION_YR = 0.7e10 / 1.65   # cold gas over the star formation rate [R0422, R2680]
AL_MASS_FRACTION = 6.0e-5        # a solar-scale gas; the generator passes the star's own


def test_the_yields_per_mass_formed_are_the_gamma_ray_mass_over_rate_and_mean_life():
    assert birth.yield_per_solar_mass_formed("Al26") == pytest.approx(2.0 / (1.65 * 1.035e6))
    al = birth.yield_per_solar_mass_formed("Al26")
    fe = birth.yield_per_solar_mass_formed("Fe60")
    assert 3e-7 < fe < al < 2e-6
    with pytest.raises(ValueError):
        birth.yield_per_solar_mass_formed("Ca41")


def test_the_mean_ratio_is_about_ten_times_under_the_solar_systems():
    mean = birth.mean_ratio("Al26", MILKY_WAY_DEPLETION_YR, AL_MASS_FRACTION)
    assert 3e-6 < mean < 8e-6
    assert 5 < SOLAR_SYSTEM_AL26 / mean < 20


def test_a_faster_depleting_region_has_a_higher_ratio_and_a_poorer_gas_a_higher_one_too():
    base = birth.mean_ratio("Al26", 4e9, 6e-5)
    assert birth.mean_ratio("Al26", 2e9, 6e-5) == pytest.approx(2 * base)
    assert birth.mean_ratio("Al26", 4e9, 3e-5) == pytest.approx(2 * base)
    with pytest.raises(ValueError):
        birth.mean_ratio("Al26", 0.0, 6e-5)


def test_the_draws_keep_the_mean_and_are_spiky_for_the_short_lived_isotope():
    mean = birth.mean_ratio("Al26", MILKY_WAY_DEPLETION_YR, AL_MASS_FRACTION)
    stream = draw.Stream("test-al26")
    values = [birth.sample_ratio("Al26", mean, stream) for _ in range(3000)]
    assert statistics.mean(values) == pytest.approx(mean, rel=0.2)
    assert statistics.median(values) < 0.2 * mean            # most stars formed away from a recent addition
    fe_mean = birth.mean_ratio("Fe60", MILKY_WAY_DEPLETION_YR, 1.08e-3)
    fe = [birth.sample_ratio("Fe60", fe_mean, stream) for _ in range(3000)]
    assert statistics.median(fe) / fe_mean > 3 * statistics.median(values) / mean   # 60Fe averages over more events


def test_the_solar_system_ratio_is_in_the_upper_tail_of_the_distribution():
    mean = birth.mean_ratio("Al26", MILKY_WAY_DEPLETION_YR, AL_MASS_FRACTION)
    above = birth.fraction_above("Al26", mean, SOLAR_SYSTEM_AL26)
    assert 0.002 < above < 0.05
    assert 0.08 < birth.fraction_above("Al26", mean, 0.1 * SOLAR_SYSTEM_AL26) < 0.4     # the share rich enough to matter (row o6: 0.04 to 0.66)


def test_the_solar_systems_iron_60_is_not_far_from_the_distribution():
    mean = birth.mean_ratio("Fe60", MILKY_WAY_DEPLETION_YR, 1.08e-3)
    stream = draw.Stream("test-fe60")
    low = sum(birth.sample_ratio("Fe60", mean, stream) < SOLAR_SYSTEM_FE60 for _ in range(3000)) / 3000
    assert 0.02 < low < 0.4


def test_a_draw_is_deterministic_on_its_stream_and_a_zero_mean_gives_zero():
    a = birth.sample_ratio("Al26", 5e-6, draw.Stream("same"))
    b = birth.sample_ratio("Al26", 5e-6, draw.Stream("same"))
    assert a == b
    assert birth.sample_ratio("Al26", 0.0, draw.Stream("same")) == 0.0
    with pytest.raises(ValueError):
        birth.sample_ratio("Al26", -1.0, draw.Stream("same"))
