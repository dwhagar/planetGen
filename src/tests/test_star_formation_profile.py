# tests/test_star_formation_profile.py

"""GEN.134: the recent star-formation factor is derived from the shape, normalised, bounded, and plausible.

The Milky Way figures here (a rate peak between 3 and 6 kpc, 84% of the star formation inside the solar circle, 1%
beyond 13.5 kpc, a fall of about 0.28 dex per kpc, a young tail of 3% to 23% in the bulge) are test targets that
show the derived profile is plausible. None of them is an input of the generator."""

import math

import numpy as np
import pytest

from planetgen import tuning
from planetgen.galaxy import density, star_formation

SHAPE = density.build_galaxy_shape(
    disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0, bulge_amplitude=3.11,
    arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4,
)
RADII = np.linspace(20.0, 25_000.0, 6000)


def _rings(shape):
    """The star formation in each ring of a grid of radii: the rate per area times the ring's circumference."""
    return RADII * star_formation.star_formation_rate_profile(RADII, shape)


def test_the_factor_averages_one_over_the_disk():
    """The disk makes as many stars as before: stars weighted by the disk's own surface density average to 1."""
    weight = RADII * np.exp(-RADII / SHAPE.disk_scale_length_pc)
    mean = np.trapezoid(weight * star_formation.recent_factor(RADII, SHAPE), RADII) / np.trapezoid(weight, RADII)
    assert mean == pytest.approx(1.0, abs=2e-3)


def test_the_profile_rings_at_the_bars_end():
    """With the bar gap the rate per area peaks near the gap's edge, not at the centre (target: the Milky Way's
    molecular ring, 4 to 5 kpc)."""
    peak = RADII[np.argmax(star_formation.star_formation_rate_profile(RADII, SHAPE))]
    assert 3000.0 <= peak <= 6000.0


def test_most_star_formation_is_inside_the_solar_circle_and_little_is_far_out():
    rings = _rings(SHAPE)
    total = np.trapezoid(rings, RADII)
    solar = tuning.REFERENCE_RADIUS_SCALE_LENGTHS * SHAPE.disk_scale_length_pc
    inside = np.trapezoid(rings[RADII <= solar], RADII[RADII <= solar]) / total
    far = np.trapezoid(rings[RADII >= 13_500.0], RADII[RADII >= 13_500.0]) / total
    assert 0.75 <= inside <= 0.95         # observed 84%
    assert far <= 0.03                    # observed 1%


def test_the_outer_fall_is_steeper_than_the_stars_and_near_the_observed_slope():
    """Beyond the peak the rate falls as the gas law steepens the gas disk: n / L, 0.23 dex per kpc here, against
    the observed 0.28."""
    rate = star_formation.star_formation_rate_profile(RADII, SHAPE)
    peak = RADII[np.argmax(rate)]
    window = (RADII > peak + 2000.0) & (RADII < peak + 8000.0)
    slope = np.polyfit(RADII[window] / 1000.0, np.log10(rate[window]), 1)[0]
    assert -0.35 <= slope <= -0.15
    stars = -math.log10(math.e) / (SHAPE.disk_scale_length_pc / 1000.0)
    assert slope < stars


def test_the_recent_share_falls_outward_and_in_the_gap():
    solar = tuning.REFERENCE_RADIUS_SCALE_LENGTHS * SHAPE.disk_scale_length_pc
    gap = star_formation.gap_radius_pc(SHAPE)
    edge = star_formation.recent_factor(1.4 * gap, SHAPE)
    assert star_formation.recent_factor(300.0, SHAPE) < 0.2 * edge
    assert star_formation.recent_factor(solar, SHAPE) < edge
    assert star_formation.recent_factor(3 * solar, SHAPE) < star_formation.recent_factor(solar, SHAPE)


def _without_thin_gas_steepening(monkeypatch):
    """The plain law (GEN.134): a saturation density so low that every ring counts as dense gas."""
    monkeypatch.setattr(tuning, "STAR_FORMATION_SATURATION_GAS_MSUN_PC2", 1e-12)
    star_formation._normalisation.cache_clear()


def test_the_factor_scales_with_the_galaxy(monkeypatch):
    """A galaxy twice as large in every length has the same factor at twice the radius (no Milky Way length), once the
    gas density that sets where the index rises (GEN.450, a surface density that a bigger galaxy does not share) is
    left out."""
    _without_thin_gas_steepening(monkeypatch)
    big = density.build_galaxy_shape(
        disk_scale_length_pc=5200.0, disk_scale_height_pc=600.0, bulge_scale_radius_pc=3160.0, bulge_amplitude=3.11,
        arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4,
    )
    radii = np.array([200.0, 1500.0, 3500.0, 9000.0, 20000.0])
    assert star_formation.recent_factor(2 * radii, big) == pytest.approx(star_formation.recent_factor(radii, SHAPE))


def test_the_law_exponent_sets_the_outer_fall(monkeypatch):
    """A steeper law (a larger exponent) makes the recent share fall faster outward."""
    star_formation._normalisation.cache_clear()
    low = star_formation.recent_factor(15_000.0, SHAPE) / star_formation.recent_factor(6_000.0, SHAPE)
    monkeypatch.setattr(tuning, "STAR_FORMATION_LAW_EXPONENT", 1.65)
    star_formation._normalisation.cache_clear()
    high = star_formation.recent_factor(15_000.0, SHAPE) / star_formation.recent_factor(6_000.0, SHAPE)
    star_formation._normalisation.cache_clear()
    assert high < low


@pytest.mark.parametrize("near, far", [(10.0, 100.0), (200.0, 900.0), (300.0, 6000.0), (4000.0, 5200.0),
                                       (6000.0, 30_000.0), (1.0, 2.0e6), (20_000.0, 90_000.0)])
def test_the_range_holds_every_value_of_the_interval(near, far):
    low, high = star_formation.recent_factor_range(near, far, SHAPE)
    values = star_formation.recent_factor(np.geomspace(near, far, 20_000), SHAPE)
    assert low <= values.min() * (1 + 1e-9)
    assert high >= values.max() * (1 - 1e-9)
    assert low >= values.min() * (1 - 1e-3)
    assert high <= values.max() * (1 + 1e-3)


def test_the_range_works_on_arrays_of_bands():
    near = np.array([10.0, 3000.0, 20_000.0])
    far = np.array([500.0, 7000.0, 40_000.0])
    low, high = star_formation.recent_factor_range(near, far, SHAPE)
    assert low.shape == high.shape == (3,)
    assert np.all(low <= high)


def test_the_populations_still_sum_to_the_density():
    """The factor moves stars between populations, never changes a cell's total."""
    for point in ((250.0, 0.0, 0.0), (3500.0, 800.0, 40.0), (8200.0, 0.0, 0.0), (15_000.0, 3000.0, 100.0)):
        total = sum(density.population_densities(point, SHAPE).values())
        assert total == pytest.approx(density.relative_density(point, SHAPE), rel=1e-9)


def test_the_array_form_matches_the_point_form():
    x = np.array([250.0, 1000.0, 3500.0, 8200.0, 15_000.0])
    y = np.array([0.0, 400.0, -800.0, 100.0, 3000.0])
    z = 25.0
    found = density.population_densities_array(x, y, z, SHAPE)
    for i in range(len(x)):
        one = density.population_densities((x[i], y[i], z), SHAPE)
        for name in ("young", "intermediate", "old", "bulge"):
            assert found[name][i] == pytest.approx(one[name], rel=1e-9)


def test_the_bulge_keeps_a_small_young_tail_from_the_disk_inside_it():
    """Stars under 5 Gyr in the bar's core: the thin disk's young and intermediate stars, and the part of the old
    disk's range under 5 Gyr, over all stars there (target: 3% from HST to 16 to 23% from microlensing)."""
    ages = tuning.STELLAR_POPULATION_AGE_RANGES_GY
    mix = density.population_densities((250.0, 0.0, 0.0), SHAPE)
    total = sum(max(v, 0.0) for v in mix.values())
    old_low, old_high = ages["old"]
    young = mix["young"] + mix["intermediate"] + mix["old"] * (5.0 - old_low) / (old_high - old_low)
    assert 0.03 <= young / total <= 0.25


def _outer_slope(shape=SHAPE):
    rate = star_formation.star_formation_rate_profile(RADII, shape)
    peak = RADII[np.argmax(rate)]
    window = (RADII > peak + 2000.0) & (RADII < peak + 8000.0)
    return np.polyfit(RADII[window] / 1000.0, np.log10(rate[window]), 1)[0]


def test_the_index_rising_where_the_gas_thins_steepens_the_outer_fall(monkeypatch):
    """GEN.450: with the index rising at low gas density the fall beyond the peak is steeper than the plain law's
    (-0.23 dex per kpc) and nearer the observed -0.28; the observed figure is a target, not an input."""
    with_rise = _outer_slope()
    _without_thin_gas_steepening(monkeypatch)
    plain = _outer_slope()
    star_formation._normalisation.cache_clear()
    assert plain == pytest.approx(-0.23, abs=0.01)
    assert with_rise < plain - 0.02
    assert with_rise == pytest.approx(-0.27, abs=0.03)


def test_the_thin_gas_factor_is_one_in_dense_gas_and_a_power_of_the_gas_where_it_thins():
    gas = np.exp(star_formation._log_gas(np.array([0.1, 8.0]), SHAPE)[0])
    factor = star_formation._thin_gas_factor(np.array([0.1, 8.0]), SHAPE)[0]
    saturation = tuning.STAR_FORMATION_SATURATION_GAS_MSUN_PC2
    power = tuning.STAR_FORMATION_LOW_DENSITY_EXPONENT - tuning.STAR_FORMATION_LAW_EXPONENT
    assert gas[0] > 2.0 * saturation and factor[0] > 0.75     # dense gas: the factor is near 1 (at twice the saturation it is (2/3)^0.7 = 0.75)
    assert factor[1] == pytest.approx((gas[1] / (gas[1] + saturation)) ** power)
    assert factor[1] < 0.3     # thin gas: the factor falls well below 1 (0.21 once the history keeps the thin gas the law gives, GEN.999)


def test_the_closed_form_log_slope_matches_the_slope_of_the_factor():
    radii = np.array([1500.0, 4000.0, 8000.0, 12_000.0, 20_000.0])
    step = 1.0
    numeric = (np.log(star_formation.recent_factor(radii + step, SHAPE))
               - np.log(star_formation.recent_factor(radii - step, SHAPE))) / (2.0 * step)
    assert star_formation._log_slope(radii, SHAPE) == pytest.approx(numeric, rel=2e-3, abs=1e-7)
