# tests/test_seasonal_swing.py

"""
GEN.1039: the seasonal swing of the temperature from the tilt and the orbit (`planetgen.physics.seasonal_swing`): the light of the year (its sphere mean, its eccentric
orbit), the annual means, the Earth, a Mars-like dry thin-air planet, a Uranus-like tilt, the linear limit against the closed-form harmonic, the freezing and wet-bulb
helpers, and the inputs refused.
"""

import math

import numpy as np
import pytest

from planetgen.physics import habitability, regional_temperature
from planetgen.physics import seasonal_swing as ss

EARTH = {"eccentricity": 0.0167, "mean_k": 288.0, "stellar_flux_w_m2": 1361.0, "bond_albedo": 0.30, "pressure_pa": 1.013e5, "gravity_ms2": 9.81,
         "period_years": 1.0, "perihelion_longitude_deg": 283.0}
MARS = {"mean_k": 210.0, "stellar_flux_w_m2": 590.0, "bond_albedo": 0.25, "pressure_pa": 610.0, "gravity_ms2": 3.71, "ocean_fraction": 0.0,
        "period_years": 1.881, "perihelion_longitude_deg": 251.0}
SIGMA = 5.670374419e-8


def _equal_area_latitudes(n=36):
    return np.degrees(np.arcsin(-1.0 + (np.arange(n) + 0.5) / n * 2.0))


def test_the_time_mean_of_the_light_on_an_eccentric_orbit_is_the_published_one():
    for e in [0.0, 0.0167, 0.0934, 0.3, 0.6]:
        sin_longitude, inverse_square = ss.orbit_year(e, 40.0)
        assert inverse_square.mean() == pytest.approx(1.0 / math.sqrt(1.0 - e * e), rel=2e-4)
        assert np.all(np.abs(sin_longitude) <= 1.0)


def test_a_planet_is_closest_to_its_star_at_the_start_of_the_year_and_the_light_there_is_highest():
    _, inverse_square = ss.orbit_year(0.0934)
    assert int(np.argmax(inverse_square)) in (0, ss.STEPS - 1)
    assert inverse_square.max() == pytest.approx(1.0 / (1.0 - 0.0934) ** 2, rel=2e-3)
    assert inverse_square.min() == pytest.approx(1.0 / (1.0 + 0.0934) ** 2, rel=2e-3)


def test_the_light_of_the_year_averages_to_a_quarter_of_the_flux_over_the_sphere_at_any_tilt_and_eccentricity():
    latitudes = _equal_area_latitudes()
    for tilt in [0.0, 23.44, 54.7, 90.0, 98.0, 135.0, 180.0]:
        for e in [0.0, 0.0934, 0.5]:
            mean = ss.daily_insolation_year(latitudes, tilt, e, 70.0).mean()
            assert mean == pytest.approx(ss.sphere_mean_insolation(e), rel=2e-3), (tilt, e)


def test_the_seasons_of_a_tilt_and_of_its_supplement_are_the_same():
    a = ss.daily_insolation_year([10.0, 60.0], 98.0, 0.1, 30.0)
    b = ss.daily_insolation_year([10.0, 60.0], 82.0, 0.1, 30.0)
    assert np.allclose(a, b)


def test_an_untilted_round_orbit_has_no_seasons_and_an_eccentric_one_has_one():
    flat = ss.seasonal_swing([0.0, 45.0], 0.0, 0.0, 288.0, 1361.0, 0.3, 1.013e5, 9.81, 0.7, 1.0)
    assert all(s.swing_k < 1e-6 for s in flat)
    eccentric = ss.seasonal_swing([0.0], 0.0, 0.0934, 288.0, 1361.0, 0.3, 1.013e5, 9.81, 0.0, 1.0)[0]
    assert eccentric.swing_k > 1.0


def test_the_annual_means_average_to_the_worlds_mean_and_a_high_tilt_warms_the_pole_over_the_equator():
    latitudes = _equal_area_latitudes()
    for tilt in [0.0, 23.44, 98.0]:
        means = ss.annual_mean_temperatures(latitudes, tilt, 0.0, 288.0, 0.05)
        assert means.mean() == pytest.approx(288.0, abs=0.3)
    flat = ss.annual_mean_temperatures([0.0, 90.0], 0.0, 0.0, 288.0, 0.05)
    tilted = ss.annual_mean_temperatures([0.0, 90.0], 98.0, 0.0, 288.0, 0.05)
    assert flat[1] < flat[0] and tilted[1] > tilted[0]


def test_the_earth_has_small_tropical_seasons_and_mid_latitude_seasons_of_the_published_order():
    tropics = ss.seasonal_swing([0.0], 23.44, ocean_fraction=0.71, **EARTH)[0]
    assert tropics.swing_k < 3.0
    mid = ss.seasonal_swing([60.0, -60.0], 23.44, ocean_fraction=0.4, **EARTH)
    # the zonal mean seasonal range at 60 degrees on the Earth is about 25 K (the simple model gives about 17)
    assert all(10.0 < s.swing_k < 40.0 for s in mid)
    assert mid[0].swing_k > 5.0 * tropics.swing_k


def test_the_hemispheres_have_opposite_seasons():
    north, south = ss.seasonal_swing([60.0, -60.0], 23.44, ocean_fraction=0.4, **EARTH)
    offset = abs(north.warmest_fraction - south.warmest_fraction)
    assert min(offset, 1.0 - offset) == pytest.approx(0.5, abs=0.08)


def test_a_sea_damps_the_seasons_and_dry_ground_does_not():
    kwargs = dict(EARTH)
    swings = [ss.seasonal_swing([60.0], 23.44, ocean_fraction=f, **kwargs)[0].swing_k for f in (0.0, 0.4, 0.8, 1.0)]
    assert swings == sorted(swings, reverse=True) and swings[0] > 3.0 * swings[-1]


def test_a_deep_mixed_layer_damps_and_delays_the_seasons_more():
    shallow = ss.seasonal_swing([60.0], 23.44, ocean_fraction=1.0, mixed_layer_m=20.0, **EARTH)[0]
    deep = ss.seasonal_swing([60.0], 23.44, ocean_fraction=1.0, mixed_layer_m=200.0, **EARTH)[0]
    assert deep.swing_k < shallow.swing_k and deep.warmest_fraction > shallow.warmest_fraction


def test_a_dry_thin_air_planet_follows_the_light_as_the_quasi_static_fourth_root():
    for e in (0.0, 0.0934):
        seasons = ss.seasonal_swing([0.0, 45.0], 25.2, e, **MARS)
        for season in seasons:
            q = ss.daily_insolation_year([season.latitude_deg], 25.2, e, 251.0)[0]
            static = season.mean_k * (q / q.mean()) ** 0.25
            assert season.swing_k == pytest.approx(float(static.max() - static.min()), rel=0.08), (e, season.latitude_deg)


def test_mars_eccentric_orbit_adds_to_the_equators_season_and_its_swing_is_the_published_order():
    round_orbit = ss.seasonal_swing([0.0], 25.2, 0.0, **MARS)[0]
    mars = ss.seasonal_swing([0.0], 25.2, 0.0934, **MARS)[0]
    assert mars.swing_k > 3.0 * round_orbit.swing_k
    assert 15.0 < mars.swing_k < 35.0           # the equator's surface temperature varies about 20 to 30 K between the seasons


def test_a_uranus_like_tilt_puts_the_hottest_summer_at_the_pole():
    wet = ss.seasonal_swing([0.0, 45.0, 90.0], 98.0, 0.0, 288.0, 1361.0, 0.3, 1.013e5, 9.81, 0.71, 1.0)
    equator, _, pole = wet
    assert pole.mean_k > equator.mean_k
    assert pole.max_k > equator.max_k and pole.swing_k > equator.swing_k
    dry = ss.seasonal_swing([0.0, 90.0], 98.0, 0.0, 288.0, 1361.0, 0.3, 1.013e5, 9.81, 0.0, 1.0)
    assert dry[1].max_k > dry[0].max_k + 50.0 and dry[1].min_k < dry[1].mean_k - 50.0
    # the pole has a long dark winter: well under freezing for a good part of the year on a dry world
    assert dry[1].warm_fraction < 0.8


def test_a_small_tilt_is_the_closed_form_harmonic_response():
    for fraction in (1.0, 0.3):
        latitude = [45.0]
        season = ss.seasonal_swing(latitude, 3.0, 0.0, 288.0, 1361.0, 0.3, 1.013e5, 9.81, fraction, 1.0)[0]
        column = regional_temperature.column_in_earths(1.013e5, 9.81)
        share = regional_temperature.transport_share(column)
        capacity = ss.heat_capacity_j_m2_k(1.013e5, 9.81, fraction, 1.0)
        emissivity = 0.7 * 1361.0 / 4.0 / (SIGMA * 288.0 ** 4)
        mean = ss.annual_mean_temperatures(latitude, 3.0, 0.0, 288.0, share)[0]
        slope = 4.0 * emissivity * SIGMA * mean ** 3 / (1.0 - share)
        insolation = ss.daily_insolation_year(latitude, 3.0, 0.0)[0]
        forcing = 0.7 * 1361.0 * (insolation - insolation.mean())
        first = np.fft.rfft(forcing)[1] / ss.STEPS * 2.0
        omega = 2.0 * math.pi / (365.25 * 86400.0)
        amplitude = abs(first) / abs(slope + 1j * capacity * omega)
        assert season.swing_k == pytest.approx(2.0 * amplitude, rel=0.03), fraction


def test_the_year_closes_and_its_mean_is_the_places_mean():
    for season in ss.seasonal_swing([0.0, 45.0, 80.0], 40.0, 0.2, 280.0, 1000.0, 0.3, 5.0e4, 9.0, 0.5, 1.2):
        assert len(season.series_k) == ss.STEPS
        assert sum(season.series_k) / ss.STEPS == pytest.approx(season.mean_k, abs=1e-6)
        assert season.min_k == pytest.approx(min(season.series_k)) and season.max_k == pytest.approx(max(season.series_k))
        assert 0.0 <= season.warm_fraction <= 1.0 and 0.0 < season.warmest_fraction < 1.0


def test_a_thick_air_carries_the_seasons_away():
    thin = ss.seasonal_swing([70.0], 50.0, 0.0, 288.0, 1361.0, 0.3, 1.0e5, 9.81, 0.0, 1.0)[0]
    thick = ss.seasonal_swing([70.0], 50.0, 0.0, 288.0, 1361.0, 0.3, 9.0e6, 9.81, 0.0, 1.0)[0]
    assert thick.swing_k < 0.3 * thin.swing_k


def test_the_freezing_share_of_the_year_follows_the_temperature():
    cold = ss.seasonal_swing([60.0], 23.44, ocean_fraction=0.4, **{**EARTH, "mean_k": 250.0})[0]
    mild = ss.seasonal_swing([60.0], 23.44, ocean_fraction=0.4, **EARTH)[0]
    tropics = ss.seasonal_swing([0.0], 23.44, ocean_fraction=0.4, **EARTH)[0]
    assert cold.warm_fraction < mild.warm_fraction < tropics.warm_fraction
    assert tropics.warm_fraction == 1.0 and cold.warm_fraction < 0.1


def test_the_peak_wet_bulb_is_the_stull_fit_at_the_warmest_hour_and_rises_with_the_humidity():
    season = ss.seasonal_swing([20.0], 23.44, ocean_fraction=0.0, **{**EARTH, "mean_k": 290.0})[0]
    assert season.max_k - ss.FREEZING_K < 50.0
    dry, wet = ss.peak_wet_bulb_c(season, 0.2), ss.peak_wet_bulb_c(season, 0.9)
    assert dry < wet < season.max_k - ss.FREEZING_K
    assert wet == pytest.approx(habitability.wet_bulb_c(season.max_k - ss.FREEZING_K, 0.9))


def test_inputs_the_model_cannot_use_are_refused():
    base = (23.44, 0.0167, 288.0, 1361.0, 0.3, 1.013e5, 9.81, 0.7, 1.0)
    with pytest.raises(ss.SeasonalSwingError):
        ss.seasonal_swing([], *base)
    with pytest.raises(ss.SeasonalSwingError):
        ss.seasonal_swing([95.0], *base)
    with pytest.raises(ss.SeasonalSwingError):
        ss.seasonal_swing([0.0], 200.0, *base[1:])
    with pytest.raises(ss.SeasonalSwingError):
        ss.seasonal_swing([0.0], 23.44, 1.0, *base[2:])
    with pytest.raises(ss.SeasonalSwingError):
        ss.seasonal_swing([0.0], 23.44, 0.1, 288.0, 1361.0, 1.0, 1.013e5, 9.81, 0.7, 1.0)
    with pytest.raises(ss.SeasonalSwingError):
        ss.seasonal_swing([0.0], 23.44, 0.1, 288.0, 1361.0, 0.3, 1.013e5, 9.81, 1.5, 1.0)
