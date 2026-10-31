# tests/test_latitude_bands.py

"""
GEN.1108 part A: the climate of a body in seven latitude bands (`physics/latitude_bands.py`), a pure function of stored columns. Earth, Uranus, Venus and a locked planet
are test targets, never inputs.
"""

import math

import pytest

from planetgen.physics import latitude_bands as lb

EARTH = dict(surface_temperature_k=288.0, atmospheric_pressure_pa=101325.0, gravity_g=1.0, rotation_period_hours=23.934, axial_tilt_deg=23.44, period_years=1.0)
URANUS_LIKE = dict(EARTH, axial_tilt_deg=97.77, rotation_period_hours=17.24, period_years=84.0)
VENUS = dict(surface_temperature_k=737.0, atmospheric_pressure_pa=9.2e6, gravity_g=0.9, rotation_period_hours=5832.5, axial_tilt_deg=177.4, period_years=0.6152)
LOCKED = dict(surface_temperature_k=270.0, atmospheric_pressure_pa=50000.0, gravity_g=1.0, rotation_period_hours=8766.0 * 0.05, axial_tilt_deg=0.0, period_years=0.05)
ALPHA_EARTH_DEG = 0.267


def test_the_seven_bands_tile_the_sphere_with_the_equator_band_in_the_middle():
    bands = lb.compute(EARTH, ALPHA_EARTH_DEG)
    assert [b.name for b in bands] == list(lb.NAMES) and len(bands) == 7
    assert bands[0].from_deg == 90.0 and bands[-1].to_deg == -90.0
    assert bands[3].from_deg == pytest.approx(-bands[3].to_deg) == pytest.approx(lb.BAND_WIDTH_DEG / 2.0)
    for upper, lower in zip(bands, bands[1:]):
        assert upper.to_deg == pytest.approx(lower.from_deg)
    assert sum(b.area_share for b in bands) == pytest.approx(1.0)
    assert [b.area_share for b in bands] == pytest.approx([b.area_share for b in reversed(bands)])


def test_earth_has_twelve_hours_of_day_a_polar_day_of_about_half_a_year_and_cold_poles():
    bands = {b.name: b for b in lb.compute(EARTH, ALPHA_EARTH_DEG)}
    equator, pole = bands["equator"], bands["north_polar"]
    assert equator.day_hours == pytest.approx(12.0, abs=0.15) and equator.night_hours == pytest.approx(12.0, abs=0.15)
    assert equator.longest_day_days == 0.0 and equator.longest_night_days == 0.0 and equator.state == lb.CYCLE
    assert pole.longest_day_days == pytest.approx(184.0, abs=3.0)       # 186 with refraction, which the model leaves out
    assert pole.longest_night_days == pytest.approx(179.0, abs=3.0)
    assert pole.state == lb.POLAR
    assert pole.mean_k < bands["north_middle"].mean_k < bands["north_tropic"].mean_k < equator.mean_k
    assert equator.mean_k - pole.mean_k == pytest.approx(50.0, abs=15.0)
    assert bands["north_polar"].mean_k == pytest.approx(bands["south_polar"].mean_k)      # a circular orbit


def test_the_day_and_night_temperatures_straddle_the_mean_and_the_swing_falls_toward_the_poles():
    bands = lb.compute(EARTH, ALPHA_EARTH_DEG)
    for band in bands:
        assert band.night_k < band.mean_k < band.day_k
        assert band.day_k - band.night_k == pytest.approx(band.swing_k)
    swings = [b.swing_k for b in bands[:4]]
    assert swings == sorted(swings)


def test_uranus_has_polar_days_and_nights_of_years_in_every_band():
    bands = lb.compute(URANUS_LIKE, 0.0)
    assert {b.state for b in bands} == {lb.POLAR}
    assert bands[0].longest_day_days == pytest.approx(84.0 * 365.25 / 2.0, rel=0.01)       # the pole: half the 84-year orbit
    assert bands[3].longest_day_days > 365.0                                                  # even the equator band
    assert bands[0].mean_k > bands[3].mean_k       # the poles take more of the year's light than the equator at a tilt past 90 degrees


def test_venus_turns_slowly_backwards_with_a_solar_day_of_117_days():
    bands = {b.name: b for b in lb.compute(VENUS, 0.5)}
    equator = bands["equator"]
    assert equator.day_hours + equator.night_hours == pytest.approx(116.75 * 24.0, rel=0.01)
    assert equator.state == lb.CYCLE
    assert equator.mean_k > bands["north_polar"].mean_k > 650.0      # a thick air column barely lets the poles cool


def test_a_locked_planet_has_rings_from_a_day_side_to_a_night_side():
    bands = lb.compute(LOCKED, 0.5)
    assert [b.name for b in bands] == list(lb.RING_NAMES)
    assert [b.state for b in bands] == [lb.DAY_SIDE] * 3 + [lb.TERMINATOR] + [lb.NIGHT_SIDE] * 3
    assert sum(b.area_share for b in bands) == pytest.approx(1.0)
    temperatures = [b.mean_k for b in bands]
    assert temperatures == sorted(temperatures, reverse=True)
    assert all(b.day_hours is None and b.longest_day_days is None and b.swing_k == 0.0 for b in bands)
    assert min(temperatures) >= 20.0
    assert sum(b.area_share * b.mean_k for b in bands) == pytest.approx(LOCKED["surface_temperature_k"], rel=0.02)


def test_the_air_carries_heat_so_a_thin_air_world_has_more_contrast_than_a_thick_one():
    thin = lb.compute(dict(EARTH, atmospheric_pressure_pa=1000.0), ALPHA_EARTH_DEG)
    thick = lb.compute(dict(EARTH, atmospheric_pressure_pa=9.0e6), ALPHA_EARTH_DEG)
    assert thin[3].mean_k - thin[0].mean_k > thick[3].mean_k - thick[0].mean_k


def test_the_mean_of_the_bands_is_the_mean_of_the_world():
    bands = lb.compute(EARTH, ALPHA_EARTH_DEG)
    assert sum(b.area_share * b.mean_k for b in bands) == pytest.approx(EARTH["surface_temperature_k"], rel=0.01)


def test_a_body_without_a_temperature_or_a_spin_has_no_bands():
    assert lb.compute(dict(EARTH, surface_temperature_k=None)) == []
    assert lb.compute(dict(EARTH, rotation_period_hours=None)) == []


def test_the_polar_day_share_follows_the_tilt():
    assert lb.polar_day_share(60.0, 0.0, 0.0) == 0.0
    assert lb.polar_day_share(90.0, 23.44, 0.0) == pytest.approx(0.5)
    assert lb.polar_day_share(60.0, 23.44, 0.0) == 0.0       # south of the Arctic circle at 66.56
    assert lb.polar_day_share(70.0, 23.44, 0.0) == pytest.approx(0.17, abs=0.02)
    assert lb.polar_day_share(80.0, 23.44, 0.0) == pytest.approx(0.35, abs=0.03)
    assert lb.polar_night_share(90.0, 23.44, 0.0) == pytest.approx(0.5)
    assert lb.polar_day_share(90.0, 23.44, 0.00466) > 0.5 > lb.polar_night_share(90.0, 23.44, 0.00466)


def test_the_angular_radius_of_the_star_comes_back_from_the_stored_share_of_time_it_is_up():
    """GEN.1108: the stored mean day share is `(1 + sin alpha) / 2`, so the bands of a stored body see the star the generator saw."""
    assert lb.alpha_deg_of({}) == 0.0 and lb.alpha_deg_of({"day_share_mean": 0.5}) == 0.0
    share = (1.0 + math.sin(math.radians(ALPHA_EARTH_DEG))) / 2.0
    assert lb.alpha_deg_of({"day_share_mean": share}) == pytest.approx(ALPHA_EARTH_DEG)
    assert lb.RING_STATES == (lb.DAY_SIDE, lb.TERMINATOR, lb.NIGHT_SIDE)
