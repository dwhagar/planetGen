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
    assert equator.longest_day_days == pytest.approx(0.5, abs=0.05) and equator.longest_night_days == pytest.approx(0.5, abs=0.05) and equator.state == lb.CYCLE      # half a solar day
    assert pole.longest_day_days == pytest.approx(184.0, abs=3.0)       # 186 with refraction, which the model leaves out
    assert pole.longest_night_days == pytest.approx(181.0, abs=3.0)       # the star's disk keeps the day the longer; 179 for a point star
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


def test_the_angular_radius_of_the_star_comes_back_from_the_stored_share_of_time_it_is_up():
    """GEN.1108: the stored mean day share is `(1 + sin alpha) / 2`, so the bands of a stored body see the star the generator saw."""
    assert lb.alpha_deg_of({}) == 0.0 and lb.alpha_deg_of({"day_share_mean": 0.5}) == 0.0
    share = (1.0 + math.sin(math.radians(ALPHA_EARTH_DEG))) / 2.0
    assert lb.alpha_deg_of({"day_share_mean": share}) == pytest.approx(ALPHA_EARTH_DEG)
    assert lb.RING_STATES == (lb.DAY_SIDE, lb.TERMINATOR, lb.NIGHT_SIDE)


RING = dict(ring_inner_km=1.2 * 6371.0, ring_outer_km=2.3 * 6371.0, ring_optical_depth=1.0, radius_km=6371.0)


def test_a_ring_shades_the_bands_behind_it_and_leaves_the_poles_alone():
    """GEN.1055 R8: a tilted world with a broad ring keeps less sunlight in its tropic and equator bands; a pole sees the Sun above the ring."""
    shade = lb.ring_shade(dict(EARTH, **RING), 26.7)
    assert len(shade) == 7 and shade == pytest.approx(list(reversed(shade)))
    assert shade[0] == pytest.approx(1.0) and shade[-1] == pytest.approx(1.0)
    assert shade[3] < shade[2] < shade[1] <= 1.0 and 0.6 < shade[3] < 0.9


def test_the_ring_cools_the_shaded_bands_and_the_mean_of_the_world_stays():
    plain = lb.compute(dict(EARTH, axial_tilt_deg=26.7), ALPHA_EARTH_DEG)
    ringed = lb.compute(dict(EARTH, axial_tilt_deg=26.7, **RING), ALPHA_EARTH_DEG)
    assert ringed[3].mean_k < plain[3].mean_k and ringed[0].mean_k > plain[0].mean_k
    mean = lambda bands: sum(b.area_share * b.mean_k for b in bands)
    assert mean(ringed) == pytest.approx(mean(plain), rel=1e-6)


def test_a_thicker_ring_shades_more_and_no_ring_or_no_tilt_shades_nothing():
    thin = lb.ring_shade(dict(EARTH, **dict(RING, ring_optical_depth=0.1)), 26.7)
    thick = lb.ring_shade(dict(EARTH, **RING), 26.7)
    assert thick[3] < thin[3] < 1.0
    assert lb.ring_shade(EARTH, 26.7) is None
    assert lb.ring_shade(dict(EARTH, **dict(RING, ring_optical_depth=0.0)), 26.7) is None
    assert lb.ring_shade(dict(EARTH, **RING), 0.0) is None
    # a world that turns once an orbit has rings of its own (the day side to the night side), and its ring is not asked for
    assert [b.mean_k for b in lb.compute(dict(LOCKED, **RING), 0.0)] == [b.mean_k for b in lb.compute(LOCKED, 0.0)]


def test_the_climate_values_of_a_body_carry_its_ring():
    from planetgen.physics import habitability_world

    class Body:
        ring_inner_km, ring_outer_km, ring_optical_depth = 8000.0, 14000.0, 0.5
    values = habitability_world.values_from_body(Body())
    assert (values["ring_inner_km"], values["ring_outer_km"], values["ring_optical_depth"]) == (8000.0, 14000.0, 0.5)


def longest_spells(latitude, tilt, alpha, orbital_hours, solar_day_hours, **options):
    """`(longest day, longest night, step)` at one latitude from `latitude_bands.spells`."""
    found = lb.spells([latitude], tilt, alpha, orbital_hours, solar_day_hours, **options)
    return float(found.longest_day[0]), float(found.longest_night[0]), found.step


# --- GEN.1260: the longest unbroken day and night of a slow spinner -------------------------------------------------------------

NITRIKIN_LIKE = dict(surface_temperature_k=280.0, atmospheric_pressure_pa=101325.0, gravity_g=1.0, rotation_period_hours=27.9 * 24.0, axial_tilt_deg=127.0,
                     period_years=28.5 / 365.25)
"""dict: A body that turns backwards in 27.9 days round a 28.5-day orbit, so its solar day is 14.1 days (the shape of Nitrikin III, which is a target and never an input)."""


def test_a_slow_spinner_has_an_equator_day_of_half_a_solar_day_and_a_pole_day_of_half_a_year():
    orbital, solar = 28.5 * 24.0, 24.0 / (1.0 / 27.9 + 1.0 / 28.5)
    assert solar == pytest.approx(14.1 * 24.0, rel=0.002)
    day, night, _ = longest_spells(0.0, 127.0, 0.0, orbital, solar)
    assert day == pytest.approx(7.05, abs=0.2) and night == pytest.approx(7.05, abs=0.2)       # about 7.1 days, not 0
    day, night, _ = longest_spells(25.7, 127.0, 0.0, orbital, solar)
    assert day == pytest.approx(8.4, abs=0.3) and night == pytest.approx(8.4, abs=0.3)
    day, night, _ = longest_spells(38.6, 127.0, 0.0, orbital, solar)
    assert day == pytest.approx(17.0, abs=1.5) and night == pytest.approx(17.0, abs=1.5) and day > solar / 24.0      # the star stays up past a solar day
    pole_day, pole_night, _ = longest_spells(90.0, 127.0, 0.0, orbital, solar)
    assert pole_day == pytest.approx(14.25, abs=0.2) and pole_night == pytest.approx(14.25, abs=0.2)       # half the 28.5-day year


def test_the_bands_of_a_slow_spinner_show_the_spell_and_not_the_polar_day_alone():
    bands = {b.name: b for b in lb.compute(NITRIKIN_LIKE, 0.0)}
    solar_days = 14.1
    equator = bands["equator"]
    assert 7.0 < equator.longest_day_days < 8.2 and 7.0 < equator.longest_night_days < 8.2      # read at the band's edge, 12.9 degrees
    assert equator.state == lb.CYCLE
    assert bands["north_tropic"].longest_day_days > solar_days and bands["north_tropic"].state == lb.POLAR       # its poleward edge, 38.6 degrees, sees a day past the solar day
    assert bands["north_middle"].longest_day_days > solar_days and bands["north_middle"].state == lb.POLAR
    assert bands["north_polar"].longest_day_days == pytest.approx(14.25, abs=0.2) and bands["north_polar"].longest_night_days == pytest.approx(14.25, abs=0.2)
    assert all(b.longest_day_days > 0.0 and b.longest_night_days > 0.0 for b in bands.values())


def _closed_form_polar_days(latitude_deg, tilt_deg, alpha, orbital_days):
    """The textbook polar day and night for a fast spinner: the star is up all day while its declination is at least `90 - |latitude| - alpha` (an oracle, not the search)."""
    sin_tilt = math.sin(math.radians(tilt_deg))
    day_ratio = math.cos(math.radians(latitude_deg) + alpha) / sin_tilt
    night_ratio = math.cos(math.radians(latitude_deg) - alpha) / sin_tilt
    day = 0.0 if day_ratio >= 1.0 else (math.pi - 2.0 * math.asin(day_ratio)) / (2.0 * math.pi)
    night = 0.0 if night_ratio >= 1.0 else (math.pi - 2.0 * math.asin(night_ratio)) / (2.0 * math.pi)
    return day * orbital_days, night * orbital_days


def test_the_search_finds_the_textbook_polar_day_and_night_without_a_case_for_it():
    alpha = math.radians(0.267)
    for latitude in (70.0, 80.0, 90.0):
        day, night, step = longest_spells(latitude, 23.44, alpha, 8766.0, 24.0)
        expected = _closed_form_polar_days(latitude, 23.44, alpha, 365.25)
        assert day == pytest.approx(expected[0], abs=2.0) and night == pytest.approx(expected[1], abs=2.0)       # the search adds the day before and the day after the polar one: the sun is up from the sunrise before it to the sunset after


def test_a_body_with_a_year_of_many_thousand_days_is_searched_by_the_same_function():
    alpha = math.radians(0.267)
    day, night, step = longest_spells(90.0, 97.77, 0.0, 84.0 * 8766.0, 17.24)
    assert step > 17.24 / 24.0 / lb._MIN_SAMPLES_PER_DAY      # the steps hold whole solar days, so the search takes the average through a step
    assert day == pytest.approx(84.0 * 365.25 / 2.0, rel=0.01) and night == pytest.approx(84.0 * 365.25 / 2.0, rel=0.01)
    equator = longest_spells(12.9, 97.77, alpha, 84.0 * 8766.0, 17.24)
    assert equator[0] > 365.0       # even the equator has a day of years at that tilt
    short = longest_spells(0.0, 23.44, alpha, 84.0 * 8766.0, 17.24)
    assert short[0] == pytest.approx(17.24 / 48.0, rel=0.1)       # a tilt of 23 degrees leaves the equator its half solar day


def test_the_day_and_night_of_a_band_are_the_mean_length_of_a_stretch_not_half_a_solar_day_at_every_latitude():
    bands = {b.name: b for b in lb.compute(NITRIKIN_LIKE, 0.0)}
    equator = bands["equator"]
    assert equator.day_hours == pytest.approx(14.1 * 24.0 / 2.0, rel=0.03) and equator.night_hours == pytest.approx(14.1 * 24.0 / 2.0, rel=0.03)      # about 7 days
    polar = bands["north_polar"]
    assert polar.day_hours > 1.2 * equator.day_hours and polar.night_hours > 1.2 * equator.night_hours      # the long stretches near the pole pull the mean up
    earth = {b.name: b for b in lb.compute(EARTH, ALPHA_EARTH_DEG)}
    assert earth["north_polar"].day_hours > 20.0 and earth["equator"].day_hours == pytest.approx(12.0, abs=0.15)      # a polar day of months is one stretch


def test_the_mean_is_the_time_up_over_the_number_of_stretches():
    orbital, solar = 28.5 * 24.0, 24.0 / (1.0 / 27.9 + 1.0 / 28.5)
    found = lb.spells([0.0, 25.7, 90.0], 127.0, 0.0, orbital, solar)
    assert list(found.mean_day) == pytest.approx([7.05, 6.9, 14.25], abs=0.4) and list(found.mean_night) == pytest.approx([7.05, 7.0, 14.25], abs=0.4)
    assert all(found.mean_day <= found.longest_day + 1e-9) and all(found.mean_night <= found.longest_night + 1e-9)
