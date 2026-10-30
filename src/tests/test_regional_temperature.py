# tests/test_regional_temperature.py

"""
GEN.656: regional temperatures and the day and night swing (`physics/regional_temperature.py`), checked against
Earth, Mars, the Moon, Mercury and Venus, and for the rules that hold for any world.
"""

import pytest

from planetgen.physics import regional_temperature as rt

EARTH = dict(surface_temperature_k=288.0, atmospheric_pressure_pa=101_325.0, gravity_g=1.0, rotation_period_hours=24.0,
             axial_tilt_deg=23.4, period_years=1.0)
MARS = dict(surface_temperature_k=210.0, atmospheric_pressure_pa=610.0, gravity_g=0.38, rotation_period_hours=24.6,
            axial_tilt_deg=25.2, period_years=1.88)
MOON = dict(surface_temperature_k=250.0, atmospheric_pressure_pa=0.0, gravity_g=0.165, rotation_period_hours=655.7,
            axial_tilt_deg=1.5, period_years=0.0748)
MERCURY = dict(surface_temperature_k=440.0, atmospheric_pressure_pa=0.0, gravity_g=0.38, rotation_period_hours=1407.6,
               axial_tilt_deg=0.03, period_years=0.241)
VENUS = dict(surface_temperature_k=737.0, atmospheric_pressure_pa=9.2e6, gravity_g=0.9, rotation_period_hours=5832.0,
             axial_tilt_deg=177.0, period_years=0.615)
LATITUDE_FIELDS = ("temp_polar_k", "temp_north_k", "temp_equator_k", "temp_south_k", "temp_swing_k")
LOCKED_FIELDS = ("temp_day_side_k", "temp_terminator_k", "temp_night_side_k")


def test_the_band_areas_are_the_surface_shares():
    areas = rt.band_areas()
    assert areas["polar"] == pytest.approx(0.134, abs=0.001)
    assert areas["north"] == pytest.approx(0.262, abs=0.001) and areas["south"] == pytest.approx(0.262, abs=0.001)
    assert areas["equator"] == pytest.approx(0.342, abs=0.001)
    assert sum(areas.values()) == pytest.approx(1.0)


def test_a_sphere_averages_the_insolation_it_receives():
    insolation = rt.band_insolation(23.4)
    areas = rt.band_areas()
    assert sum(areas[name] * insolation[name] for name in areas) / sum(areas.values()) == pytest.approx(1.0, abs=0.03)
    # no tilt: the poles get the least, the equator the most; a flipped axis (90 degrees) puts the pole on top
    assert insolation["polar"] < insolation["north"] < insolation["equator"]
    assert rt.band_insolation(90.0)["polar"] > rt.band_insolation(90.0)["equator"]


def test_earth_has_poles_35_to_45_k_under_its_mean_and_a_warm_equator():
    out = rt.compute(EARTH)
    assert 288.0 - out["temp_polar_k"] == pytest.approx(38.0, abs=8.0)
    assert 8.0 < out["temp_equator_k"] - 288.0 < 20.0
    assert out["temp_north_k"] == out["temp_south_k"]
    assert all(out[name] is None for name in LOCKED_FIELDS)
    assert 0.5 < out["temp_swing_k"] < 15.0


def test_the_air_flattens_the_contrast():
    thin = rt.compute({**EARTH, "atmospheric_pressure_pa": 1_000.0})
    thick = rt.compute({**EARTH, "atmospheric_pressure_pa": 9.2e6})
    assert (thin["temp_equator_k"] - thin["temp_polar_k"]) > (thick["temp_equator_k"] - thick["temp_polar_k"]) > 0.0
    venus = rt.compute(VENUS)
    assert venus["temp_equator_k"] - venus["temp_polar_k"] < 60.0
    assert venus["temp_swing_k"] < 100.0


def test_mars_swings_by_tens_of_kelvin():
    out = rt.compute(MARS)
    assert 20.0 < out["temp_swing_k"] < 100.0
    assert out["temp_polar_k"] < MARS["surface_temperature_k"] < out["temp_equator_k"]


def test_the_moon_is_locked_with_a_hot_day_side_and_a_cold_night():
    out = rt.compute(MOON)
    assert all(out[name] is None for name in LATITUDE_FIELDS)
    assert 250.0 < out["temp_day_side_k"] < 320.0
    assert out["temp_night_side_k"] == rt.NIGHT_FLOOR_K
    assert out["temp_night_side_k"] < out["temp_terminator_k"] < out["temp_day_side_k"]


def test_a_thick_air_column_shares_the_heat_of_a_locked_world():
    locked = dict(EARTH, rotation_period_hours=0.1 * 8766.0, period_years=0.1, axial_tilt_deg=0.0)
    out = rt.compute(locked)
    assert out["temp_day_side_k"] - out["temp_night_side_k"] < 80.0
    assert out["temp_night_side_k"] > 200.0


def test_mercury_is_not_locked_and_swings_most_of_its_mean_twice_over():
    out = rt.compute(MERCURY)
    assert out["temp_swing_k"] == pytest.approx(1.414 * 440.0, rel=0.05)
    assert out["temp_equator_k"] > out["temp_north_k"] > out["temp_polar_k"]


def test_lock_is_a_rotation_period_equal_to_the_orbital_one():
    assert rt.is_locked(8766.0, 1.0) and rt.is_locked(8766.0 * 1.0005, 1.0)
    assert not rt.is_locked(24.0, 1.0) and not rt.is_locked(None, 1.0) and not rt.is_locked(24.0, None)


def test_a_body_without_a_temperature_or_a_spin_has_no_regions():
    assert not any(rt.compute({**EARTH, "surface_temperature_k": None}).values())
    assert not any(rt.compute({**EARTH, "rotation_period_hours": None}).values())


def test_the_bands_average_back_to_the_mean_before_the_air_moves_heat():
    out = rt.compute({**EARTH, "atmospheric_pressure_pa": 0.0})
    areas = rt.band_areas()
    mean = sum(areas[name] * out[f"temp_{name}_k"] for name in areas) / sum(areas.values())
    assert mean == pytest.approx(288.0, rel=0.002)


def _world(temperature_c):
    from planetgen.physics import habitability
    return habitability.World(pressure_kpa=101.3, temperature_c=temperature_c, gases={"O2": 0.21})


def test_the_sentence_names_the_equator_the_poles_and_the_swing():
    from planetgen.physics import protection
    text = protection.regional_sentence(_world(15.0), rt.compute(EARTH))
    assert text.startswith("The equator averages 28.7 C and the poles -17.9 C.")
    assert "warm clothing is recommended" in text and text.count("recommended") == 1
    assert "Day and night" not in text  # an Earth-like swing of a few kelvin is not worth a sentence


def test_the_sentence_of_a_locked_world_names_its_three_zones():
    from planetgen.physics import protection
    locked = dict(EARTH, rotation_period_hours=0.1 * 8766.0, period_years=0.1, axial_tilt_deg=0.0)
    text = protection.regional_sentence(_world(7.0), rt.compute(locked))
    assert "day side" in text and "terminator" in text and "night side" in text


def test_a_body_with_no_regions_has_no_sentence_and_no_summary():
    from planetgen.physics import protection
    assert protection.regional_sentence(_world(15.0), None) == ""
    assert protection.regional_sentence(_world(15.0), {name: None for name in rt.FIELDS}) == ""
    assert rt.summary({}) is None
    assert rt.summary(rt.compute(MOON)) == {"day_side": pytest.approx(rt.compute(MOON)["temp_day_side_k"]),
                                            "terminator": pytest.approx(rt.compute(MOON)["temp_terminator_k"]),
                                            "night_side": rt.NIGHT_FLOOR_K,
                                            # GEN.848: a locked world has no day, so its places swing by nothing
                                            "day_side_swing": 0.0, "terminator_swing": 0.0, "night_side_swing": 0.0,
                                            "mean_swing": 0.0}


def test_the_daily_range_has_the_stored_value_exactly_in_the_middle():
    """GEN.795, GEN.796: the low and high are the value minus and plus half the peak-to-peak swing, in kelvin."""
    low, high = rt.daily_range(300.0, 9.07)
    assert (low, high) == pytest.approx((295.465, 304.535)) and (low + high) / 2.0 == pytest.approx(300.0)
    assert rt.daily_range(None, 9.0) is None and rt.daily_range(300.0, None) is None and rt.daily_range(300.0, 0.0) is None


def test_the_sentence_gives_the_swing_as_plus_or_minus_and_the_equators_low_and_high():
    from planetgen.physics import protection
    regions = {"temp_polar_k": 250.0, "temp_equator_k": 300.0, "temp_swing_k": 20.0}
    text = protection.regional_sentence(_world(15.0), regions)
    assert "from the coldest hour to the hottest (±10 K)" in text
    assert "so the equator runs from 16.9 C to 36.9 C over a day." in text


BAND_SWINGS = ("temp_polar_swing_k", "temp_north_swing_k", "temp_south_swing_k")
ZONE_SWINGS = ("temp_day_side_swing_k", "temp_terminator_swing_k", "temp_night_side_swing_k")


def test_gen848_every_band_of_a_turning_world_has_a_swing_that_falls_toward_the_pole():
    out = rt.compute(EARTH)
    assert out["temp_swing_k"] > out["temp_north_swing_k"] > out["temp_polar_swing_k"] > 0.0
    assert out["temp_north_swing_k"] == out["temp_south_swing_k"]
    assert all(out[name] is None for name in ZONE_SWINGS)
    # the world's mean swing is the area-weighted mean of its bands: between the polar and the equatorial swing
    assert out["temp_polar_swing_k"] < out["temp_mean_swing_k"] < out["temp_swing_k"]
    areas = rt.band_areas()
    expected = sum(areas[name] * out[column] for name, column in
                   (("polar", "temp_polar_swing_k"), ("north", "temp_north_swing_k"), ("equator", "temp_swing_k"),
                    ("south", "temp_south_swing_k"))) / sum(areas.values())
    assert out["temp_mean_swing_k"] == pytest.approx(expected)


def test_gen848_the_bands_swings_are_the_equators_cut_by_the_cosine_of_the_latitude():
    ratios = rt._band_swing_ratios()
    assert ratios["equator"] == 1.0
    assert 0.2 < ratios["polar"] < 0.5 < ratios["north"] < 0.95 and ratios["north"] == ratios["south"]
    out = rt.compute(MARS)
    assert out["temp_polar_swing_k"] == pytest.approx(out["temp_swing_k"] * ratios["polar"])


def test_gen848_a_locked_world_has_no_day_so_its_places_swing_by_nothing():
    out = rt.compute(MOON)
    assert all(out[name] == 0.0 for name in ZONE_SWINGS) and out["temp_mean_swing_k"] == 0.0
    assert all(out[name] is None for name in BAND_SWINGS) and out["temp_swing_k"] is None


def test_gen848_every_stored_temperature_names_the_column_of_its_swing():
    assert set(rt.SWING_OF.values()) <= set(rt.FIELDS)
    assert set(rt.RANGE_FIELDS) == set(rt.SWING_OF) and "surface_temperature_k" in rt.RANGE_FIELDS
    assert rt.SWING_OF["temp_equator_k"] == "temp_swing_k"
    assert len(set(rt.SWING_OF.values())) == len(rt.SWING_OF)


def test_gen848_no_air_and_no_rotation_leave_the_swings_empty():
    assert all(value is None for value in rt.compute(dict(EARTH, rotation_period_hours=0.0)).values())
