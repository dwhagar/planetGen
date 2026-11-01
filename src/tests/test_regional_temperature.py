# tests/test_regional_temperature.py

"""
GEN.656, GEN.1108: the climate of a body in seven bands of latitude and the day and night swing (`physics/regional_temperature.py`, `physics/latitude_bands.py`), checked against
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
LOCKED = dict(EARTH, rotation_period_hours=0.1 * 8766.0, period_years=0.1, axial_tilt_deg=0.0)


def _means(out):
    return [out[rt.column(number, "mean_k")] for number in range(1, rt.BAND_COUNT + 1)]


def _swings(out):
    return [out[rt.column(number, "swing_k")] for number in range(1, rt.BAND_COUNT + 1)]


def test_a_body_stores_seven_bands_of_seven_figures_and_a_mean_swing():
    assert len(rt.BAND_COLUMNS) == 49 and len(set(rt.BAND_COLUMNS)) == 49
    assert rt.FIELDS == (*rt.BAND_COLUMNS, "temp_mean_swing_k")
    assert rt.column(1, "mean_k") == "band1_mean_k" and rt.column(7, "state") == "band7_state"
    assert set(rt.compute(EARTH)) == set(rt.FIELDS)


def test_earth_has_poles_35_to_45_k_under_its_mean_and_a_warm_equator():
    out = rt.compute(EARTH)
    means = _means(out)
    assert 288.0 - means[0] == pytest.approx(38.0, abs=8.0)
    assert 8.0 < means[3] - 288.0 < 20.0
    assert means[:3] == pytest.approx(means[:-4:-1])
    assert 0.5 < out[rt.column(4, "swing_k")] < 15.0


def test_the_air_flattens_the_contrast():
    thin = _means(rt.compute({**EARTH, "atmospheric_pressure_pa": 1_000.0}))
    thick = _means(rt.compute({**EARTH, "atmospheric_pressure_pa": 9.2e6}))
    assert (thin[3] - thin[0]) > (thick[3] - thick[0]) > 0.0
    venus = rt.compute(VENUS)
    assert _means(venus)[3] - _means(venus)[0] < 60.0
    assert venus[rt.column(4, "swing_k")] < 100.0


def test_mars_swings_by_tens_of_kelvin():
    out = rt.compute(MARS)
    assert 20.0 < out[rt.column(4, "swing_k")] < 100.0
    assert _means(out)[0] < MARS["surface_temperature_k"] < _means(out)[3]


def test_the_moon_is_locked_with_a_hot_day_side_and_a_cold_night():
    out = rt.compute(MOON)
    means = _means(out)
    assert [out[rt.column(number, "state")] for number in (1, 4, 7)] == ["permanent day side", "terminator band", "permanent night side"]
    assert means[0] > 250.0 and means[6] == rt.NIGHT_FLOOR_K
    assert means == sorted(means, reverse=True)
    assert all(out[rt.column(number, field)] is None for number in range(1, 8) for field in ("day_hours", "night_hours", "longest_day_days", "longest_night_days"))


def test_a_thick_air_column_shares_the_heat_of_a_locked_world():
    means = _means(rt.compute(LOCKED))
    assert means[0] - means[6] < 80.0
    assert means[6] > 200.0


def test_mercury_is_not_locked_and_swings_most_of_its_mean_twice_over():
    out = rt.compute(MERCURY)
    assert out[rt.column(4, "swing_k")] == pytest.approx(1.414 * 440.0, rel=0.5)
    assert out[rt.column(4, "state")] == "day and night"
    means = _means(out)
    assert means[3] > means[2] > means[1] > means[0]


def test_lock_is_a_rotation_period_equal_to_the_orbital_one():
    assert rt.is_locked(8766.0, 1.0) and rt.is_locked(8766.0 * 1.0005, 1.0)
    assert not rt.is_locked(24.0, 1.0) and not rt.is_locked(None, 1.0) and not rt.is_locked(24.0, None)


def test_a_body_without_a_temperature_or_a_spin_has_no_regions():
    assert not any(rt.compute({**EARTH, "surface_temperature_k": None}).values())
    assert not any(rt.compute({**EARTH, "rotation_period_hours": None}).values())
    assert not any(rt.compute(dict(EARTH, rotation_period_hours=0.0)).values())


def test_the_bands_average_back_to_the_mean():
    from planetgen.physics import latitude_bands
    out = rt.compute(EARTH)
    bands = latitude_bands.compute(EARTH, 0.0)
    mean = sum(band.area_share * out[rt.column(number, "mean_k")] for number, band in enumerate(bands, 1))
    assert mean == pytest.approx(288.0, rel=0.002)


def test_the_mean_swing_is_the_area_weighted_mean_of_the_bands_swings():
    from planetgen.physics import latitude_bands
    out = rt.compute(MARS)
    bands = latitude_bands.compute(MARS, 0.0)
    assert out["temp_mean_swing_k"] == pytest.approx(sum(band.area_share * out[rt.column(number, "swing_k")] for number, band in enumerate(bands, 1)))
    assert min(_swings(out)) < out["temp_mean_swing_k"] < max(_swings(out))
    assert rt.compute(MOON)["temp_mean_swing_k"] == 0.0


def test_the_swing_falls_toward_the_pole():
    swings = _swings(rt.compute(EARTH))
    assert swings[3] > swings[2] > swings[1] > swings[0] > 0.0 and swings[:3] == pytest.approx(swings[:-4:-1])


def _world(temperature_c):
    from planetgen.physics import habitability
    return habitability.World(pressure_kpa=101.3, temperature_c=temperature_c, gases={"O2": 0.21})


def test_the_sentence_names_the_equator_the_poles_and_the_swing():
    from planetgen.physics import protection
    text = protection.regional_sentence(_world(15.0), rt.compute(EARTH))
    assert text.startswith("The equator averages 29.7 C and the poles -20.5 C.")
    assert "At the north pole, where it is -20.5 C, an insulated suit is recommended." in text and text.count("recommended") == 1
    assert "Day and night" not in text  # an Earth-like swing of a few kelvin is not worth a sentence


def test_the_sentence_of_a_locked_world_names_its_three_zones():
    from planetgen.physics import protection
    text = protection.regional_sentence(_world(7.0), rt.compute(LOCKED))
    assert "day side" in text and "terminator" in text and "night side" in text


def test_a_body_with_no_regions_has_no_sentence_and_no_summary():
    from planetgen.physics import protection
    assert protection.regional_sentence(_world(15.0), None) == ""
    assert protection.regional_sentence(_world(15.0), {name: None for name in rt.FIELDS}) == ""
    assert rt.summary({}) is None


def test_the_summary_lists_the_seven_bands_and_the_mean_swing():
    from planetgen.physics import latitude_bands
    earth = rt.summary(rt.compute(EARTH))
    assert [band["name"] for band in earth["bands"]] == list(latitude_bands.NAMES)
    assert earth["mean_swing_k"] == pytest.approx(rt.compute(EARTH)["temp_mean_swing_k"])
    assert set(earth["bands"][0]) == {"name", *rt.BAND_FIELDS}
    moon = rt.summary(rt.compute(MOON))
    assert [band["name"] for band in moon["bands"]] == list(latitude_bands.RING_NAMES) and moon["mean_swing_k"] == 0.0
    assert all("day_hours" not in band for band in moon["bands"])


def test_the_daily_range_has_the_stored_value_exactly_in_the_middle():
    """GEN.795, GEN.796: the low and high are the value minus and plus half the peak-to-peak swing, in kelvin."""
    low, high = rt.daily_range(300.0, 9.07)
    assert (low, high) == pytest.approx((295.465, 304.535)) and (low + high) / 2.0 == pytest.approx(300.0)
    assert rt.daily_range(None, 9.0) is None and rt.daily_range(300.0, None) is None and rt.daily_range(300.0, 0.0) is None


def test_the_sentence_gives_the_swing_as_plus_or_minus_and_the_equators_low_and_high():
    from planetgen.physics import protection
    regions = {rt.column(1, "mean_k"): 250.0, rt.column(4, "mean_k"): 300.0, rt.column(7, "mean_k"): 250.0, rt.column(4, "swing_k"): 20.0}
    text = protection.regional_sentence(_world(15.0), regions)
    assert "from the coldest hour to the hottest (±10 K)" in text
    assert "so the equator runs from 16.9 C to 36.9 C over a day." in text
