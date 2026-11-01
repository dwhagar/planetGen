# tests/test_outside_hours.py

"""GEN.768: the hours a person can be outside in a day, the smallest of the daylight, the dose allowance, the heat and cold
a worn garment copes with, the breathing gas and the ultraviolet, and the limit that binds."""

import itertools

import pytest

from planetgen.physics import constants, equipment, habitability, outside_hours
from planetgen.physics.outside_hours import Temperatures
from planetgen.physics.protection_limits import (
    GARMENTS,
    LIMIT_ROWS,
    limit,
    supply_hours,
)

AIR = {"N2": 0.78, "O2": 0.21, "Ar": 0.01}
C = constants.CELSIUS_ZERO_K


def _world(**fields):
    base = {"pressure_kpa": 101.3, "temperature_c": 15.0, "gases": AIR, "relative_humidity": 0.5, "uv_index": 1.0}
    base.update(fields)
    return habitability.World(**base)


def _day(world, daylight=12.0, days=365.25, swing=0.0):
    """`(hours, limit)` for `world` with its listed equipment and a day and night `swing` kelvin wide."""
    cold, hot = habitability.temperature_extremes_c(world)
    return outside_hours.outside_hours_day(world, equipment.compute(world), daylight, days,
                                           Temperatures(cold + C - swing / 2, hot + C + swing / 2))


def test_earth_gets_the_full_daylight():
    hours, binding = _day(_world(), daylight=14.0, swing=10.0)
    assert (hours, binding) == (14.0, "daylight")


def test_a_day_without_light_is_none_and_the_light_binds():
    assert _day(_world(), daylight=0.0) == (0.0, "daylight")


def test_a_hot_desert_world_is_cut_by_the_heat_near_the_top_of_its_garment():
    """120 C (a liquid-cooled suit is listed from 45 C): with a 40 K swing the afternoon reaches 140 C."""
    dry = _world(temperature_c=120.0, relative_humidity=0.1)
    hours, binding = _day(dry, swing=40.0)
    assert binding == "heat" and 1.0 < hours < 12.0
    assert _day(dry, swing=0.0) == (12.0, "daylight")


def test_the_hotter_the_afternoon_the_fewer_the_hours_in_the_same_suit():
    garment = GARMENTS[(3, "heat")]
    hours = [outside_hours.exposure_hours(garment, garment, 101.3, Temperatures(250.0, 273.15 + t), 0.1)[0] for t in (300, 350, 390, 410)]
    assert all(later < earlier for earlier, later in itertools.pairwise(hours))


def test_a_cold_world_is_cut_by_the_cold_near_the_bottom_of_its_garment():
    """-49 C: an insulated suit is listed down to -50 C and a person working in it loses more than the body makes."""
    hours, binding = _day(_world(temperature_c=-49.0), swing=0.0)
    assert binding == "cold" and 1.0 < hours < 12.0


def test_a_heated_suit_keeps_the_core_warm_so_air_binds_in_the_deep_cold():
    assert _day(_world(temperature_c=-150.0)) == (supply_hours(2), "air")


def test_a_vacuum_world_is_limited_by_the_gas_in_the_suit():
    hours, binding = _day(_world(pressure_kpa=0.0, temperature_c=0.0, gases={}), swing=40.0)
    assert (hours, binding) == (supply_hours(2), "air") and hours < 12.0


def test_a_high_dose_world_has_its_yearly_hours_spread_over_the_work_days():
    world = _world(surface_dose_sv_yr=5.0, dose_gcr_sv_yr=4.0, dose_sep_sv_yr=0.0, dose_ground_sv_yr=1.0)
    levels = equipment.compute(world)
    assert levels["outside_hours"] is not None
    hours, binding = _day(world, days=300.0)
    assert binding == "dose" and hours == pytest.approx(levels["outside_hours"] / 300.0)


def test_a_dose_that_is_not_a_limit_does_not_bind():
    world = _world(surface_dose_sv_yr=0.5, dose_gcr_sv_yr=0.4, dose_sep_sv_yr=0.0, dose_ground_sv_yr=0.1)
    assert equipment.compute(world)["outside_hours"] is None
    assert _day(world)[1] == "daylight"


def test_a_light_that_is_past_the_garments_rating_cuts_the_hours_in_proportion():
    """A full-cover suit is rated to the lethal edge (300 times Earth's): at 500 the day is 300/500 of the daylight."""
    hours, binding = _day(_world(uv_index=500.0), daylight=10.0)
    assert binding == "uv" and hours == pytest.approx(10.0 * limit("uv_lethal_index") / 500.0)


def test_a_sealed_suit_is_rated_for_any_light():
    assert outside_hours.uv_edge(3) == float("inf") and outside_hours.uv_edge(1) == limit("uv_suit_index")


def test_a_body_that_needs_a_vehicle_has_no_hours():
    world = _world(temperature_c=600.0, pressure_kpa=101.3)
    assert equipment.compute(world)["eq_body"] == 6
    assert _day(world, swing=0.0) == (0.0, "vehicle")


def test_the_tie_goes_to_the_first_limit_in_the_order():
    assert outside_hours.LIMITS[:2] == ("daylight", "dose")


def test_a_chemical_barrier_suit_that_cannot_dump_heat_is_cut_by_the_heat_even_at_earths_temperature():
    """A sealed, uncooled suit (Level A) holds a working person's heat: about an hour before the core limit."""
    garment = GARMENTS[(3, "chemical")]
    heat, cold = outside_hours.exposure_hours(garment, garment, 101.3, Temperatures(288.0, 288.0), 0.5)
    assert cold is None and 0.5 < heat < 3.0


def test_sweat_has_less_to_carry_in_thin_air_and_none_in_a_sealed_suit_or_wet_air():
    open_suit, sealed = GARMENTS[(2, "heat")], GARMENTS[(3, "heat")]
    thick = outside_hours._sweat_w(open_suit, 101.3, 330.0, 0.1)
    thin = outside_hours._sweat_w(open_suit, 5.0, 330.0, 0.1)
    assert 0.0 < thin < thick <= limit("sweat_max_w")
    assert outside_hours._sweat_w(sealed, 101.3, 330.0, 0.1) == 0.0
    assert outside_hours._sweat_w(open_suit, 101.3, 330.0, 1.0) == 0.0


def test_the_saturation_pressure_of_water_matches_its_known_points():
    assert outside_hours.saturation_kpa(373.15) == pytest.approx(101.3, rel=0.05)
    assert outside_hours.saturation_kpa(293.15) == pytest.approx(2.34, rel=0.02)
    assert outside_hours.saturation_kpa(263.15) == pytest.approx(0.26, rel=0.05)


def test_the_body_at_a_temperature_is_the_ladders_own_choice():
    assert equipment.body_at(_world(), 15.0) == (0, "none")
    assert equipment.body_at(_world(), -49.0) == (2, "cold")
    assert equipment.body_at(_world(), 120.0) == (2, "heat")
    assert equipment.body_at(_world(), -90.0)[0] == 3


def test_every_garment_row_is_a_listed_body_name_and_has_a_power_only_where_it_is_active():
    for (level, kind), garment in GARMENTS.items():
        assert (level, kind) in equipment.BODY_NAMES or (level, kind) == (0, "none")
        assert garment.insulation_clo > 0 and garment.cooling_w >= 0 and garment.heating_w >= 0
        if not garment.sealed:
            assert garment.cooling_w <= 600.0
    for level in (3, 4, 5):
        assert all(garment.sealed for (lv, _), garment in GARMENTS.items() if lv == level)


def test_the_supply_hours_follow_the_shift_and_ambient_air_has_none():
    assert supply_hours(0) is None
    assert supply_hours(1) == supply_hours(4) == limit("shift_hours")
    assert 0 < supply_hours(2) < supply_hours(3) < supply_hours(4)


def test_the_work_days_of_a_year_are_the_solar_days():
    assert outside_hours.work_days_per_year(23.9345, 1.0) == pytest.approx(365.2572, abs=0.01)
    assert outside_hours.work_days_per_year(24.0, 1.0, locked=True) == pytest.approx(365.25)
    assert outside_hours.work_days_per_year(None, 2.0) == pytest.approx(730.5)
    assert outside_hours.work_days_per_year(10.0, 0.0001) == pytest.approx(0.0001 * 8766.0 / 24.0)


def test_the_new_limit_rows_have_a_value_a_unit_and_a_meaning():
    rows = {row.name: row for row in LIMIT_ROWS}
    for name in ("body_mass_kg", "body_area_m2", "core_hot_limit_c", "core_cold_limit_c", "work_metabolism_w", "shift_hours",
                 "vent_factor", "sweat_max_w"):
        assert rows[name].value > 0 and rows[name].unit and rows[name].meaning
    assert limit("core_cold_limit_c") < limit("core_normal_c") < limit("core_hot_limit_c")


# --- GEN.768 part 2: the hours of a body from its stored columns ---------------------------------------------------------------

def _body(world, **columns):
    """`for_body` of `world` with the columns of a body that turns in 24 h over a year, and `columns` over them."""
    values = {"rotation_period_hours": 24.0, "period_years": 1.0, "daylight_hours": 12.0, "band4_swing_k": 0.0}
    values.update(columns)
    return outside_hours.for_body(world, equipment.compute(world), values)


def test_earth_from_its_columns_has_the_full_daylight():
    assert _body(_world(), daylight_hours=12.1, band4_swing_k=10.0) == (12.1, "daylight")


def test_a_hot_desert_from_its_columns_is_cut_by_the_heat():
    hours, binding = _body(_world(temperature_c=120.0, relative_humidity=0.1), band4_swing_k=40.0)
    assert binding == "heat" and 1.0 < hours < 12.0


def test_a_cold_world_from_its_columns_is_cut_by_the_cold():
    hours, binding = _body(_world(temperature_c=-49.0))
    assert binding == "cold" and 1.0 < hours < 12.0


def test_a_high_dose_world_from_its_columns_has_the_yearly_hours_over_its_solar_days():
    world = _world(surface_dose_sv_yr=5.0, dose_gcr_sv_yr=4.0, dose_sep_sv_yr=0.0, dose_ground_sv_yr=1.0)
    yearly = equipment.compute(world)["outside_hours"]
    hours, binding = _body(world)
    assert binding == "dose" and hours == pytest.approx(yearly / outside_hours.work_days_per_year(24.0, 1.0))


def test_a_vacuum_world_from_its_columns_is_limited_by_the_gas_in_the_suit():
    assert _body(_world(pressure_kpa=0.0, temperature_c=0.0, gases={}), band4_swing_k=40.0) == (supply_hours(2), "air")


def test_a_locked_world_works_the_terminator_belt_for_a_whole_shift_at_the_belts_temperature():
    """A world that turns once an orbit has no stored daylight and no day of its own: the belt's 24 hours of light, and the belt's temperature."""
    belt = 288.15
    assert _body(_world(), rotation_period_hours=8766.0, period_years=1.0, daylight_hours=None, band4_mean_k=belt) == \
        (outside_hours.LOCKED_BELT_LIGHT_HOURS, "daylight")
    hours, binding = _body(_world(temperature_c=-49.0), rotation_period_hours=8766.0, period_years=1.0, daylight_hours=None,
                           band4_mean_k=224.15)
    assert binding == "cold" and hours < outside_hours.LOCKED_BELT_LIGHT_HOURS


def test_a_locked_world_without_a_terminator_temperature_has_no_answer():
    assert _body(_world(), rotation_period_hours=8766.0, period_years=1.0, daylight_hours=None) == (None, None)


def test_a_body_without_stored_daylight_or_rotation_has_no_answer_and_none_is_not_zero():
    assert _body(_world(), daylight_hours=None) == (None, None)
    assert _body(_world(), rotation_period_hours=None) == (None, None)
    assert _body(_world(), daylight_hours=0.0) == (0.0, "daylight")


def test_a_body_that_needs_a_vehicle_has_no_hours_from_its_columns():
    world = _world(pressure_kpa=0.0, temperature_c=900.0, gases={})
    assert equipment.compute(world)["eq_body"] >= 6, "the test world must be one that needs a vehicle"
    assert _body(world) == (0.0, "vehicle")
