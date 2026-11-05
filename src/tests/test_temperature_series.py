# tests/test_temperature_series.py

"""
GEN.1263: the temperature of every place at every time (`planetgen.physics.seasonal_swing.temperature_series`, on `surface_cycle`'s cells): the stored mean kept exactly, the worlds of
the research note (a thick-air locked world livable everywhere, a 0.1 bar one about a quarter, Nitrikin-like at 1 and 0.3 bar), the ocean's damping, the ice laid by place, the
filter against a time-stepped integration of the same equation, and the periodic span.
"""

import math

import numpy as np
import pytest

from planetgen.physics import climate_branches, constants, regional_temperature, surface_cycle
from planetgen.physics import seasonal_swing as ss

YEAR_H = 8766.0
EARTH = (24.0, YEAR_H, 23.4, 288.0, 1361.0, 0.30, 101325.0, 9.8, 0.7)
NITRIKIN = (27.9 * 24.0, 28.5 * 24.0, 128.3, 288.0, 1361.0, 0.30, 101325.0, 9.8, 0.0)
SIGMA = constants.STEFAN_BOLTZMANN_CONSTANT


def locked(pressure_pa):
    return (720.0, 720.0, 0.0, 285.0, 1000.0, 0.30, pressure_pa, 9.8, 0.0)


def livable_share(series, low_c=-20.0, high_c=45.0):
    """The share of area and time in the Green band."""
    inside = (series.temperature_k - 273.15 >= low_c) & (series.temperature_k - 273.15 <= high_c)
    return float(np.sum(inside.mean(axis=(1, 2)) * series.grid.row_weight))


def livable_always(series, low_c=-20.0, high_c=45.0):
    inside = (series.temperature_k - 273.15 >= low_c) & (series.temperature_k - 273.15 <= high_c)
    return float(np.sum((inside.mean(axis=2) >= 0.999) * series.grid.row_weight[:, None]) / series.temperature_k.shape[1])


@pytest.mark.parametrize("world", [EARTH, NITRIKIN, locked(101325.0), locked(10132.0), (24.6, 687.0 * 24.0, 25.2, 235.0, 590.0, 0.25, 600.0, 3.7, 0.0),
                                   (12.0, YEAR_H, 30.0, 288.0, 1361.0, 0.30, 30000.0, 9.8, 0.0)])
def test_the_area_and_time_mean_is_the_stored_mean_to_rounding(world):
    series = ss.temperature_series(*world)
    assert series.area_time_mean_k() == pytest.approx(world[3], abs=1.0e-6)
    assert series.mean_k.shape == (12, 12)
    assert series.temperature_k.min() >= ss.FLOOR_K - 1.0e-9 + series.shift_k * 0 - abs(series.shift_k)
    assert series.converged


def test_a_thick_air_locked_world_is_livable_everywhere_and_a_thin_one_about_a_quarter():
    thick, thin, bare = (ss.temperature_series(*locked(p)) for p in (101325.0, 10132.0, 1013.0))
    assert livable_share(thick) == pytest.approx(1.0, abs=0.02)
    assert 0.12 < livable_share(thin) < 0.35
    assert livable_share(bare) < 0.15
    assert thick.share > 0.9 and thin.share < 0.6
    assert thick.share == pytest.approx(regional_temperature.transport_share(1.0, locked=True), abs=1.0e-3)


def test_a_locked_world_is_hot_under_its_star_and_cold_away_from_it():
    series = ss.temperature_series(*locked(10132.0))
    columns = series.mean_k.mean(axis=0)
    hottest = int(np.argmax(columns))
    assert hottest in (0, 11)
    assert columns[hottest] - columns.min() > 100.0
    assert columns.min() == pytest.approx(columns[5:7].min(), abs=1.0e-6) or int(np.argmin(columns)) in (5, 6)


def test_nitrikin_like_at_one_bar_is_livable_and_at_a_third_of_a_bar_is_not_always_so():
    full, thin = ss.temperature_series(*NITRIKIN), ss.temperature_series(*NITRIKIN[:6], 30000.0, *NITRIKIN[7:])
    assert livable_share(full) > 0.9
    assert 0.4 < livable_share(thin) < 0.9
    assert livable_always(thin) < 0.2
    assert full.solar_days == pytest.approx(14.25, abs=0.3)
    assert full.span_days == pytest.approx(28.5)


def test_earth_has_cold_poles_and_a_warm_equator_and_the_ocean_damps_its_cycles():
    wet = ss.temperature_series(*EARTH)
    dry = ss.temperature_series(*EARTH[:8], 0.0)
    equator = slice(5, 7)
    assert wet.mean_k[0:2].mean() < wet.mean_k[equator].mean() - 20.0
    assert wet.mean_k[10:12].mean() < wet.mean_k[equator].mean() - 20.0
    swing_wet = np.ptp(wet.temperature_k[equator], axis=2).mean()
    swing_dry = np.ptp(dry.temperature_k[equator], axis=2).mean()
    assert swing_wet < 0.5 * swing_dry
    assert livable_share(wet) > 0.8


def test_a_world_with_no_tilt_has_no_seasons_and_one_that_turns_slowly_swings_more_than_one_that_turns_fast():
    flat = ss.temperature_series(24.0, YEAR_H, 0.0, 288.0, 1361.0, 0.3, 101325.0, 9.8, 0.0)
    slow = ss.temperature_series(60.0 * 24.0, YEAR_H, 5.0, 288.0, 1361.0, 0.3, 101325.0, 9.8, 0.0)
    fast = ss.temperature_series(12.0, YEAR_H, 5.0, 288.0, 1361.0, 0.3, 101325.0, 9.8, 0.0)
    daily = np.ptp(flat.temperature_k[5], axis=1)
    assert np.all(daily < 30.0)
    assert np.ptp(slow.temperature_k[5:7], axis=2).mean() > 2.0 * np.ptp(fast.temperature_k[5:7], axis=2).mean()


def test_ice_is_laid_by_place_to_the_stored_share_and_brightens_the_cold_cells():
    law = climate_branches.IceAlbedoLaw(ice_free=0.25, ice_covered=0.60)
    iced = ss.temperature_series(*EARTH, ice_fraction=0.2, ice_law=law)
    assert float(np.sum(iced.ice_cover.mean(axis=1) * iced.grid.row_weight)) == pytest.approx(0.2, abs=0.01)
    assert iced.ice_cover[0].mean() > iced.ice_cover[5:7].mean()
    assert iced.albedo[0].mean() > iced.albedo[5:7].mean()
    bare = ss.temperature_series(*EARTH, ice_fraction=0.0, ice_law=law)
    assert bare.ice_cover.max() == 0.0 and np.allclose(bare.albedo, 0.25)
    no_law = ss.temperature_series(*EARTH)
    assert no_law.ice_cover.max() == 0.0 and np.allclose(no_law.albedo, EARTH[5])
    assert iced.mean_k[0].mean() < no_law.mean_k[0].mean()


def _stepped(world, capacity, years, max_steps, emissivity=None):
    """A time-stepped implicit integration of the same equation with one heat capacity: the reference the filter is held to."""
    rotation, orbital, tilt, mean_k, flux, albedo, pressure, gravity, _ = world
    light = surface_cycle.series(rotation, orbital, tilt, None, None, samples_per_day=24, max_steps=max_steps, periodic=True)
    sunlight = np.maximum(light.cosine, 0.0)
    rows, columns, steps = sunlight.shape
    weights = np.repeat(light.grid.row_weight / columns, columns).reshape(rows, columns)
    absorbed_emission = min((1.0 - albedo) * flux / 4.0 / (SIGMA * mean_k ** 4), 1.0) * SIGMA
    emission = absorbed_emission if emissivity is None else emissivity * SIGMA
    share = min(regional_temperature.transport_share(regional_temperature.column_in_earths(pressure, gravity),
                                                     regional_temperature.is_locked(rotation, orbital / (365.25 * 24.0))), ss.MAX_SHARE)
    carried = 4.0 * emission * mean_k ** 3 * share / (1.0 - share)
    forcing = absorbed_emission * mean_k ** 4 / float(np.sum(sunlight.mean(axis=2) * weights)) * sunlight
    dt = light.step_days * 86400.0
    state = np.full((rows, columns), mean_k)
    record = np.empty_like(sunlight)
    for _ in range(years):
        for k in range(steps):
            world_mean = float(np.sum(state * weights))
            guess = state.copy()
            for _newton in range(4):
                residual = capacity * (guess - state) / dt - forcing[:, :, k] + emission * guess ** 4 + carried * (guess - world_mean)
                guess = guess - residual / (capacity / dt + 4.0 * emission * guess ** 3 + carried)
            state = guess
            record[:, :, k] = state
    return record + (mean_k - float(np.sum(record.mean(axis=2) * weights)))


def test_the_filter_agrees_with_a_time_stepped_integration_of_the_same_equation(monkeypatch):
    capacity = 2.0e5
    monkeypatch.setattr(ss, "capacity_by_harmonic", lambda span_s, harmonics, *a, **k: np.full(harmonics, capacity))
    world = (27.9 * 24.0, 28.5 * 24.0, 128.3, 288.0, 1361.0, 0.30, 30000.0, 9.8, 0.0)
    series = ss.temperature_series(*world, max_steps=1500)
    reference = _stepped(world, capacity, 6, 1500, series.emissivity)
    assert series.temperature_k.shape == reference.shape
    assert float(np.sqrt(np.mean((series.temperature_k - reference) ** 2))) < 1.0
    assert float(np.max(np.abs(series.temperature_k - reference))) < 4.0


def test_a_cycle_reaches_a_thin_layer_of_the_sea_by_day_and_a_thick_one_by_year():
    year = YEAR_H * 3600.0
    capacity = ss.capacity_by_harmonic(year, 400, 0.0, 9.8, 1.0)
    assert capacity[1] == pytest.approx(ss.WATER_HEAT_CAPACITY_J_M3_K * ss.OCEAN_MIXED_LAYER_M, rel=0.02)
    day = int(round(365.25))
    assert capacity[day] == pytest.approx(ss.WATER_HEAT_CAPACITY_J_M3_K * math.sqrt(ss.OCEAN_EDDY_DIFFUSIVITY_M2_S * 86400.0 / math.pi), rel=0.02)
    assert capacity[day] < 0.1 * capacity[1]
    ground = ss.capacity_by_harmonic(year, 400, 0.0, 9.8, 0.0)
    assert ground[1] == pytest.approx(regional_temperature.GROUND_INERTIA * math.sqrt(year / math.pi))
    assert ss.capacity_by_harmonic(year, 5, 101325.0, 9.8, 0.0)[1] > ground[1]


def test_the_span_holds_whole_orbits_and_whole_solar_days():
    span, day = surface_cycle.periodic_span(28.5, 14.2)
    assert span == pytest.approx(28.5) and day == pytest.approx(14.25)
    span, day = surface_cycle.periodic_span(365.25, 1.0007)
    assert span == pytest.approx(365.25) and 365.25 / day == pytest.approx(365.0)
    assert surface_cycle.periodic_span(30.0, math.inf) == (30.0, math.inf)
    span, day = surface_cycle.periodic_span(30.0, 30.0 * 1000.0)
    assert span == 30.0 and day == math.inf
    span, day = surface_cycle.periodic_span(10.0, 25.0)
    assert span == pytest.approx(50.0) and day == pytest.approx(25.0)
    span, day = surface_cycle.periodic_span(10.0, 3.7)
    assert span / 10.0 <= surface_cycle.MAX_PERIODIC_ORBITS and abs(span / day - round(span / day)) < 1.0e-9
    series = surface_cycle.series(27.9 * 24.0, 28.5 * 24.0, 128.3, periodic=True, max_steps=3000)
    assert series.cosine.shape[2] * series.step_days == pytest.approx(28.5)
    first, last = series.cosine[:, :, 0], series.cosine[:, :, -1]
    assert np.max(np.abs(first - last)) < 0.1


def test_the_series_is_deterministic_and_refuses_what_it_cannot_use():
    first, second = ss.temperature_series(*EARTH), ss.temperature_series(*EARTH)
    assert np.array_equal(first.temperature_k, second.temperature_k)
    with pytest.raises(ss.SeasonalSwingError):
        ss.temperature_series(24.0, YEAR_H, 23.4, 0.0, 1361.0, 0.3, 101325.0, 9.8, 0.7)
    with pytest.raises(ss.SeasonalSwingError):
        ss.temperature_series(24.0, YEAR_H, 23.4, 288.0, 1361.0, 0.3, 101325.0, 9.8, 1.5)
    with pytest.raises(ss.SeasonalSwingError):
        ss.temperature_series(0.0, YEAR_H, 23.4, 288.0, 1361.0, 0.3, 101325.0, 9.8, 0.5)


def test_a_body_with_no_air_has_no_transport_and_the_floor_holds():
    series = ss.temperature_series(24.6 * 24.0, 24.6 * 24.0 * 3.0, 0.0, 250.0, 1361.0, 0.1, 0.0, 1.6, 0.0)
    assert series.share == 0.0
    assert series.temperature_k.min() >= ss.FLOOR_K - abs(series.shift_k) - 1.0e-9
