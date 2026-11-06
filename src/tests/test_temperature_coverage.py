# tests/test_temperature_coverage.py

"""
GEN.1264: what is read off the temperature of every place (`planetgen.physics.temperature_coverage`): the shares of area and time in each PHI band, the livable share and its two
pieces, the extremes by place and season, the best site and its longest spell indoors, the wet bulb rule per sample, and the ten example worlds of the research note
(`reports/temperature-coverage-2026-11-04.md`, section 3).
"""

import math
from types import SimpleNamespace

import numpy as np
import pytest

from planetgen.physics import habitability, latitude_bands, surface_cycle
from planetgen.physics import seasonal_swing as ss
from planetgen.physics import temperature_coverage as tc

YEAR_H = 8766.0
WORLDS = {
    "earth": (24.0, YEAR_H, 23.4, 288.0, 1361.0, 0.30, 101325.0, 9.8, 0.7),
    "earth tilt 0": (24.0, YEAR_H, 0.0, 288.0, 1361.0, 0.30, 101325.0, 9.8, 0.7),
    "spin 12 h thin": (12.0, YEAR_H, 30.0, 288.0, 1361.0, 0.30, 30000.0, 9.8, 0.0),
    "locked 1 bar": (720.0, 720.0, 0.0, 285.0, 1000.0, 0.30, 101325.0, 9.8, 0.0),
    "locked 0.1 bar": (720.0, 720.0, 0.0, 285.0, 1000.0, 0.30, 10132.0, 9.8, 0.0),
    "locked 0.01 bar": (720.0, 720.0, 0.0, 285.0, 1000.0, 0.30, 1013.0, 9.8, 0.0),
    "nitrikin 1 bar": (27.9 * 24.0, 28.5 * 24.0, 128.3, 288.0, 1361.0, 0.30, 101325.0, 9.8, 0.0),
    "nitrikin 0.3 bar": (27.9 * 24.0, 28.5 * 24.0, 128.3, 288.0, 1361.0, 0.30, 30000.0, 9.8, 0.0),
    "slow 60 d spin": (60.0 * 24.0, YEAR_H, 5.0, 288.0, 1361.0, 0.30, 101325.0, 9.8, 0.0),
    "tilt 98": (24.0, YEAR_H, 98.0, 288.0, 1361.0, 0.30, 101325.0, 9.8, 0.7),
}


@pytest.fixture(scope="module")
def reads():
    return {name: tc.coverage(ss.temperature_series(*world)) for name, world in WORLDS.items()}


def synthetic(values, steps=None):
    """A series for a 12 by 12 grid: `values` is `[row, column, step]` kelvin (or a function of the grid)."""
    grid = surface_cycle.cell_grid()
    temperature = np.asarray(values, dtype=float)
    count = temperature.shape[2]
    return SimpleNamespace(grid=grid, temperature_k=temperature, time_days=(np.arange(count) + 0.5), step_days=1.0, orbital_days=float(count), solar_days=1.0)


def uniform(celsius, steps=8):
    return np.full((12, 12, steps), celsius + 273.15)


def test_the_bands_are_the_phi_temperature_bands_and_the_shares_add_to_one():
    for celsius, band in ((15.0, "blue"), (-10.0, "green"), (40.0, "green"), (-40.0, "yellow"), (100.0, "yellow"), (-100.0, "red"), (150.0, "red")):
        read = tc.coverage(synthetic(uniform(celsius)))
        assert getattr(read, band) == pytest.approx(1.0), (celsius, band)
        assert read.blue + read.green + read.yellow + read.red == pytest.approx(1.0)
        assert read.livable_share == pytest.approx(1.0 if band in ("blue", "green") else 0.0)
    cold, hot = tc.coverage(synthetic(uniform(-30.0))), tc.coverage(synthetic(uniform(60.0)))
    assert cold.cold == pytest.approx(1.0) and cold.hot == 0.0
    assert hot.hot == pytest.approx(1.0) and hot.cold == 0.0


def test_the_edges_of_the_bands_are_the_habitability_edges():
    edges = habitability.TEMPERATURE_BAND
    for celsius, band in ((edges[3], "blue"), (edges[4], "blue"), (edges[2], "green"), (edges[5], "green"), (edges[1], "yellow"), (edges[6], "yellow")):
        assert getattr(tc.coverage(synthetic(uniform(celsius))), band) == pytest.approx(1.0)


def test_shares_weigh_area_and_time_alike():
    values = uniform(15.0, 10)
    values[:, :, :3] = 60.0 + 273.15
    read = tc.coverage(synthetic(values))
    assert read.blue == pytest.approx(0.7) and read.hot == pytest.approx(0.3)
    values = uniform(15.0)
    values[:2] = 60.0 + 273.15
    weight = 2 * surface_cycle.cell_grid().row_weight[0]
    assert tc.coverage(synthetic(values)).hot == pytest.approx(weight)


def test_livable_time_per_place_and_the_areas_it_gives():
    values = uniform(15.0, 10)
    values[:3, :, :4] = 60.0 + 273.15          # three rows livable 60 percent of the time
    values[3:5, :, :] = 60.0 + 273.15          # two rows never livable
    values[5:7, :, :1] = 60.0 + 273.15         # two rows livable 90 percent of the time
    read = tc.coverage(synthetic(values))
    twelfth = 1.0 / 12.0
    assert read.area_always == pytest.approx(5 * twelfth)
    assert read.area_most == pytest.approx(7 * twelfth)
    assert read.livable_area == pytest.approx(10 * twelfth)
    assert read.area_ever == pytest.approx(10 * twelfth)
    expected = (5 * 1.0 + 2 * 0.9 + 3 * 0.6) * twelfth
    assert read.livable_share == pytest.approx(expected)
    assert read.livable_time == pytest.approx(expected / (10 * twelfth))
    assert read.place_time[0, 0] == pytest.approx(0.6) and read.place_time[3, 0] == 0.0


def test_the_wet_bulb_rule_is_applied_to_each_hot_sample():
    ceiling = tc.wet_bulb_ceiling_c(0.8)
    assert ceiling < 45.0 and habitability.wet_bulb_c(ceiling, 0.8) == pytest.approx(habitability.WET_BULB_CRIT_C, abs=1.0e-6)
    assert tc.wet_bulb_ceiling_c(None) == math.inf
    assert tc.wet_bulb_ceiling_c(0.05) == math.inf or tc.wet_bulb_ceiling_c(0.05) > tc.wet_bulb_ceiling_c(0.8)
    values = uniform(15.0, 10)
    values[:, :, :4] = ceiling + 1.0 + 273.15
    plain, humid = tc.coverage(synthetic(values)), tc.coverage(synthetic(values), 0.8)
    assert plain.livable_share == pytest.approx(1.0) and plain.humid == 0.0
    assert humid.livable_share == pytest.approx(0.6) and humid.humid == pytest.approx(0.4) and humid.yellow == pytest.approx(0.4)
    assert humid.ceiling_c == pytest.approx(ceiling)
    dry = tc.coverage(synthetic(values), 0.05)
    assert dry.humid <= humid.humid
    hotter = uniform(55.0)
    assert tc.coverage(synthetic(hotter), 0.8).humid == 0.0


def test_the_extremes_are_the_samples_and_say_where_and_when():
    values = uniform(15.0, 8)
    values[2, 3, 5] = 250.0
    values[9, 7, 1] = 330.0
    read = tc.coverage(synthetic(values))
    grid = surface_cycle.cell_grid()
    assert (read.coldest.temperature_k, read.coldest.row, read.coldest.column, read.coldest.time_days) == (250.0, 2, 3, 5.5)
    assert read.coldest.latitude_deg == grid.latitude_deg[2] and read.coldest.longitude_deg == grid.longitude_deg[3]
    assert (read.hottest.temperature_k, read.hottest.row, read.hottest.column, read.hottest.time_days) == (330.0, 9, 7, 1.5)
    assert read.by_column[3][0] == 250.0 and read.by_column[7][1] == 330.0 and read.by_column[0] == (288.15, 288.15)
    bands = [tc.band_of_latitude(lat) for lat in grid.latitude_deg]
    assert read.by_band[bands[2]][0] == 250.0 and read.by_band[bands[9]][1] == 330.0
    assert len(read.by_band) == latitude_bands.BAND_COUNT
    assert read.by_quarter[2][0] == 250.0 and read.by_quarter[0][1] == 330.0
    assert read.sectors is None


def test_the_bands_of_latitude_are_the_seven_of_the_stored_bands():
    assert [tc.band_of_latitude(lat) for lat in (90.0, 70.0, 50.0, 20.0, 0.0, -20.0, -50.0, -70.0, -90.0)] == [0, 0, 1, 2, 3, 4, 5, 6, 6]


def test_the_longest_spell_out_of_the_band_joins_the_end_of_the_series_to_its_start():
    assert tc.longest_out_days([1, 1, 1, 1], 2.0) == 0.0
    assert tc.longest_out_days([0, 0, 0], 2.0) == 6.0
    assert tc.longest_out_days([1, 0, 0, 1, 0, 1, 1, 1], 0.5) == 1.0
    assert tc.longest_out_days([0, 1, 1, 1, 1, 0, 0, 0], 1.0) == 4.0
    assert tc.longest_out_days([0, 0, 1, 0, 0, 0, 1, 0], 1.0) == 3.0


def test_the_best_site_is_the_most_livable_cell_and_a_tie_goes_to_the_more_moderate():
    values = uniform(60.0, 10)
    values[4, 4, :] = 15.0 + 273.15                 # always livable, steady
    values[6, 6, :] = 15.0 + 273.15
    values[6, 6, 0] = 40.0 + 273.15                 # always livable, but swings more
    values[2, 2, :] = 15.0 + 273.15
    values[2, 2, :2] = 70.0 + 273.15                # livable 80 percent of the time
    site = tc.coverage(synthetic(values)).best_site
    assert (site.row, site.column) == (4, 4) and site.livable_time == 1.0 and site.longest_out_days == 0.0
    assert site.coldest_k == site.hottest_k == pytest.approx(288.15)
    values[4, 4, :] = 60.0 + 273.15
    site = tc.coverage(synthetic(values)).best_site
    assert (site.row, site.column) == (6, 6)
    values[6, 6, :] = 60.0 + 273.15
    site = tc.coverage(synthetic(values)).best_site
    assert (site.row, site.column) == (2, 2) and site.livable_time == pytest.approx(0.8)
    assert site.longest_out_days == 2.0 and site.hottest_k == pytest.approx(343.15)


def test_a_world_with_no_livable_cell_still_names_its_best_cell_with_no_livable_time():
    read = tc.coverage(synthetic(uniform(150.0)))
    assert read.livable_share == 0.0 and read.best_site.livable_time == 0.0
    assert read.best_site.longest_out_days == 8.0
    assert read.area_ever == 0.0 and read.livable_time == 0.0


def test_a_locked_world_reports_its_day_side_terminator_and_night_side():
    series = ss.temperature_series(*WORLDS["locked 0.1 bar"])
    read = tc.coverage(series)
    assert set(read.sectors) == set(tc.SECTORS)
    day, terminator, night = (read.sectors[name] for name in tc.SECTORS)
    assert day[2] > terminator[2] > night[2] and day[1] > terminator[1] >= night[1]    # the dark cells all share the air's one temperature
    assert terminator[0] > day[0] and terminator[0] > night[0]
    masks = tc.sector_masks(series)
    areas = [float(np.sum(np.repeat(series.grid.row_weight / 12, 12).reshape(12, 12)[mask])) for mask in masks.values()]
    assert areas == pytest.approx([0.25, 0.5, 0.25], abs=0.05)
    assert sum(areas) == pytest.approx(1.0)
    assert tc.coverage(ss.temperature_series(*WORLDS["earth"])).sectors is None


def test_the_ten_worlds_of_the_note_read_as_the_note_says(reads):
    livable = {name: read.livable_share for name, read in reads.items()}
    assert livable["locked 1 bar"] == pytest.approx(1.0, abs=0.02) and reads["locked 1 bar"].area_always > 0.98
    assert 0.12 < livable["locked 0.1 bar"] < 0.35
    assert livable["locked 0.01 bar"] < livable["locked 0.1 bar"] < livable["locked 1 bar"]
    assert livable["nitrikin 1 bar"] > 0.9 > livable["nitrikin 0.3 bar"] > 0.4
    assert reads["nitrikin 1 bar"].area_always > 0.9 and reads["nitrikin 0.3 bar"].area_always < 0.2
    assert livable["earth"] > 0.8 and reads["earth"].blue > 0.4
    assert livable["earth tilt 0"] > 0.8 and reads["earth tilt 0"].cold > 0.1
    assert 0.5 < livable["spin 12 h thin"] < 0.8 and reads["spin 12 h thin"].hot > 0.1 and reads["spin 12 h thin"].cold > 0.1
    assert 0.5 < livable["slow 60 d spin"] < 0.8 and reads["slow 60 d spin"].area_always < 0.1
    assert livable["spin 12 h thin"] < livable["earth"] and livable["slow 60 d spin"] < livable["earth"]
    for read in reads.values():
        assert read.blue + read.green + read.yellow + read.red == pytest.approx(1.0)
        assert read.livable_share == pytest.approx(read.blue + read.green)
        assert read.area_always <= read.area_most <= read.livable_area <= read.area_ever <= 1.0 + 1e-9
        assert read.livable_area * read.livable_time <= read.livable_share + 1e-9
        assert read.coldest.temperature_k <= min(low for low, _ in read.by_band) + 1e-9
        assert read.hottest.temperature_k >= max(high for _, high in read.by_band) - 1e-9
        assert read.best_site.livable_time == pytest.approx(read.place_time.max())


def test_a_thin_air_locked_world_is_livable_only_near_the_terminator():
    read = tc.coverage(ss.temperature_series(*WORLDS["locked 0.1 bar"]))
    assert read.sectors["terminator belt"][0] > 0.35
    assert read.sectors["night side"][0] < 0.05 and read.sectors["day side"][0] < 0.05
    site = read.best_site
    assert abs(site.longitude_deg - 90.0) < 40.0 or abs(site.longitude_deg - 270.0) < 40.0


def test_the_best_site_of_a_slow_spinner_has_to_go_indoors_for_days_and_that_of_a_mild_world_does_not(reads):
    slow = reads["slow 60 d spin"]
    assert 0.7 < slow.best_site.livable_time < 1.0
    assert 5.0 < slow.best_site.longest_out_days < 40.0
    for name in ("earth", "nitrikin 1 bar", "nitrikin 0.3 bar", "locked 1 bar"):
        assert reads[name].best_site.livable_time == 1.0 and reads[name].best_site.longest_out_days == 0.0
