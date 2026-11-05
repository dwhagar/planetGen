# tests/test_surface_cycle.py

"""
GEN.1262: the sunlight on every cell of a body through its year and its turns (`physics/surface_cycle.py`). Earth, a locked world, a backwards spinner and a 28-day spinner near
resonance are test targets, never inputs, and every one of them goes through the same functions.
"""

import math

import numpy as np
import pytest

from planetgen.physics import surface_cycle as sc


def test_the_grid_has_rows_of_equal_area_and_columns_of_equal_longitude():
    grid = sc.cell_grid()
    assert len(grid.latitude_deg) == sc.CELL_ROWS and len(grid.longitude_deg) == sc.CELL_COLUMNS
    assert grid.row_weight.sum() == pytest.approx(1.0)
    edges = np.sin(np.radians(grid.latitude_deg))
    assert np.allclose(np.diff(edges), 2.0 / sc.CELL_ROWS)       # equal steps in the sine of the latitude are equal areas
    assert grid.latitude_deg[0] == pytest.approx(-grid.latitude_deg[-1]) and np.allclose(np.diff(grid.longitude_deg), 360.0 / sc.CELL_COLUMNS)


def test_the_solar_day_is_one_formula_for_prograde_retrograde_and_locked_bodies():
    assert sc.solar_day_days(23.934, 8766.0, 23.44) == pytest.approx(1.0, abs=0.001)
    assert sc.solar_day_days(27.9 * 24.0, 28.5 * 24.0, 128.3) == pytest.approx(14.1, abs=0.05)       # backwards: the day is shorter than either period
    assert sc.solar_day_days(27.9 * 24.0, 28.5 * 24.0, 20.0) > 1000.0       # forwards at the same rates: the star hardly moves
    assert sc.solar_day_days(720.0, 720.0, 0.0) == math.inf       # once an orbit: no day
    assert sc.solar_day_days(5832.5, 5400.0, 177.4) == pytest.approx(1.0 / (1.0 / 243.02 + 1.0 / 225.0), rel=0.01)       # Venus turns backwards


def test_earth_gets_a_quarter_of_the_light_on_average_and_its_noon_sun_follows_the_declination():
    series = sc.series(23.934, 8766.0, 23.44)
    assert series.step_days <= 1.0 / sc.SAMPLES_PER_SOLAR_DAY + 1.0e-9
    lit = np.maximum(series.cosine, 0.0)
    assert float(np.mean(lit * series.grid.row_weight[:, None, None])) * sc.CELL_ROWS == pytest.approx(0.25, abs=0.005)       # a sphere takes a quarter of the beam on average
    june = series.declination_rad.argmax()
    north = int(np.argmax(series.grid.latitude_deg > 60.0))
    assert float(series.cosine[north, :, june - 24:june + 24].max()) == pytest.approx(math.cos(math.radians(series.grid.latitude_deg[north] - 23.44)), abs=0.01)       # noon sun
    assert float(series.declination_rad.max()) == pytest.approx(math.radians(23.44), abs=0.01)


def test_a_tilt_past_90_degrees_has_the_seasons_of_its_supplement_and_a_backwards_day():
    forward = sc.series(23.934, 8766.0, 30.0)
    backward = sc.series(23.934, 8766.0, 150.0)
    assert float(backward.declination_rad.max()) == pytest.approx(float(forward.declination_rad.max()))
    assert backward.solar_days == pytest.approx(1.0 / (1.0 / 0.99726 + 1.0 / 365.25), rel=1e-3)
    assert forward.solar_days == pytest.approx(1.0 / (1.0 / 0.99726 - 1.0 / 365.25), rel=1e-3)


def test_a_locked_body_has_a_star_that_stands_still_over_every_cell():
    series = sc.series(720.0, 720.0, 0.0)
    assert series.solar_days == math.inf
    assert np.allclose(series.cosine, series.cosine[:, :, :1])       # nothing changes with time
    substellar = float(series.cosine[:, :, 0].max())
    assert substellar > 0.9
    dayside = int((series.cosine[:, :, 0] > 0.0).sum())
    assert dayside == pytest.approx(0.5 * sc.CELL_ROWS * sc.CELL_COLUMNS, abs=sc.CELL_ROWS)       # half the cells face the star
    assert series.cosine.shape[2] >= 2 and series.step_days == pytest.approx(30.0 / series.cosine.shape[2])


def test_a_28_day_spinner_near_resonance_is_stepped_through_a_solar_day_of_14_days_and_a_year_of_28_5():
    series = sc.series(27.9 * 24.0, 28.5 * 24.0, 128.3)
    assert series.solar_days == pytest.approx(14.1, abs=0.05) and series.orbital_days == pytest.approx(28.5)
    assert series.step_days <= 14.1 / sc.SAMPLES_PER_SOLAR_DAY + 1.0e-9
    assert series.time_days[-1] < 28.5 and series.time_days[-1] > 28.0
    pole = int(np.argmax(series.grid.latitude_deg))
    column = series.cosine[pole, 0, :]
    assert column.max() > 0.5 and column.min() < -0.5       # the north is lit and dark in turn


def test_the_step_is_the_least_of_the_solar_day_the_thermal_time_and_the_year_within_the_cap():
    step, count = sc.time_step_days(365.25, 1.0)
    assert count == int(math.ceil(365.25 * 48.0)) and step == pytest.approx(1.0 / 48.0, rel=1e-3)       # a 48th of a solar day, 17,532 steps
    step, count = sc.time_step_days(365.25 * 84.0, 0.72)
    assert count == sc.MAX_STEPS and step == pytest.approx(365.25 * 84.0 / sc.MAX_STEPS)       # a year of 42,000 solar days is cut at the cap
    step, count = sc.time_step_days(28.5, 14.1)
    assert step <= 14.1 / 48.0 and count == sc.STEPS_PER_YEAR and step == pytest.approx(28.5 / sc.STEPS_PER_YEAR)       # the year over 3000 is the least of the three limits
    step, _ = sc.time_step_days(28.5, 14.1, thermal_days=1.0)
    assert step <= 0.25 + 1e-9
    step, count = sc.time_step_days(10.0, math.inf)
    assert step == pytest.approx(10.0 / count) and step <= 10.0 / sc.STEPS_PER_YEAR + 1e-12


def test_a_ring_shade_takes_light_off_the_lit_cells_only():
    plain = sc.series(23.934, 8766.0, 23.44, max_steps=4000)
    shaded = sc.series(23.934, 8766.0, 23.44, transmission=lambda latitudes, declinations: np.full((len(latitudes), len(declinations)), 0.5), max_steps=4000)
    lit = plain.cosine > 0.0
    assert np.allclose(shaded.cosine[lit], 0.5 * plain.cosine[lit]) and np.allclose(shaded.cosine[~lit], plain.cosine[~lit])
