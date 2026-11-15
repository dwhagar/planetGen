"""GEN.553 birth step, for GEN.1611: the luminosity of a star while its disk was there (`physics/birth_luminosity.py`)."""

import math

import pytest

from planetgen.datasets import baraffe_isochrones, fetch
from planetgen.physics import birth_luminosity as bl, stellar_evolution
from planetgen.util import draw

DISK_AGE_GY = 0.01      # GEN.1611's FORMATION_DISK_AGE_GY, the age the callers pass


@pytest.fixture(scope="module")
def grid():
    try:
        return baraffe_isochrones.isochrones()
    except fetch.DataMissing:
        pytest.skip("the BHAC15 isochrones are not downloaded")


@pytest.mark.parametrize("mass, low, high", [
    (0.3, 0.05, 0.065),      # BHAC15 at 10 Myr: 0.0589 (the main sequence holds 0.02 Sun)
    (1.0, 0.40, 0.50),       # 0.457 at 10 Myr; the Sun is 0.7 of its present light on arrival at the main sequence
])
def test_a_contracting_star_is_read_from_the_bhac15_row_not_from_the_main_sequence(grid, mass, low, high):
    value = bl.birth_luminosity_sol(mass, 0.0, DISK_AGE_GY)
    assert low < value < high
    assert value == pytest.approx(10.0 ** baraffe_isochrones.interpolate(mass, DISK_AGE_GY, grid)["log_luminosity_solar"])


def test_an_m_dwarf_is_several_times_brighter_at_the_disk_age_than_on_the_main_sequence(grid):
    on_main_sequence = bl.birth_luminosity_sol(0.3, 0.0, 5.0)
    assert bl.birth_luminosity_sol(0.3, 0.0, DISK_AGE_GY) > 2.5 * on_main_sequence


def test_a_star_of_four_solar_masses_is_near_its_zero_age_main_sequence_at_the_disk_age():
    value = bl.birth_luminosity_sol(4.0, 0.0, DISK_AGE_GY)
    zams = stellar_evolution.evolve_star(4.0, 0.001, rng=draw.Stream(0), fe_h=0.0)
    assert value == pytest.approx(zams["luminosity_sol"], rel=0.15)
    assert 100.0 < value < 400.0


def test_the_light_rises_with_mass_at_the_disk_age(grid):
    masses = [0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.2, 1.4, 1.41, 2.0, 4.0, 10.0, 30.0]
    values = [bl.birth_luminosity_sol(m, 0.0, DISK_AGE_GY) for m in masses]
    assert all(b > a for a, b in zip(values, values[1:])), list(zip(masses, values))


def test_below_the_lightest_row_of_the_table_the_lightest_row_stands(grid):
    lightest = float(grid[0].mass_solar[0])
    assert bl.birth_luminosity_sol(lightest / 4.0, 0.0, DISK_AGE_GY) == bl.birth_luminosity_sol(lightest, 0.0, DISK_AGE_GY)


def test_an_age_younger_than_the_grid_reads_the_first_age(grid):
    assert bl.birth_luminosity_sol(0.5, 0.0, 1.0e-5) == bl.birth_luminosity_sol(0.5, 0.0, grid[0].age_gyr)


def test_above_the_table_the_star_is_the_one_stellar_evolution_gives(grid):
    state = stellar_evolution.evolve_star(2.0, DISK_AGE_GY, rng=draw.Stream(0), fe_h=0.0)
    assert bl.birth_luminosity_sol(2.0, 0.0, DISK_AGE_GY) == state["luminosity_sol"]


def test_it_draws_nothing(grid):
    draw.set_run_seed(7)
    first = draw.random()
    draw.set_run_seed(7)
    bl.birth_luminosity_sol(2.0, 0.0, DISK_AGE_GY)
    bl.birth_luminosity_sol(0.5, 0.0, DISK_AGE_GY)
    assert draw.random() == first


@pytest.mark.parametrize("mass, age", [(0.0, 0.01), (-1.0, 0.01), (1.0, 0.0)])
def test_a_mass_or_age_that_is_not_above_zero_is_refused(mass, age):
    with pytest.raises(ValueError):
        bl.birth_luminosity_sol(mass, 0.0, age)


def test_the_step_at_the_top_of_the_table_is_the_one_the_docstring_names(grid):
    below, above = bl.birth_luminosity_sol(1.4, 0.0, DISK_AGE_GY), bl.birth_luminosity_sol(1.4000001, 0.0, DISK_AGE_GY)
    assert 1.5 < above / below < 2.0 and math.isfinite(above)
