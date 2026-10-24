# tests/test_main_sequence.py

"""
A main-sequence star's luminosity, temperature and radius follow its age and birth [Fe/H] from the MIST isochrones
(GEN.652), so the Sun is faint at birth and near one solar luminosity at its present age, and the habitable zone it
gives moves outward as the star ages.
"""

import math

from planetgen.physics import constants, main_sequence, orbits
from planetgen.physics import stellar_evolution as evolution

SUN_AGE_GY = 4.57


def test_the_sun_is_near_one_solar_luminosity_now_and_about_thirty_percent_fainter_at_birth():
    now = main_sequence.state(1.0, SUN_AGE_GY, 0.0)
    birth = main_sequence.state(1.0, 0.0, 0.0)
    assert 0.95 <= now.luminosity_sol <= 1.15
    assert 0.65 <= birth.luminosity_sol <= 0.8
    assert 0.65 <= birth.luminosity_sol / now.luminosity_sol <= 0.75
    assert 5500.0 <= now.temperature_k <= 6000.0 and 0.95 <= now.radius_sol <= 1.1


def test_a_main_sequence_star_brightens_and_swells_as_it_ages():
    rows = [main_sequence.state(1.0, age, 0.0) for age in (0.1, 1.0, 3.0, 5.0, 7.0, 9.0)]
    for earlier, later in zip(rows, rows[1:]):
        assert later.luminosity_sol > earlier.luminosity_sol
        assert later.radius_sol > earlier.radius_sol


def test_a_metal_poor_star_is_hotter_and_brighter_than_a_metal_rich_one_of_the_same_mass_and_age():
    poor, rich = main_sequence.state(1.0, 4.57, -1.0), main_sequence.state(1.0, 4.57, 0.3)
    assert poor.temperature_k > rich.temperature_k
    assert poor.luminosity_sol > rich.luminosity_sol * 0.9


def test_the_stars_own_record_reads_the_same_table():
    for age in (0.05, 4.57, 8.0):
        record = evolution.evolve_star(1.0, age, fe_h=0.0)
        row = main_sequence.state(1.0, age, 0.0)
        assert math.isclose(record["luminosity_sol"], row.luminosity_sol, rel_tol=1e-6)
        assert math.isclose(record["temperature_k"], row.temperature_k, rel_tol=1e-6)


def test_the_habitable_zone_of_a_sun_like_star_moves_outward_with_its_age():
    def zone(age):
        record = evolution.evolve_star(1.0, age, fe_h=0.0)
        return orbits.calculate_habitable_zone(record["luminosity_sol"] * constants.SOLAR_LUMINOSITY,
                                               record["temperature_k"])
    young, now, old = zone(0.05), zone(SUN_AGE_GY), zone(9.0)
    for edge in (0, 1):
        assert young[edge] < now[edge] < old[edge]
    assert 0.9 <= now[0] <= 1.05 and 1.6 <= now[1] <= 1.85
    assert young[1] < 1.6
