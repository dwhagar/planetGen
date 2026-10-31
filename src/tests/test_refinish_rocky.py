# tests/test_refinish_rocky.py

"""
GEN.992 (found by GEN.986): a rocky planet and its moons are made again at the star's final age (`planets.refinish_rocky`), from streams of
their own, and `StarSystem` judges the forced requests after it.
"""

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.physics import constants, planets
from planetgen.util import draw
from tests import class_support

EARTH_KG = constants.EARTH_MASS_TO_KG


@pytest.fixture()
def host_star():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    with draw.bound("refinish:star"):
        star = Star(cfg)
    star.age = 0.5
    return star


def _rocky(star, distance=1.0, moon_count=0, mass_earth=1.0, seeds=40):
    for seed in range(seeds):
        planet = class_support.draw_planet(star, distance, seed, SystemConfig(), mass=mass_earth * EARTH_KG, moon_count=moon_count)
        if planet.body_type == "t":
            return planet
    pytest.fail("no rocky planet in the draws")


def test_a_planet_made_at_the_final_age_is_not_made_again(host_star):
    planet = _rocky(host_star)
    assert planets.age_moved(planet) is False
    assert planets.refinish_rocky(planet) == 0


def test_a_small_rise_of_the_age_leaves_the_planet(host_star):
    planet = _rocky(host_star)
    host_star.age = 0.5 * (1.0 + planets.AGE_SETTLE_TOLERANCE * 0.5)
    assert planets.age_moved(planet) is False
    host_star.age = 0.5 * (1.0 + planets.AGE_SETTLE_TOLERANCE * 2.0)
    assert planets.age_moved(planet) is True


def test_a_rocky_planet_is_made_again_at_the_final_age_once(host_star):
    planet = _rocky(host_star)
    dose = planet.dose_ground_msv_yr
    host_star.age = 8.0
    assert planets.refinish_rocky(planet) == 1
    assert planets.age_moved(planet) is False and planet._age_drawn == 8.0
    assert planet.planet_class is not None
    if dose:
        assert planet.dose_ground_msv_yr != dose
    assert planets.refinish_rocky(planet) == 0


def test_the_same_planet_is_made_again_the_same_way_and_the_run_stream_stays(host_star):
    first = []
    for _ in range(2):
        planet = _rocky(host_star)
        host_star.age = 8.0
        draw.set_run_seed(31)
        planets.refinish_rocky(planet)
        first.append((planet.planet_class, planet.atmospheric_pressure, planet.surface_temperature, planet.rotation_period_hours,
                      planet.magnetic_moment_a_m2, draw.random()))
        host_star.age = 0.5
    assert first[0] == first[1]
    draw.set_run_seed(31)
    assert first[0][-1] == draw.random()


def test_the_moons_of_a_giant_are_made_again_when_the_age_rose(host_star):
    for seed in range(60):
        giant = class_support.draw_planet(host_star, 5.2, seed, SystemConfig(), mass=300.0 * EARTH_KG, zone_override="c", moon_count=4)
        if giant.moons:
            break
    else:
        pytest.fail("no giant with a moon in 60 draws")
    host_star.age = 8.0
    assert planets.refinish_rocky(giant) == len(giant.moons)
    assert all(moon._age_drawn == 8.0 for moon in giant.moons)


def test_a_made_again_planet_and_its_moons_keep_the_birth_spin_they_were_born_with(host_star):
    # GEN.961: the history is run again from the birth spin (`generate_spin(..., birth=)`), not from a new draw out of the settle stream
    planet = _rocky(host_star, distance=0.5, moon_count=2)
    born = {id(body): body._birth_spin for body in (planet, *planet.moons) if getattr(body, "_birth_spin", None) is not None}
    assert born
    host_star.age = 8.0
    assert planets.refinish_rocky(planet) >= 1
    for body in (planet, *planet.moons):
        if id(body) in born:
            assert body._birth_spin == born[id(body)]


def test_a_moon_the_planets_history_replaced_is_made_again_at_the_final_age_too():
    # GEN.992 after GEN.959: the planet's spin history can swap a moon for a new body, which has to be made again like the rest (seed 38 of the final-age sweep)
    from planetgen.generation.system import StarSystem
    with draw.bound(38):
        system = StarSystem(system_config=SystemConfig())
    for planet in system.all_planets:
        if planet.body_type == "a":
            continue
        for moon in planet.moons:
            if moon.body_type != "g":
                assert not planets.age_moved(moon), (moon.name, moon._age_drawn, moon.star.age)
