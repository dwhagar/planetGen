# tests/test_giant_settle.py

"""
GEN.961 and GEN.976: a giant sized from the radius tables is drawn when the star's age is provisional, and sized again when the age and the
orbit are final (`planets.settle_giant`, the "giant settle pass" of `StarSystem`). The unit rules of `giant_draw.settle` are in
test_giant_draw.py; here the planet it is wired to and the finished system.
"""

import math
from types import SimpleNamespace

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.generation.system import StarSystem
from planetgen.physics import constants, giant_draw, giant_radius, interior_heat, orbit_elements, planets, tidal_heating
from planetgen.util import draw
from tests import class_support

GIANT_KG = 100.0 * constants.EARTH_MASS_TO_KG
SYSTEMS = 40
SEED = 7


@pytest.fixture()
def host_star():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    with draw.bound("giant-settle:star"):
        star = Star(cfg)
    star.age = 1.0       # drawn young, so a raise of a few Gy stays inside the radius tables' age range
    return star


def _table_giant(star, seeds=40):
    for seed in range(seeds):
        planet = class_support.draw_planet(star, 0.4, seed, SystemConfig(), mass=GIANT_KG, zone_override="h", moon_count=0)
        if getattr(planet, giant_draw.SETTLE_ATTRIBUTE, None) is not None:
            return planet
    pytest.fail("no giant sized from the tables in the draws")


def _table_radius_km(planet):
    drawn = getattr(planet, giant_draw.SETTLE_ATTRIBUTE)
    mass_earth = planet.mass / constants.EARTH_MASS_TO_KG
    core = giant_radius.table_core_mass_earth(mass_earth, drawn.heavy_earth, drawn.envelope_z)
    distance_km = orbit_elements.flux_distance_au(planet.distance, planet.eccentricity) * constants.AU_TO_KM
    light = planet.star.luminosity / (4.0 * math.pi * distance_km ** 2) / 1e6
    return giant_radius.radius_km(mass_earth, core, light / giant_radius.flux_w_m2(1.0), planet.star.age)


def test_settling_a_giant_sets_its_radius_volume_gravity_and_air_together(host_star):
    planet = _table_giant(host_star)
    young = planet.radius
    host_star.age = 7.0
    assert planets.settle_giant(planet) is True
    assert planet.radius < young
    assert planet.radius == pytest.approx(_table_radius_km(planet), rel=1e-9)
    assert planet.volume == pytest.approx(4 / 3 * math.pi * planet.radius ** 3, rel=1e-12)
    assert planet.gravity == pytest.approx(constants.G * planet.mass / (planet.radius * 1000) ** 2 / constants.EARTH_GRAVITY, rel=1e-6)
    assert planet.density == pytest.approx(planet.mass / (4 / 3 * math.pi * (planet.radius * 1000) ** 3) / 1000, rel=1e-9)
    assert planet.surface_temperature is not None and planet.scale_height is not None
    assert planet.planet_class is not None


def test_a_settled_giant_is_settled_for_good(host_star):
    planet = _table_giant(host_star)
    host_star.age = 7.0
    planets.settle_giant(planet)
    stored = (planet.radius, planet.gravity, planet.scale_height, planet.composition)
    assert planets.settle_giant(planet) is False
    assert (planet.radius, planet.gravity, planet.scale_height, planet.composition) == stored


def test_settling_draws_nothing_from_the_run_stream(host_star):
    planet = _table_giant(host_star)
    host_star.age = 7.0
    draw.set_run_seed(5)
    planets.settle_giant(planet)
    after = draw.random()
    draw.set_run_seed(5)
    assert after == draw.random()


@pytest.mark.slow
def test_every_table_giant_of_a_finished_system_has_the_radius_of_its_final_age_and_orbit():
    draw.set_run_seed(SEED)
    checked = 0
    for _ in range(SYSTEMS):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        for planet in StarSystem(system_config=cfg).all_planets:
            if getattr(planet, giant_draw.SETTLE_ATTRIBUTE, None) is None:
                continue
            checked += 1
            assert planet.radius == pytest.approx(_table_radius_km(planet), rel=giant_draw.SETTLE_TOLERANCE * 1.01)
    assert checked, "the sweep made no giant sized from the tables"


def test_a_moon_the_smaller_giant_cannot_hold_is_dropped(host_star, monkeypatch):
    """The moons here are bare stand-ins, so the history of the moons (`moon_history.apply`, GEN.959) is left out: it would carry them
    through the tide and could remove the kept one for a reason of its own, and what is tested is the size limit of the smaller giant."""
    monkeypatch.setattr(planets.moon_history, "apply", lambda planet, primary_mass_kg: None)
    planet = _table_giant(host_star)
    limit_km, _ = planets.moon_size_limits(planet)
    kept = SimpleNamespace(radius=0.1 * limit_km, mass=1e-6 * planet.mass, distance=planet.hill_radius * 0.1 / constants.AU_TO_KM,
                           eccentricity=0.0, position_x=0.0, position_y=0.0, position_z=0.0)
    too_big = SimpleNamespace(radius=0.999 * limit_km, mass=1e-6 * planet.mass, distance=planet.hill_radius * 0.12 / constants.AU_TO_KM,
                              eccentricity=0.0, position_x=0.0, position_y=0.0, position_z=0.0)
    planet.moons = [kept, too_big]
    host_star.age = 7.0
    assert planets.settle_giant(planet) is True
    assert planet.moons == [kept]


def test_a_moons_tide_does_not_read_its_parents_radius(host_star):
    """GEN.987: the tide a moon feels is the host's mass, the moon's distance and its eccentricity (`tidal_heating.forcing_of`) and the
    heat it makes is a function of the moon's own mass and radius, so a giant that settles to a smaller radius leaves its moons' heat as it was."""
    for seed in range(60):
        planet = class_support.draw_planet(host_star, 5.2, seed, SystemConfig(), mass=GIANT_KG, zone_override="c", moon_count=6)
        moons = [moon for moon in planet.moons if (moon.eccentricity or 0.0) > 0.0]
        if moons:
            break
    else:
        pytest.fail("no giant with an eccentric moon in 60 draws")
    moon = moons[0]
    before = (tidal_heating.forcing_of(moon), interior_heat.update_heat(moon))
    planet.radius *= 0.9
    assert (tidal_heating.forcing_of(moon), interior_heat.update_heat(moon)) == before


# --- GEN.985: the tides of the orbit are taken again at the final age -----------------------------------------------------------------------

def _close_in_planet(star, mass_kg=5.0 * constants.EARTH_MASS_TO_KG, distance=0.1, moon_count=0):
    """A planet drawn at a star of 0.5 Gyr that kept an eccentricity of at least 0.05 and that 13.8 Gyr would wear down."""
    star.age = 0.5
    for seed in range(80):
        planet = class_support.draw_planet(star, distance, seed, SystemConfig(), mass=mass_kg, zone_override="h", moon_count=moon_count)
        birth = getattr(planet, "_orbit_birth", None)
        if birth is None or planet.eccentricity < 0.05:
            continue
        time_yr = orbit_elements.circularisation_time_yr(planet, planet.primary_mass_kg, planet.distance)
        if planet.eccentricity - orbit_elements.circularised(birth, 13.8e9, time_yr) > 0.01:
            return planet
    pytest.fail("no close-in planet that 13.8 Gyr would circularise in 80 draws")


def test_a_close_in_planet_is_circularised_for_the_final_age(host_star):
    planet = _close_in_planet(host_star)
    drawn = planet.eccentricity
    host_star.age = 13.8
    assert planets.settle_orbit(planet) is True
    time_yr = orbit_elements.circularisation_time_yr(planet, planet.primary_mass_kg, planet.distance)
    assert planet.eccentricity == pytest.approx(orbit_elements.circularised(planet._orbit_birth, 13.8e9, time_yr), rel=1e-9)
    assert planet.eccentricity < drawn
    assert planet.orbit_from_vector()["eccentricity"] == pytest.approx(planet.eccentricity, abs=1e-7)
    assert planet.orbit_from_vector()["semi_major_axis"] == pytest.approx(planet.distance, rel=1e-9)


def test_a_planet_whose_age_did_not_move_keeps_its_orbit(host_star):
    planet = _close_in_planet(host_star)
    stored = (planet.eccentricity, planet.position_x, planet.position_y, planet.position_z)
    assert planets.settle_orbit(planet) is False
    assert (planet.eccentricity, planet.position_x, planet.position_y, planet.position_z) == stored


def test_the_moons_of_a_settled_planet_are_told_its_new_eccentricity(host_star):
    planet = _close_in_planet(host_star)
    planet.moons = [SimpleNamespace(host_eccentricity=planet.eccentricity), SimpleNamespace(host_eccentricity=planet.eccentricity)]
    host_star.age = 13.8
    assert planets.settle_orbit(planet) is True
    assert all(moon.host_eccentricity == planet.eccentricity for moon in planet.moons)
