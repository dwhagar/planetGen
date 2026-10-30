"""
Moon stability helpers (TODO TEST.33).

Direct tests of `planets.moon_orbit_bounds_km`,
`planets.drop_unstable_moons` and `planets.update_hill_sphere`,
and of the ValueErrors `generate_planet_properties` raises when a requested
radius or mass is out of bounds or implies an impossible density. Behavior checks only; reference values live in
the physics reference gate.
"""
import math
from types import SimpleNamespace

import pytest

from planetgen.physics import constants, planets
from planetgen import tuning
from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.star import Star

from tests.fuzz_support import deterministic_entropy

CUBE_ROOT_10 = 10 ** (1 / 3)
GAS_GIANT_KG = 300 * constants.EARTH_MASS_TO_KG


@pytest.fixture(autouse=True)
def _seeded():
    with deterministic_entropy(3301):
        yield


@pytest.fixture
def star():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    return Star(cfg)


def cold_distance(star):
    return star.habitable_zone[1] * 3


def bare_planet(radius=10000.0, scale_height=None, hill_radius=1e6, moons=(), mass=6e24):
    """A stand-in carrying only what the moon helpers read."""
    return SimpleNamespace(radius=radius, scale_height=scale_height, hill_radius=hill_radius, moons=list(moons),
                           mass=mass)


def moon_at(distance_km, radius=100.0, mass=1e18):
    return SimpleNamespace(distance=distance_km / constants.AU_TO_KM, radius=radius, mass=mass, eccentricity=0.0)


# --- moon_orbit_bounds_km -------------------------------------------------

def test_bounds_without_atmosphere_use_a_100_km_margin():
    low, high = planets.moon_orbit_bounds_km(bare_planet(radius=6000.0, scale_height=None, hill_radius=2e6))
    assert low == pytest.approx(6000.0 + 6000.0 / CUBE_ROOT_10 + 100)
    assert high == pytest.approx(2e6 * 0.4895)


@pytest.mark.parametrize("scale_height", [0, 0.0])
def test_zero_scale_height_counts_as_no_atmosphere(scale_height):
    with_zero = planets.moon_orbit_bounds_km(bare_planet(scale_height=scale_height))
    without = planets.moon_orbit_bounds_km(bare_planet(scale_height=None))
    assert with_zero == without


def test_bounds_with_atmosphere_use_15_scale_heights():
    low, _ = planets.moon_orbit_bounds_km(bare_planet(radius=6000.0, scale_height=8.5))
    assert low == pytest.approx(6000.0 + 6000.0 / CUBE_ROOT_10 + 15 * 8.5)


def test_inner_bound_clears_the_planet_and_the_rigid_roche_limit():
    for radius in (50.0, 6000.0, 70000.0):
        low, _ = planets.moon_orbit_bounds_km(bare_planet(radius=radius))
        assert low > 1.26 * radius


def test_outer_bound_grows_with_the_hill_radius():
    highs = [planets.moon_orbit_bounds_km(bare_planet(hill_radius=h))[1] for h in (1e3, 1e5, 1e7)]
    assert highs == sorted(highs) and highs[0] < highs[-1]


def test_a_tiny_hill_sphere_leaves_no_room():
    low, high = planets.moon_orbit_bounds_km(bare_planet(radius=6000.0, hill_radius=1000.0))
    assert low >= high


def test_real_gas_giant_bounds_are_ordered_and_hold_its_moons(star):
    planet = Planet(star.system_config, star, star.habitable_zone, cold_distance(star), mass=GAS_GIANT_KG, moon_count=6)
    low, high = planets.moon_orbit_bounds_km(planet)
    assert isinstance(low, float) and isinstance(high, float)
    assert planet.radius < low < high
    assert planet.moons
    for moon in planet.moons:
        assert low <= moon.distance * constants.AU_TO_KM <= high


# --- drop_unstable_moons --------------------------------------------------

def test_drop_with_no_moons_returns_zero():
    planet = bare_planet()
    assert planets.drop_unstable_moons(planet) == 0
    assert planet.moons == []


def test_drop_removes_moons_outside_the_bounds_and_counts_them():
    planet = bare_planet(radius=6000.0, hill_radius=1e6)
    low, high = planets.moon_orbit_bounds_km(planet)
    inside = moon_at((low + high) / 2)
    too_close = moon_at(low * 0.5)
    too_far = moon_at(high * 1.5)
    planet.moons[:] = [too_close, inside, too_far]
    moons_list = planet.moons
    assert planets.drop_unstable_moons(planet) == 2
    assert planet.moons == [inside]
    assert planet.moons is moons_list  # trimmed in place


def test_drop_keeps_moons_exactly_on_either_bound():
    planet = bare_planet(radius=6000.0, hill_radius=1e6)
    low, high = planets.moon_orbit_bounds_km(planet)
    edge_moons = [moon_at(low * (1 + 1e-12)), moon_at(high * (1 - 1e-12))]
    planet.moons[:] = edge_moons
    assert planets.drop_unstable_moons(planet) == 0
    assert planet.moons == edge_moons


def test_drop_removes_a_moon_too_large_for_the_planet():
    planet = bare_planet(radius=6000.0, hill_radius=1e6)
    low, high = planets.moon_orbit_bounds_km(planet)
    max_radius = 6000.0 / CUBE_ROOT_10
    fits = moon_at((low + high) / 2, radius=max_radius)
    too_big = moon_at((low + high) / 2, radius=max_radius * 1.01)
    planet.moons[:] = [fits, too_big]
    assert planets.drop_unstable_moons(planet) == 1
    assert planet.moons == [fits]


def test_drop_removes_a_moon_too_heavy_for_the_planet():
    planet = bare_planet(radius=6000.0, hill_radius=1e6, mass=6e24)
    low, high = planets.moon_orbit_bounds_km(planet)
    fits = moon_at((low + high) / 2, mass=6e23)
    too_heavy = moon_at((low + high) / 2, mass=6.01e23)
    planet.moons[:] = [fits, too_heavy]
    assert planets.drop_unstable_moons(planet) == 1
    assert planet.moons == [fits]


def test_drop_everything_when_there_is_no_room():
    planet = bare_planet(radius=6000.0, hill_radius=1000.0)
    planet.moons[:] = [moon_at(d) for d in (500.0, 10000.0, 1e6)]
    assert planets.drop_unstable_moons(planet) == 3
    assert planet.moons == []


def test_drop_is_idempotent_on_a_freshly_generated_planet(star):
    planet = Planet(star.system_config, star, star.habitable_zone, cold_distance(star), mass=GAS_GIANT_KG, moon_count=6)
    count = len(planet.moons)
    assert planets.drop_unstable_moons(planet) == 0
    assert len(planet.moons) == count


def test_shrinking_the_hill_sphere_drops_the_outer_moons_first(star):
    # (an eccentric planet holds its moons closer in, so the first draws can leave fewer than two: draw again)
    for _attempt in range(12):
        planet = Planet(star.system_config, star, star.habitable_zone, cold_distance(star), mass=GAS_GIANT_KG, moon_count=6)
        distances = sorted(moon.distance for moon in planet.moons)
        if len(distances) >= 2:
            break
    assert len(distances) >= 2
    planet.distance /= 50  # drag it far inward: its Hill sphere shrinks 50x
    planets.update_hill_sphere(planet)
    dropped = planets.drop_unstable_moons(planet)
    assert dropped >= 1
    kept = sorted(moon.distance for moon in planet.moons)
    assert kept == distances[:len(kept)]


# --- update_hill_sphere ---------------------------------------------------

def hill_stub(distance_au, mass_kg, star_mass_kg=constants.SOLAR_MASS_TO_KG):
    return SimpleNamespace(distance=distance_au, mass=mass_kg, star=SimpleNamespace(mass=star_mass_kg))


def test_update_hill_sphere_sets_both_fields_consistently():
    planet = hill_stub(1.0, 6e24)
    planets.update_hill_sphere(planet)
    assert planet.hill_radius > 0
    assert planet.min_orbit_distance == pytest.approx(5 * planet.hill_radius / constants.AU_TO_KM)


def test_hill_radius_is_linear_in_distance():
    near, far = hill_stub(1.0, 6e24), hill_stub(4.0, 6e24)
    planets.update_hill_sphere(near)
    planets.update_hill_sphere(far)
    assert far.hill_radius == pytest.approx(4 * near.hill_radius)
    assert far.min_orbit_distance == pytest.approx(4 * near.min_orbit_distance)


def test_hill_radius_scales_with_the_cube_root_of_mass_ratio():
    light, heavy = hill_stub(1.0, 1e24), hill_stub(1.0, 8e24)
    heavy_star = hill_stub(1.0, 8e24, constants.SOLAR_MASS_TO_KG * 8)
    for planet in (light, heavy, heavy_star):
        planets.update_hill_sphere(planet)
    # The pair's total mass (GEN.138) differs from the star's alone by parts per million.
    assert heavy.hill_radius == pytest.approx(2 * light.hill_radius, rel=1e-4)
    assert heavy_star.hill_radius == pytest.approx(light.hill_radius, rel=1e-4)


def test_update_hill_sphere_tracks_a_moved_planet(star):
    planet = Planet(star.system_config, star, star.habitable_zone, cold_distance(star), mass=GAS_GIANT_KG, moon_count=0)
    before = planet.hill_radius
    planet.distance *= 2
    assert planet.hill_radius == before  # stale until refreshed
    planets.update_hill_sphere(planet)
    assert planet.hill_radius == pytest.approx(2 * before)


def test_update_hill_sphere_at_zero_distance_is_zero():
    planet = hill_stub(0.0, 6e24)
    planets.update_hill_sphere(planet)
    assert planet.hill_radius == 0 and planet.min_orbit_distance == 0


# --- errors for a radius or mass that cannot be a planet -----------------

@pytest.mark.parametrize("kwargs, message", [
    ({"radius": 1e9}, "The given radius"),
    ({"radius": 0.1}, "The given radius"),
    ({"mass": 1.0}, "The given mass"),
    ({"mass": 1e40}, "The given mass"),
    ({"radius": 6000.0, "mass": 2e28}, "density"),
    ({"radius": 200000.0, "mass": 1e22}, "density"),
])
def test_a_radius_or_mass_that_cannot_be_a_planet_raises(star, kwargs, message):
    with pytest.raises(ValueError, match=message):
        Planet(star.system_config, star, star.habitable_zone, cold_distance(star), **kwargs)


def test_a_radius_alone_is_kept_and_gets_a_class_derived_from_the_body(star):
    planet = Planet(star.system_config, star, star.habitable_zone, cold_distance(star), radius=60000.0)
    assert planet.radius == 60000.0
    assert planet.body_type == "g"


def test_mass_only_planet_gets_a_class_and_a_radius(star):
    planet = Planet(star.system_config, star, star.habitable_zone, cold_distance(star), mass=1e27)
    assert planet.mass == 1e27
    assert planet.body_type == "g"
    assert planet.planet_class in tuning.PLANET_CLASSES
    assert math.isfinite(planet.hill_radius) and planet.hill_radius > 0
    assert 10000.0 <= planet.radius <= 200000.0
