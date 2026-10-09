"""
Moon stability helpers (TODO TEST.33).

Direct tests of `planets.moon_orbit_bounds_km`,
`planets.drop_unstable_moons` and `planets.update_hill_sphere`,
and of the "No valid planet class ..." / "Invalid planet class for this
zone" errors `generate_planet_properties` raises when a requested radius,
mass or class fits nothing. Behavior checks only; reference values live in
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
    return SimpleNamespace(distance=distance_km / constants.AU_TO_KM, radius=radius, mass=mass)


# --- moon_orbit_bounds_km -------------------------------------------------

def test_bounds_without_atmosphere_use_a_100_km_margin():
    low, high = planets.moon_orbit_bounds_km(bare_planet(radius=6000.0, scale_height=None, hill_radius=2e6))
    assert low == pytest.approx(6000.0 + 6000.0 / CUBE_ROOT_10 + 100)
    assert high == pytest.approx(2e6 * tuning.MOON_PROGRADE_STABLE_HILL_FRACTION)


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
    planet = Planet(star.system_config, star, star.habitable_zone, cold_distance(star), planet_class="J", moon_count=6)
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
    planet = Planet(star.system_config, star, star.habitable_zone, cold_distance(star), planet_class="J", moon_count=6)
    count = len(planet.moons)
    assert planets.drop_unstable_moons(planet) == 0
    assert len(planet.moons) == count


def test_shrinking_the_hill_sphere_drops_the_outer_moons_first(star):
    planet = Planet(star.system_config, star, star.habitable_zone, cold_distance(star), planet_class="J", moon_count=6)
    distances = sorted(moon.distance for moon in planet.moons)
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
    planet = Planet(star.system_config, star, star.habitable_zone, cold_distance(star), planet_class="J", moon_count=0)
    before = planet.hill_radius
    planet.distance *= 2
    assert planet.hill_radius == before  # stale until refreshed
    planets.update_hill_sphere(planet)
    assert planet.hill_radius == pytest.approx(2 * before)


def test_update_hill_sphere_at_zero_distance_is_zero():
    planet = hill_stub(0.0, 6e24)
    planets.update_hill_sphere(planet)
    assert planet.hill_radius == 0 and planet.min_orbit_distance == 0


# --- "no valid planet class" errors ---------------------------------------

@pytest.mark.parametrize("kwargs, message", [
    ({"radius": 1e9}, "No valid planet class for the given radius in this zone"),
    ({"radius": 1.0}, "No valid planet class for the given radius in this zone"),
    ({"mass": 1.0}, "No valid planet class for the given mass in this zone"),
    ({"mass": 1e40}, "No valid planet class for the given mass in this zone"),
    ({"radius": 6000.0, "mass": 1e30}, "No valid planet class for the given radius/mass in this zone"),
    ({"radius": 200000.0, "mass": 1e22}, "No valid planet class for the given radius/mass in this zone"),
])
def test_no_valid_class_for_given_radius_or_mass(star, kwargs, message):
    with pytest.raises(ValueError, match=message):
        Planet(star.system_config, star, star.habitable_zone, cold_distance(star), **kwargs)


def test_no_valid_class_when_only_habitable_classes_fit_and_they_are_barred(star):
    """12,000 km in the ecosphere fits only class V, which is habitable;
    with HABITABLE_WORLD=False nothing is left."""
    star.system_config.HABITABLE_WORLD = False
    inner, outer = star.habitable_zone
    with pytest.raises(ValueError, match="No valid planet class for the given radius in this zone"):
        Planet(star.system_config, star, star.habitable_zone, (inner + outer) / 2, radius=12000.0)


def test_radius_that_fits_a_class_elsewhere_but_not_in_this_zone(star):
    """52,000 km fits classes J and T; T is cold-zone only, so in the hot
    zone only J is left -- but 12,000 km (class V, ecosphere only) fits
    nothing hot at all."""
    hot = star.habitable_zone[0] * 0.5
    assert Planet(star.system_config, star, star.habitable_zone, hot, radius=52000.0).planet_class == "J"
    with pytest.raises(ValueError, match="No valid planet class for the given radius in this zone"):
        Planet(star.system_config, star, star.habitable_zone, hot, radius=12000.0)


@pytest.mark.parametrize("planet_class", ["Z", "M", "E"])
def test_invalid_class_for_the_zone(star, planet_class):
    """Unknown class 'Z', and ecosphere-only classes in the cold zone."""
    with pytest.raises(ValueError, match="Invalid planet class for this zone"):
        Planet(star.system_config, star, star.habitable_zone, cold_distance(star), planet_class=planet_class)


def test_class_and_radius_mismatch_raises(star):
    with pytest.raises(ValueError, match="Invalid radius for planet class"):
        Planet(star.system_config, star, star.habitable_zone, cold_distance(star), planet_class="J", radius=1000.0)


def test_class_and_mass_mismatch_raises(star):
    with pytest.raises(ValueError, match="Invalid mass for planet class"):
        Planet(star.system_config, star, star.habitable_zone, cold_distance(star), planet_class="J", mass=1.0)


def test_mass_only_planet_gets_a_class_and_a_radius_in_range(star):
    mass = 1e27  # only class J's mass range reaches this
    planet = Planet(star.system_config, star, star.habitable_zone, cold_distance(star), mass=mass)
    low, high = tuning.PLANET_CLASSES[planet.planet_class]["radius_range"]
    assert planet.planet_class == "J"
    assert low <= planet.radius <= high
    assert math.isfinite(planet.hill_radius) and planet.hill_radius > 0
