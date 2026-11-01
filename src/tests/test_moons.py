"""
Moon-generation regression tests.

Generates moons across a spread of host planets (terrestrial and gas giant, of the masses a caller gives, since
planetPhysics.generate_moons doesn't restrict what *type* of planet can host moons) and checks the result against the invariants
generate_moons is meant to guarantee: each moon's mass is drawn from the moon size law (GEN.646) between the lightest moon and
its parent's ceiling, and the moon is small enough for its host to hold.

Run with: pytest tests/test_moons.py
"""
import math

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.star import Star
from planetgen import tuning as prog_c
from planetgen.physics import constants, planet_kind

EARTH = constants.EARTH_MASS_TO_KG

# A spread of hosts: light and heavy terrestrial bodies and ice and gas giants, each in a zone: (mass in Earth masses, zone).
HOSTS_AND_ZONES = [
    (1.0, "e"), (0.2, "h"), (0.05, "c"), (300.0, "c"), (25.0, "c"), (3.0, "h"),
]


@pytest.fixture(scope="module")
def host_star():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    return Star(cfg)


def distance_for_zone(star, zone):
    inner, outer = star.habitable_zone
    if zone == 'h':
        return inner * 0.5
    if zone == 'c':
        return outer * 2.0
    return (inner + outer) / 2.0


def make_host_planet(host_star, mass_earth, zone, moon_count=None):
    cfg = SystemConfig()
    distance = distance_for_zone(host_star, zone)
    return Planet(
        cfg, host_star, host_star.habitable_zone, distance,
        mass=mass_earth * EARTH, moon_count=moon_count,
    )


def assert_moons_are_valid(host):
    for i, moon in enumerate(host.moons):
        assert moon.is_moon
        if moon.planet_class is None:      # GEN.1160: an Asteroid, which has no class letter
            assert moon.roundness_ratio < 1.0 and moon.body_type == "t"
        else:
            assert moon.planet_class in prog_c.PLANET_CLASSES
            shown = prog_c.PLANET_CLASSES[moon.planet_class]
            assert moon.body_type in (shown["type"], shown.get("also_type"))  # GEN.1142: a rocky moon with a hydrogen air is class T
        assert moon.zone == host.zone
        assert math.isfinite(moon.mass) and moon.mass > 0
        ceiling = planet_kind.moon_mass_ceiling_kg(host.mass, host.body_type == "g")
        assert prog_c.MOON_MIN_MASS_KG * 0.999 <= moon.mass <= ceiling * 1.001
        assert moon.mass <= host.mass / 10 * 1.001
        assert moon.radius <= host.radius / (10 ** (1 / 3)) * 1.001
        if i > 0:
            prev = host.moons[i - 1]
            assert moon.distance >= prev.distance, "moons must be generated in increasing orbital distance"


@pytest.mark.parametrize("mass_earth,zone", HOSTS_AND_ZONES, ids=[f"{m}-{z}" for m, z in HOSTS_AND_ZONES])
@pytest.mark.parametrize("moon_count", [None, 0, 1, 3, 50])
def test_moon_generation_never_errors_and_respects_requested_count(host_star, mass_earth, zone, moon_count):
    for _ in range(3):
        host = make_host_planet(host_star, mass_earth, zone, moon_count=moon_count)
        if moon_count is not None:
            assert len(host.moons) <= moon_count, "generate_moons must never exceed the requested count"
        if moon_count == 0:
            assert host.moons == []
        assert_moons_are_valid(host)


@pytest.mark.parametrize("moons_flag", [True, False, None])
def test_system_config_moons_flag(host_star, moons_flag):
    cfg = SystemConfig()
    cfg.MOONS = moons_flag
    distance = distance_for_zone(host_star, "c")
    saw_moons = False
    saw_none = False
    for _ in range(20):
        planet = Planet(
            cfg, host_star, host_star.habitable_zone, distance,
            mass=300 * EARTH,
        )
        if moons_flag is False:
            assert planet.moons == []
        elif planet.moons:
            saw_moons = True
        else:
            saw_none = True
        assert_moons_are_valid(planet)

    if moons_flag is None:
        # Should see both outcomes over enough trials (50/50 chance per planet).
        assert saw_moons and saw_none


def test_tiny_host_planet_gets_no_room_for_moons_without_error(host_star):
    """
    A host planet small/light enough that its min_orbit_distance leaves no
    orbital room should simply end up with zero moons, not raise.
    """
    cfg = SystemConfig()
    distance = distance_for_zone(host_star, "h")
    for _ in range(10):
        planet = Planet(
            cfg, host_star, host_star.habitable_zone, distance,
            mass=prog_c.PLANET_MIN_MASS_KG, moon_count=5,
        )
        assert_moons_are_valid(planet)


def test_moons_cannot_themselves_have_moons(host_star):
    host = make_host_planet(host_star, 300.0, "c", moon_count=3)
    for moon in host.moons:
        assert moon.moons == []


def test_a_moons_hill_sphere_is_about_its_planet_not_the_star(host_star):
    """GEN.138: the Moon's Hill sphere about the Earth is about 61,500 km, not the 888 km the Sun's mass gave."""
    from planetgen.physics.planets import update_hill_sphere

    host = make_host_planet(host_star, 300.0, "c", moon_count=3)
    assert host.moons
    for moon in host.moons:
        assert moon.primary_mass_kg == host.mass
        expected_km = moon.distance * constants.AU_TO_M * (moon.mass / (3 * (host.mass + moon.mass))) ** (1 / 3) / 1000
        assert moon.hill_radius == pytest.approx(expected_km, rel=1e-9)
        assert moon.min_orbit_distance == pytest.approx(5 * moon.hill_radius / constants.AU_TO_KM, rel=1e-9)

    # The Earth-Moon pair, by hand: 384,400 km out, 7.342e22 kg about 5.972e24 kg.
    moon = host.moons[0]
    moon.distance, moon.mass, moon.primary_mass_kg = 384_400e3 / constants.AU_TO_M, 7.342e22, 5.972e24
    update_hill_sphere(moon)
    assert moon.hill_radius == pytest.approx(61_500, rel=0.01)
