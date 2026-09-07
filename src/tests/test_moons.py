"""
Moon-generation regression tests.

Generates moons across a spread of host planet classes (terrestrial and gas
giant, since planetPhysics.generate_moons doesn't restrict what *type* of
planet can host moons -- only what type of planet can *be* one: type 't',
not in MOON_BLACKLIST, and small enough relative to the host) and checks the
result against the invariants generate_moons is meant to guarantee.

Run with: pytest tests/test_moons.py
"""
import math

import pytest

from stellarObjects.config import SystemConfig
from stellarObjects.planetData import Planet
from stellarObjects.starData import Star
from stellarObjects import program_constants as prog_c

# A spread of host classes: small/large terrestrial and gas giants, in
# whichever zone each actually supports (see PLANET_CLASSES h/e/c flags).
HOST_CLASSES_AND_ZONES = [
    ("M", "e"), ("D", "h"), ("C", "c"), ("J", "c"), ("S", "c"), ("I", "c"), ("A", "h"),
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


def make_host_planet(host_star, cls, zone, moon_count=None):
    cfg = SystemConfig()
    distance = distance_for_zone(host_star, zone)
    return Planet(
        cfg, host_star, host_star.habitable_zone, distance,
        planet_class=cls, moon_count=moon_count,
    )


def assert_moons_are_valid(host):
    for i, moon in enumerate(host.moons):
        assert moon.is_moon
        assert moon.planet_class not in prog_c.MOON_BLACKLIST
        assert prog_c.PLANET_CLASSES[moon.planet_class]["type"] == 't'
        assert prog_c.PLANET_CLASSES[moon.planet_class][host.zone]
        min_r, max_r = prog_c.PLANET_CLASSES[moon.planet_class]["radius_range"]
        assert min_r <= moon.radius <= max_r
        assert math.isfinite(moon.mass) and moon.mass > 0
        assert moon.mass <= host.mass / 10 * 1.001
        assert moon.radius <= host.radius / (10 ** (1 / 3)) * 1.001
        if i > 0:
            prev = host.moons[i - 1]
            assert moon.distance >= prev.distance, "moons must be generated in increasing orbital distance"


@pytest.mark.parametrize("cls,zone", HOST_CLASSES_AND_ZONES, ids=[f"{c}-{z}" for c, z in HOST_CLASSES_AND_ZONES])
@pytest.mark.parametrize("moon_count", [None, 0, 1, 3, 50])
def test_moon_generation_never_errors_and_respects_requested_count(host_star, cls, zone, moon_count):
    for _ in range(3):
        host = make_host_planet(host_star, cls, zone, moon_count=moon_count)
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
            planet_class="J",
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
            planet_class="D", moon_count=5,
        )
        assert_moons_are_valid(planet)


def test_moons_cannot_themselves_have_moons(host_star):
    cfg = SystemConfig()
    distance = distance_for_zone(host_star, "c")
    host = make_host_planet(host_star, "J", "c", moon_count=3)
    for moon in host.moons:
        assert moon.moons == []
