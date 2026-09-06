"""
Full star x planet x system combinatorial matrix.

test_star_matrix.py exhaustively covers every star type alone (no planets).
test_planets.py exhaustively covers every (planet class, zone) pair, but only
ever against a single G2V host. Neither combines the two, so a planet's
physics (surface temperature, atmospheric pressure, gravity clamps) has never
actually been exercised against the full range of host luminosities/masses --
from the dimmest subdwarf to the brightest hypergiant -- only against one
Sun-like star. This module closes that gap with two complementary sweeps:

1. Every (spectral class x Yerkes class) host, at a representative subclass,
   crossed with every (planet class, zone) pair `PLANET_CLASSES` declares
   valid -- built as standalone `Planet` objects (bypassing `StarSystem`'s
   orbit loop) so each cell of the matrix is under direct, explicit control
   of exactly which star type hosts exactly which class in exactly which
   zone, with no ambiguity from the loop's own zone-snapping/placement logic.

2. A full `StarSystem` for every one of those same host star types, each
   with ASTEROID_BELT and HABITABLE_WORLD both forced True and MOONS left at
   its default (None, i.e. random per-planet) -- checking every invariant
   test_systems.py already establishes, across the full star matrix rather
   than that module's ~14-type representative sample.

The subclass dimension (0-9) is deliberately NOT swept here: it mainly
perturbs temperature/luminosity continuously within a spectral class, which
test_star_matrix.py's star-only sweep already covers at full resolution;
crossing all 10 subclasses here as well would multiply this module's already
large matrix ~10x for no meaningfully different physics coverage.

Run with: pytest tests/test_full_matrix.py
"""
import math

import pytest

from stellarObjects.config import SystemConfig
from stellarObjects.planetData import Planet
from stellarObjects.starData import Star
from stellarObjects.systemData import StarSystem
from stellarObjects import program_constants as prog_c

SPECTRAL_CLASSES = ["O", "B", "A", "F", "G", "K", "M"]
YERKES_CLASSES = ["0", "IA+", "IA", "IAB", "IB", "II", "III", "IV", "V", "VI", "VII"]
REPRESENTATIVE_SUBCLASS = 5

ALL_HOST_STAR_TYPES = [
    f"{spec}{REPRESENTATIVE_SUBCLASS}{yk}"
    for spec in SPECTRAL_CLASSES
    for yk in YERKES_CLASSES
]

ZONE_CHARS = "hec"

# (class, zone) pairs PLANET_CLASSES actually declares valid -- same
# computation as test_planets.py's VALID_CLASS_ZONE_PAIRS.
VALID_CLASS_ZONE_PAIRS = [
    (cls, zone)
    for cls, data in prog_c.PLANET_CLASSES.items()
    for zone in ZONE_CHARS
    if data[zone]
]

FULL_MATRIX = [
    (star_type, cls, zone)
    for star_type in ALL_HOST_STAR_TYPES
    for cls, zone in VALID_CLASS_ZONE_PAIRS
]
FULL_MATRIX_IDS = [f"{s}-{c}-{z}" for s, c, z in FULL_MATRIX]


def make_host_star(star_type):
    cfg = SystemConfig()
    cfg.STAR_TYPE = star_type
    return Star(cfg)


def distance_for_zone(star, zone):
    inner, outer = star.habitable_zone
    if zone == 'h':
        return inner * 0.5
    if zone == 'c':
        return outer * 2.0
    return (inner + outer) / 2.0


def make_planet(host_star, cls, zone):
    distance = distance_for_zone(host_star, zone)
    cfg = SystemConfig()
    return Planet(
        cfg, host_star, host_star.habitable_zone, distance, host_star.type[0],
        host_star.luminosity, host_star.radius, host_star.temperature, host_star.mass,
        planet_class=cls,
    )


@pytest.mark.parametrize("star_type,cls,zone", FULL_MATRIX, ids=FULL_MATRIX_IDS)
def test_planet_generates_without_error_for_every_star_type(star_type, cls, zone):
    host = make_host_star(star_type)
    make_planet(host, cls, zone)


@pytest.mark.parametrize("star_type,cls,zone", FULL_MATRIX, ids=FULL_MATRIX_IDS)
def test_planet_physical_properties_stay_sane_for_every_host_star_type(star_type, cls, zone):
    """
    Every scalar physical property a planet computes must stay finite and
    positive (or, for atmosphere-related fields on an airless class, exactly
    the declared "none" sentinel) no matter how extreme the host star is --
    from a subdwarf a hundred-thousandth as luminous as the Sun to a
    hypergiant two million times as luminous.
    """
    host = make_host_star(star_type)
    planet = make_planet(host, cls, zone)

    for attr in ("mass", "density", "gravity", "surface_temperature", "volume",
                 "period", "hill_radius", "min_orbit_distance"):
        value = getattr(planet, attr)
        assert math.isfinite(value), f"{star_type}/{cls}-{zone}: {attr}={value} is not finite"
        assert value > 0, f"{star_type}/{cls}-{zone}: {attr}={value} should be positive"

    class_data = prog_c.PLANET_CLASSES[cls]
    if class_data["atmosphere"] is None:
        assert planet.atmosphere == "None"
        assert planet.atmospheric_pressure == 0.0
    else:
        assert math.isfinite(planet.atmospheric_pressure) and planet.atmospheric_pressure >= 0, (
            f"{star_type}/{cls}-{zone}: atmospheric_pressure={planet.atmospheric_pressure}"
        )
        assert math.isfinite(planet.scale_height) and planet.scale_height > 0, (
            f"{star_type}/{cls}-{zone}: scale_height={planet.scale_height}"
        )

    # Class M is explicitly clamped to Earth-like bounds regardless of host star
    # (see planetPhysics.calculate_surface_gravity/calculate_atmospheric_conditions)
    # -- this must hold even at the extremes of the host star matrix.
    if cls == "M":
        assert 0.75 <= planet.gravity <= 1.25, f"{star_type}: M-class gravity {planet.gravity} outside clamp"
        assert 90000 <= planet.atmospheric_pressure <= 112000, (
            f"{star_type}: M-class pressure {planet.atmospheric_pressure} outside clamp"
        )
        assert 283 <= planet.surface_temperature <= 290, (
            f"{star_type}: M-class surface_temperature {planet.surface_temperature} outside clamp"
        )


@pytest.mark.parametrize("star_type", ALL_HOST_STAR_TYPES)
def test_full_system_generates_with_belt_and_habitable_world_across_every_star_type(star_type):
    """
    Companion to test_systems.py's tri-state sweep, but across the *entire*
    star matrix (every spectral x Yerkes class, at a representative
    subclass) instead of that module's ~14-type representative sample --
    both guarantees (see systemData.py's dedicated fallback slots for each)
    must hold everywhere, and moons (left at the default MOONS=None, i.e.
    random per-planet) must never break anything regardless of host star.
    """
    cfg = SystemConfig()
    cfg.STAR_TYPE = star_type
    cfg.ASTEROID_BELT = True
    cfg.HABITABLE_WORLD = True
    if cfg.HABITABLE_WORLD is True and cfg.ASTEROID_BELT is True:
        cfg.LARGE_STAR = True

    system = StarSystem(system_config=cfg)

    assert str(system).strip()
    assert system.belt_count >= 1, f"{star_type}: no asteroid belt generated"
    assert system.hab_count >= 1, f"{star_type}: no habitable world generated"

    for star in system.stars:
        if star.lifespan == float('inf'):
            continue
        assert star.age <= star.lifespan * 1.001, (
            f"{star_type}: {star.name} age {star.age} Gy exceeds lifespan {star.lifespan} Gy"
        )

    planet_count, belt_count, moon_count = system.count_objects()
    assert planet_count == system.planet_count
    assert belt_count == system.belt_count
    assert moon_count == system.moon_count
    assert planet_count + belt_count == len(system.planets)

    objects = system.planets
    for i in range(1, len(objects)):
        prev, cur = objects[i - 1], objects[i]
        gap = cur.distance - prev.upper_limit if prev.type == 'a' else cur.distance - prev.distance
        if prev.type == 'a':
            min_gap = prog_c.MIN_ASTEROID_BELT_SEPARATION
        elif cur.type == 'a':
            min_gap = prev.min_orbit_distance
        else:
            min_gap = max(cur.min_orbit_distance, prev.min_orbit_distance)
        assert gap >= min_gap - 1e-9, (
            f"{star_type}: objects {i - 1} ({prev.type}) and {i} ({cur.type}) overlap: "
            f"gap={gap}, required>={min_gap}"
        )
