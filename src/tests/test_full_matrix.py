"""
Full star x planet x system combinatorial matrix.

test_star_matrix.py exhaustively covers every star type alone (no planets).
test_planets.py covers the planet draw by zone, but only
ever against a single G2V host. Neither combines the two, so a planet's
physics (surface temperature, atmospheric pressure, gravity clamps) has never
actually been exercised against the full range of host luminosities/masses --
from the dimmest subdwarf to the brightest hypergiant -- only against one
Sun-like star. This module closes that gap with two complementary sweeps:

1. Every (spectral class x Yerkes class) host, at a representative subclass,
   crossed with each zone (the planet's class is derived from what it comes out as, GEN.646) -- built as standalone `Planet` objects (bypassing `StarSystem`'s
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

from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.star import Star
from planetgen.generation.system import StarSystem
from planetgen import tuning as prog_c
from planetgen.physics.constants import G
from tests import render_support

SPECTRAL_CLASSES = ["O", "B", "A", "F", "G", "K", "M"]
YERKES_CLASSES = ["0", "IA+", "IA", "IAB", "IB", "II", "III", "IV", "V", "VI", "VII"]
REPRESENTATIVE_SUBCLASS = 5

ALL_HOST_STAR_TYPES = [
    f"{spec}{REPRESENTATIVE_SUBCLASS}{yk}"
    for spec in SPECTRAL_CLASSES
    for yk in YERKES_CLASSES
]

ZONE_CHARS = "hec"

FULL_MATRIX = [(star_type, zone) for star_type in ALL_HOST_STAR_TYPES for zone in ZONE_CHARS]
FULL_MATRIX_IDS = [f"{s}-{z}" for s, z in FULL_MATRIX]


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


def make_planet(host_star, zone):
    distance = distance_for_zone(host_star, zone)
    cfg = SystemConfig()
    return Planet(cfg, host_star, host_star.habitable_zone, distance)


@pytest.mark.parametrize("star_type,zone", FULL_MATRIX, ids=FULL_MATRIX_IDS)
def test_planet_physical_properties_stay_sane_for_every_host_star_type(star_type, zone):
    """
    Every scalar physical property a planet computes must stay finite and
    positive (or, for atmosphere-related fields on an airless class, exactly
    the declared "none" sentinel) no matter how extreme the host star is --
    from a subdwarf a hundred-thousandth as luminous as the Sun to a
    hypergiant two million times as luminous.
    """
    host = make_host_star(star_type)
    planet = make_planet(host, zone)
    cls = planet.planet_class

    for attr in ("mass", "density", "gravity", "surface_temperature", "volume",
                 "period", "hill_radius", "min_orbit_distance"):
        value = getattr(planet, attr)
        assert math.isfinite(value), f"{star_type}/{cls}-{zone}: {attr}={value} is not finite"
        assert value > 0, f"{star_type}/{cls}-{zone}: {attr}={value} should be positive"

    assert cls in prog_c.PLANET_CLASSES
    if planet.atmospheric_pressure == 0.0:
        # an airless body, or one whose host stripped its air (GEN.411: the air is solved from what the body
        # outgassed, kept and lost)
        assert planet.atmosphere == "None"
        assert planet.atmospheric_pressure == 0.0
    else:
        # The class no longer fixes a pressure range, so the check is the one physical bound every air has:
        # it cannot outweigh the planet it sits on.  That still catches a unit-conversion bug (a kPa read as
        # Pa, a bar as Pa) while letting a thin Mars-like air, a Venus-like one or a deep steam envelope
        # through; how realistic a class's pressure is belongs to test_surface_air.py.
        pressure = planet.atmospheric_pressure
        assert math.isfinite(pressure) and pressure > 0.0, (
            f"{star_type}/{cls}-{zone}: atmospheric_pressure={pressure}"
        )
        surface_gravity = G * planet.mass / (planet.radius * 1000.0) ** 2
        air_mass = pressure * 4.0 * math.pi * (planet.radius * 1000.0) ** 2 / surface_gravity
        assert air_mass < planet.mass, (
            f"{star_type}/{cls}-{zone}: the air ({air_mass:.3g} kg) outweighs the planet ({planet.mass:.3g} kg)"
        )
        assert math.isfinite(planet.scale_height) and planet.scale_height > 0, (
            f"{star_type}/{cls}-{zone}: scale_height={planet.scale_height}"
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

    `BINARY_SYSTEM` is pinned False: left at its own default it now rolls
    real, spectral-class-dependent chance (see
    `StarSystem._companion`), and an S-type (wide) pair's own
    `a_crit_au` stability ceiling can make these two guarantees
    geometrically impossible for an especially luminous host (e.g. an
    O-type supergiant's habitable zone can sit beyond any sampled
    companion's stability limit) -- binary-vs-single is incidental to what
    this test actually checks.
    """
    cfg = SystemConfig()
    cfg.STAR_TYPE = star_type
    cfg.BINARY_SYSTEM = False
    cfg.ASTEROID_BELT = True
    cfg.HABITABLE_WORLD = True
    if cfg.HABITABLE_WORLD is True and cfg.ASTEROID_BELT is True:
        cfg.LARGE_STAR = True

    system = StarSystem(system_config=cfg)

    assert render_support.text(system).strip()
    assert system.belt_count >= 1, f"{star_type}: no asteroid belt generated"
    # GEN.646: no class is drawn to fit, so a star whose ecosphere cannot hold a world answers the request "unmet" and says so
    assert system.habitable_class_count >= 1 or system.habitable_request == "unmet", \
        f"{star_type}: HABITABLE_WORLD=True produced no habitable world and did not say so"

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
        gap = cur.distance - prev.upper_limit if prev.body_type == 'a' else cur.distance - prev.distance
        if prev.body_type == 'a':
            min_gap = prog_c.MIN_ASTEROID_BELT_SEPARATION
        elif cur.body_type == 'a':
            min_gap = prev.min_orbit_distance
        else:
            min_gap = max(cur.min_orbit_distance, prev.min_orbit_distance)
        assert gap >= min_gap - 1e-9, (
            f"{star_type}: objects {i - 1} ({prev.body_type}) and {i} ({cur.body_type}) overlap: "
            f"gap={gap}, required>={min_gap}"
        )
