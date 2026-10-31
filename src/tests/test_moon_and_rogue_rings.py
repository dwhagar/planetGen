# tests/test_moon_and_rogue_rings.py

"""
GEN.836, GEN.1057: a moon has the ring its history made by the law of GEN.732 (`planetgen.physics.rings`), a rogue planet has none (its moons and comets are not followed). A moon's
ring, when it has one, ends inside the Roche limit of its material.
"""

import math

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.generation.system import StarSystem
from planetgen.physics import constants, rings
from planetgen.util import draw


def _roche_ceiling_km(body_radius_km, body_mass_kg):
    """The outer edge no ring can pass: 1.1 of the fluid limit for pure ice, the lightest material of `rings` (so the largest limit)."""
    density = body_mass_kg / (4.0 / 3.0 * math.pi * (body_radius_km * 1000.0) ** 3)
    return rings.roche_limit_km(body_radius_km, density, rings.ICE_DENSITY_KG_M3) * rings.RANGES["outer_over_roche"][1] * 1.0001


def _rogues(count=150):
    cfg = SystemConfig()
    for seed in range(count):
        with draw.bound(seed):
            yield RoguePlanet(cfg)


def test_a_rogue_planet_has_no_ring_because_no_moon_or_comet_history_is_followed_for_it():
    for rogue in _rogues():
        assert all(getattr(rogue, name) is None for name in rings.BODY_FIELDS)


def test_the_same_seed_gives_the_same_rogue_ring():
    def ring(seed):
        with draw.bound(seed):
            rogue = RoguePlanet(SystemConfig())
        return tuple(getattr(rogue, name) for name in rings.BODY_FIELDS)
    assert [ring(seed) for seed in range(20)] == [ring(seed) for seed in range(20)]


def test_a_moon_s_ring_obeys_the_same_law_and_a_moon_has_no_moon_to_stop_it():
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = False
    ringed = moons = 0
    for seed in range(60):
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        for planet in system.planets:
            for moon in getattr(planet, "moons", ()):
                moons += 1
                if moon.ring_kind is None:
                    continue
                ringed += 1
                assert moon.radius < moon.ring_inner_km < moon.ring_outer_km <= _roche_ceiling_km(moon.radius, moon.mass)
                assert moon.ring_mass_kg <= rings.MAX_MASS_OVER_PLANET * moon.mass
                assert moon.born_beyond_snow_line or moon.ring_ice_fraction == 0.0
    assert moons > 50 and ringed < 0.1 * moons    # a moon's Roche sphere is a tiny target for a comet, and it has no moon to lose
