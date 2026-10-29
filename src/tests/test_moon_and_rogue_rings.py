# tests/test_moon_and_rogue_rings.py

"""
GEN.836: a moon and a rogue planet draw a ring system by the law of GEN.732 (`planetgen.physics.rings`). Across the mass range every ring
ends inside the Roche limit of its material, a brown dwarf has none, and a rogue giant always has one.
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


def test_a_rogue_giant_always_has_a_ring_inside_its_roche_limit():
    giants = 0
    for rogue in _rogues():
        if rogue.planet_class is None:
            assert rogue.ring_kind is None  # a brown dwarf is no planet
            continue
        if rogue.planet_type == "g":
            giants += 1
            assert rogue.ring_kind is not None
        if rogue.ring_kind is not None:
            assert rogue.ring_kind in rings.KINDS and rogue.ring_origin in rings.ORIGINS
            assert rogue.radius_km < rogue.ring_inner_km < rogue.ring_outer_km <= _roche_ceiling_km(rogue.radius_km, rogue.mass_kg)
            assert 0.0 < rogue.ring_mass_kg <= rings.MAX_MASS_OVER_PLANET * rogue.mass_kg
            assert rogue.ring_optical_depth > 0.0 and 0.0 <= rogue.ring_ice_fraction <= 1.0
    assert giants > 5


def test_a_cold_rogue_ring_keeps_its_ice_and_a_warm_one_has_none():
    seen = set()
    for rogue in _rogues(300):
        if rogue.ring_kind is None:
            continue
        cold = rogue.effective_temperature_k < constants.SNOW_LINE_TEMPERATURE_K
        assert cold or rogue.ring_ice_fraction == 0.0
        seen.add(cold)
    assert True in seen


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
    assert moons > 50 and ringed < 0.1 * moons
