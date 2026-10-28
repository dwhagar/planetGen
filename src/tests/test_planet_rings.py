# tests/test_planet_rings.py

"""
GEN.732: planetary rings. The Roche limit and the optical depth are checked against the four giants of the Solar System; a draw keeps
every ring inside the Roche limit of its own material and inside the innermost moon; a giant has a ring; a moon, a belt and a
reclassified body are handled.
"""

import math

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import moon_composition, rings
from planetgen.util import draw

# name: (radius km, mass kg, material density kg/m3, outermost edge of the observed rings km, mean optical-depth case)
GIANTS = {
    "Saturn": (58232.0, 5.683e26, 917.0, 136775.0),
    "Jupiter": (71492.0, 1.898e27, 2500.0, 129000.0),
    "Uranus": (25362.0, 8.681e25, 1000.0, 51149.0),
    "Neptune": (24622.0, 1.024e26, 1000.0, 62933.0),
}


def _density(radius_km, mass_kg):
    return mass_kg / (4.0 / 3.0 * math.pi * (radius_km * 1000.0) ** 3)


@pytest.mark.parametrize("name", sorted(GIANTS))
def test_the_four_giants_rings_end_within_a_few_percent_below_the_roche_limit_of_their_material(name):
    radius, mass, material, edge = GIANTS[name]
    limit = rings.roche_limit_km(radius, _density(radius, mass), material)
    low, high = rings.RANGES["outer_over_roche"]
    assert low <= edge / limit <= high


def test_saturns_rings_have_the_optical_depth_of_a_dense_ring_and_uranus_s_that_of_narrow_ones():
    saturn = rings.optical_depth(1.5e19, 74500.0, 140220.0, 1.0, 0.95, rings.RANGES["particle_radius_m"]["broad"])
    uranus = rings.optical_depth(1.5e16, 41837.0, 51149.0, 3.0e-3, 0.0, rings.RANGES["particle_radius_m"]["narrow"])
    assert 0.1 < saturn < 1.0
    assert 0.1 < uranus < 3.0


def test_the_mixture_density_is_volume_additive_between_ice_and_rock():
    assert rings.mixture_density_kg_m3(0.0) == pytest.approx(rings.ROCK_DENSITY_KG_M3)
    assert rings.mixture_density_kg_m3(1.0) == pytest.approx(rings.ICE_DENSITY_KG_M3)
    assert rings.ROCK_DENSITY_KG_M3 > rings.mixture_density_kg_m3(0.5) > rings.ICE_DENSITY_KG_M3


def _draws(parent, radius, mass, cold, moons=(), n=400, seed=0):
    out = []
    with draw.bound(draw.Stream(seed)):
        for _ in range(n):
            out.append(rings.draw_ring(parent, radius, mass, cold, moons))
    return out


def test_a_giant_s_rings_obey_the_physics_the_module_names():
    radius, mass = 60000.0, 5.7e26
    density = _density(radius, mass)
    drawn = [ring for ring in _draws(moon_composition.GIANT, radius, mass, True) if ring is not None]
    assert len(drawn) > 350
    assert {ring.ring_kind for ring in drawn} == set(rings.KINDS)
    for ring in drawn:
        low, high = rings.RANGES["mass_kg"][ring.ring_kind]
        assert low <= ring.ring_mass_kg <= high
        limit = rings.roche_limit_km(radius, density, rings.mixture_density_kg_m3(ring.ring_ice_fraction))
        assert ring.ring_outer_km <= limit * rings.RANGES["outer_over_roche"][1] * (1 + 1e-9)
        assert radius * rings.RANGES["inner_over_radius"][0] <= ring.ring_inner_km < ring.ring_outer_km
        assert ring.ring_optical_depth > 0.0
        assert ring.ring_origin in rings.ORIGINS


def test_a_ring_beyond_the_snow_line_keeps_its_ice_and_one_inside_has_none():
    cold = [r for r in _draws(moon_composition.GIANT, 60000.0, 5.7e26, True) if r and r.ring_kind == "broad"]
    warm = [r for r in _draws(moon_composition.GIANT, 60000.0, 5.7e26, False) if r and r.ring_kind == "broad"]
    assert cold and warm
    assert all(0.90 <= r.ring_ice_fraction <= 0.99 for r in cold)
    assert all(r.ring_ice_fraction == 0.0 for r in warm)
    assert all(r.ring_kind != "dusty" or r.ring_ice_fraction == 0.0 for r in _draws(moon_composition.GIANT, 60000.0, 5.7e26, True))


def test_a_ring_ends_inside_the_innermost_moon():
    inner_moon = 120000.0
    for ring in _draws(moon_composition.GIANT, 60000.0, 5.7e26, True, moons=[inner_moon, 900000.0]):
        assert ring is None or ring.ring_outer_km <= rings.MOON_CLEARANCE * inner_moon * (1 + 1e-9)


def test_a_moon_too_close_leaves_no_room_for_a_ring():
    assert all(ring is None for ring in _draws(moon_composition.GIANT, 60000.0, 5.7e26, True, moons=[61000.0], n=50))


def test_a_terrestrial_planet_seldom_has_a_ring_and_the_same_seed_gives_the_same_rings():
    share = sum(r is not None for r in _draws(moon_composition.TERRESTRIAL, 6371.0, 5.97e24, True, n=3000)) / 3000
    assert 0.0 < share < 0.03
    first = [r.values() if r else None for r in _draws(moon_composition.GIANT, 60000.0, 5.7e26, True, n=20, seed=7)]
    again = [r.values() if r else None for r in _draws(moon_composition.GIANT, 60000.0, 5.7e26, True, n=20, seed=7)]
    assert first == again


def _systems():
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = False
    for seed in range(40):
        with draw.bound(seed):
            yield StarSystem(system_config=cfg)


def test_generated_giants_have_rings_and_moons_and_belts_have_none():
    giants = ringed = 0
    for system in _systems():
        for planet in system.all_planets:
            if planet.body_type == "a":
                assert getattr(planet, "ring_kind", None) is None
                continue
            for moon in planet.moons:
                assert all(getattr(moon, name) is None for name in rings.BODY_FIELDS)
            if planet.body_type == "g":
                giants += 1
                ringed += planet.ring_kind is not None
                if planet.ring_kind is not None:
                    assert planet.ring_inner_km < planet.ring_outer_km
    assert giants > 10 and ringed >= 0.8 * giants


def test_a_body_regenerated_as_another_class_loses_its_ring():
    for system in _systems():
        for planet in system.planets:
            if getattr(planet, "ring_kind", None) is not None:
                rings.clear(planet)
                assert all(getattr(planet, name) is None for name in rings.BODY_FIELDS)
                return
    pytest.fail("no ringed planet in 40 systems")
