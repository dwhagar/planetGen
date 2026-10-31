# tests/test_planet_rings.py

"""
GEN.732, GEN.1057: planetary rings. The Roche limit and the optical depth are checked against the four giants of the Solar System; a built ring keeps
every ring inside the Roche limit of its own material and inside the innermost moon; a ring is the result of history (a fallen moon, a passing comet); a belt and a
reclassified body are handled (moons and rogue planets: tests/test_moon_and_rogue_rings.py).
"""

import math
import types

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import ring_origin, rings
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


MASSES_KG = (1.0e8, 1.0e10, 1.0e15, 1.0e17, 1.5e19, 1.0e21)
"""tuple: Ring masses across the three kinds, for the geometry tests."""


def _builds(radius, mass, cold, moons=(), n=100, seed=0):
    out = []
    with draw.bound(draw.Stream(seed)):
        for _ in range(n):
            for ring_mass in MASSES_KG:
                out.append(rings.build_ring("disrupted_moon", ring_mass, 2.0e8, radius, mass, cold, moons))
    return out


def test_the_kind_of_a_ring_follows_from_its_mass():
    assert [rings.kind_of_mass(mass) for mass in MASSES_KG] == ["dusty", "dusty", "narrow", "narrow", "broad", "broad"]
    assert rings.kind_of_mass(1.5e19) == "broad" and rings.kind_of_mass(1.0e16) == "narrow" and rings.kind_of_mass(1.0e9) == "dusty"


def test_a_ring_of_any_kind_obeys_the_physics_the_module_names():
    radius, mass = 60000.0, 5.7e26
    density = _density(radius, mass)
    drawn = [ring for ring in _builds(radius, mass, True) if ring is not None]
    assert len(drawn) > 500
    assert {ring.ring_kind for ring in drawn} == set(rings.KINDS)
    for ring in drawn:
        assert ring.ring_kind == rings.kind_of_mass(ring.ring_mass_kg) and ring.ring_age_yr == 2.0e8
        limit = rings.roche_limit_km(radius, density, rings.mixture_density_kg_m3(ring.ring_ice_fraction))
        assert ring.ring_outer_km <= limit * rings.RANGES["outer_over_roche"][1] * (1 + 1e-9)
        assert radius * rings.RANGES["inner_over_radius"][0] <= ring.ring_inner_km < ring.ring_outer_km
        assert ring.ring_optical_depth > 0.0
        assert ring.ring_origin in rings.ORIGINS


def test_a_ring_beyond_the_snow_line_keeps_its_ice_and_one_inside_has_none():
    cold = [r for r in _builds(60000.0, 5.7e26, True) if r and r.ring_kind == "broad"]
    warm = [r for r in _builds(60000.0, 5.7e26, False) if r and r.ring_kind == "broad"]
    assert cold and warm
    assert all(0.90 <= r.ring_ice_fraction <= 0.99 for r in cold)
    assert all(r.ring_ice_fraction == 0.0 for r in warm)
    assert all(r.ring_kind != "dusty" or r.ring_ice_fraction == 0.0 for r in _builds(60000.0, 5.7e26, True) if r)


def test_a_ring_ends_inside_the_innermost_moon():
    inner_moon = 120000.0
    for ring in _builds(60000.0, 5.7e26, True, moons=[inner_moon, 900000.0]):
        assert ring is None or ring.ring_outer_km <= rings.MOON_CLEARANCE * inner_moon * (1 + 1e-9)


def test_a_moon_too_close_leaves_no_room_for_a_ring():
    assert all(ring is None for ring in _builds(60000.0, 5.7e26, True, moons=[61000.0], n=20))


def test_the_same_seed_gives_the_same_rings():
    first = [r.values() if r else None for r in _builds(60000.0, 5.7e26, True, n=20, seed=7)]
    again = [r.values() if r else None for r in _builds(60000.0, 5.7e26, True, n=20, seed=7)]
    assert first == again


def _systems(count=40):
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = False
    for seed in range(count):
        with draw.bound(seed):
            yield StarSystem(system_config=cfg)


def test_a_ring_is_a_result_of_history_and_belts_have_none():
    giants = ringed = 0
    for system in _systems(count=200):
        if giants >= 10:      # giants are rarer since the ice boost fell to 2.17 (GEN.950): generate systems until there are enough to count
            break
        for planet in system.all_planets:
            if planet.body_type == "a":
                assert getattr(planet, "ring_kind", None) is None
                continue
            for moon in planet.moons:
                assert moon.ring_kind is None or moon.ring_inner_km < moon.ring_outer_km  # GEN.836: a moon may hold one
            if planet.body_type == "g":
                giants += 1
                ringed += planet.ring_kind is not None
                if planet.ring_kind is not None:
                    assert planet.ring_inner_km < planet.ring_outer_km
    assert giants >= 10 and 0 < ringed < giants    # GEN.1057: a giant with neither a fallen moon nor a passing comet has no ring


def test_a_body_regenerated_as_another_class_loses_its_ring():
    for system in _systems():
        for planet in system.planets:
            if getattr(planet, "ring_kind", None) is not None:
                rings.clear(planet)
                assert all(getattr(planet, name) is None for name in rings.BODY_FIELDS)
                return
    pytest.fail("no ringed planet in 40 systems")


def _planet(mass=5.97e24, radius=6371.0, distance=1.0, fallen=(), moons=(), age_gy=4.5, body_type="t", spin_h=10.0):
    star = types.SimpleNamespace(snow_line_au=3.0, mass=1.989e30, age=age_gy)
    return types.SimpleNamespace(body_type=body_type, radius=radius, mass=mass, envelope_mass_fraction=0.0, is_moon=False, born_beyond_snow_line=distance >= 3.0,
                                 distance=distance, star=star, moons=list(moons), fallen_moons=tuple(fallen), rotation_period_hours=spin_h,
                                 **{name: None for name in rings.BODY_FIELDS})


def test_a_planet_with_no_fallen_moon_and_no_comet_has_no_ring():
    with draw.bound(draw.Stream(3)):
        rings.apply(planet := _planet(distance=0.3, age_gy=0.001))
    assert all(getattr(planet, name) is None for name in rings.BODY_FIELDS)


def test_a_fallen_moon_leaves_a_disrupted_moon_ring_with_its_age():
    moon = ring_origin.FallenMoon(2.0e20, 4.0e9)
    with draw.bound(draw.Stream(3)):
        rings.apply(planet := _planet(distance=1.0, fallen=(moon,), age_gy=4.5))
    assert planet.ring_kind is not None and planet.ring_origin == "disrupted_moon"
    assert planet.ring_age_yr == pytest.approx(5.0e8)
    assert planet.ring_mass_kg == pytest.approx(ring_origin.spread(ring_origin.RING_RETAINED_FRACTION * 2.0e20, 5.0e8), rel=0.05) or planet.ring_origin == "captured_body"


def test_a_planet_inside_the_inner_limit_has_no_ring_even_after_a_moon_fell():
    with draw.bound(draw.Stream(3)):
        rings.apply(planet := _planet(distance=0.05, fallen=(ring_origin.FallenMoon(2.0e20, 4.0e9),)))
    assert planet.ring_kind is None


def test_a_dwarf_planet_may_have_a_narrow_resonance_ring_without_an_age():
    found = []
    for seed in range(400):
        with draw.bound(draw.Stream(seed)):
            rings.apply(planet := _planet(mass=1.0e21, radius=500.0, distance=40.0, body_type="d", spin_h=4.0))
        if planet.ring_kind is not None:
            found.append(planet)
    assert found and all(p.ring_kind == "narrow" and p.ring_origin == "disrupted_moon" and p.ring_age_yr is None for p in found)
    assert len(found) < 100


def _ring(kind, ice, tau, age_yr, inner_km, outer_km):
    return rings.Ring(ring_kind=kind, ring_origin="disrupted_moon", ring_inner_km=inner_km, ring_outer_km=outer_km, ring_mass_kg=1.0e19, ring_optical_depth=tau,
                      ring_ice_fraction=ice, ring_age_yr=age_yr)


def test_a_massive_ring_beyond_the_snow_line_keeps_most_of_its_ice_as_saturns_does():
    """GEN.1055 (R6): Saturn's A to C rings (ice began at 0.97, optical depth 0.5, 4.5 Gyr, 9.58 AU) are still more than 90 percent ice after the dust that fell on them."""
    ring = _ring("broad", 0.97, 0.5, 4.5e9, 74500.0, 136775.0)
    rings.darken(ring, 5.683e26, 9.58)
    assert 0.9 < ring.ring_ice_fraction < 0.97


def test_a_thin_ring_is_darkened_much_further_than_a_massive_one_and_an_older_one_further_than_a_young_one():
    thin = _ring("narrow", 0.3, 1.0e-3, 4.5e9, 40000.0, 50000.0)
    heavy = _ring("narrow", 0.3, 1.0, 4.5e9, 40000.0, 50000.0)
    young = _ring("narrow", 0.3, 1.0e-3, 1.0e6, 40000.0, 50000.0)
    for ring in (thin, heavy, young):
        rings.darken(ring, 8.68e25, 19.2)
    assert 0.3 - thin.ring_ice_fraction > 50.0 * (0.3 - heavy.ring_ice_fraction)
    assert young.ring_ice_fraction > thin.ring_ice_fraction
    assert 0.0 <= thin.ring_ice_fraction < 0.3


def test_a_ring_with_no_ice_or_no_age_is_left_as_it_is():
    dusty = _ring("dusty", 0.0, 1.0e-8, 1.0e8, 80000.0, 130000.0)
    resonance = _ring("narrow", 0.2, 0.1, None, 4000.0, 4100.0)
    for ring in (dusty, resonance):
        before = ring.ring_ice_fraction
        rings.darken(ring, 1.0e21, 40.0)
        assert ring.ring_ice_fraction == before
