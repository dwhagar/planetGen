# tests/test_dust_supply.py

"""GEN.1543: the cosmic dust a rocky body takes from its own system's belts and comets, and the V8 gate that reads it."""

import math

import pytest

from planetgen.generation import dust_supply
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants, influx, sulfur_cloud_set as sc
from planetgen.util import draw

VENUS_RADIUS_M = 6.0518e6


def test_the_gate_is_a_floor_on_the_influx_per_area():
    area = 4.0 * math.pi * VENUS_RADIUS_M ** 2
    floor = sc.number("dust_min_influx_kg_m2_yr")
    assert sc.dust_supply(1.01 * floor * area, VENUS_RADIUS_M)
    assert not sc.dust_supply(0.99 * floor * area, VENUS_RADIUS_M)
    assert not sc.dust_supply(None, VENUS_RADIUS_M) and not sc.dust_supply(0.0, VENUS_RADIUS_M) and not sc.dust_supply(1.0e9, 0.0)


def test_the_floor_is_a_tenth_of_venus_in_the_model_with_a_solar_system_asteroid_belt():
    belt = influx.Reservoir("belt", 2.4e21, 2.2, 3.2, 4.6e9)
    venus = influx.planet_influx_kg_per_yr(6051.8, 0.72333, 1.0, [belt], [], 1.0)[0]
    per_area = venus / (4.0 * math.pi * VENUS_RADIUS_M ** 2)
    assert sc.number("dust_min_influx_kg_m2_yr") == pytest.approx(0.1 * per_area, rel=0.05)
    assert sc.dust_supply(venus, VENUS_RADIUS_M)


def test_reservoirs_are_the_belts_with_a_ring_and_the_comets_with_a_nucleus():
    class Belt:
        mass_kg, lower_limit, upper_limit = 1.0e21, 2.0, 3.0

    class NoMass:
        mass_kg, lower_limit, upper_limit = None, 2.0, 3.0

    class Comet:
        mass_kg, perihelion_distance_au = 1.0e14, 0.7

    found = dust_supply.reservoirs_of([Belt(), NoMass()], [Comet()], 4.6)
    assert [r.kind for r in found] == ["belt", "comet"]
    assert found[0].age_yr == pytest.approx(4.6e9) and found[1].lower_au == 0.7


def _generated(seed):
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    config.PLANETS = True
    config.BINARY_SYSTEM = False
    draw.set_run_seed(seed)
    return StarSystem(system_config=config)


def test_a_generated_system_gives_each_rocky_body_its_dust_and_a_giant_none():
    for seed in range(1, 40):
        system = _generated(seed)
        if any(b.body_type == "a" for b in system.planets) and any(b.body_type == "t" for b in system.planets):
            break
    else:
        pytest.fail("no system with a belt and a rocky planet in 39 seeds")
    rocky = [b for p in system.all_planets if p.body_type not in ("a", "g") for b in (p, *p.moons) if b.body_type not in ("a", "g")]
    assert rocky and all(getattr(b, dust_supply.ATTRIBUTE, None) is not None for b in rocky)
    assert all(getattr(b, dust_supply.ATTRIBUTE) >= 0.0 for b in rocky)
    assert any(getattr(b, dust_supply.ATTRIBUTE) > 0.0 for b in rocky), "a belt sends something to a planet inside it"
    assert all(not hasattr(b, dust_supply.ATTRIBUTE) for p in system.all_planets if p.body_type in ("a", "g") for b in [p])


def test_a_system_without_belts_or_comets_sends_no_dust():
    class Star:
        mass = constants.SOLAR_MASS_TO_KG
        age = 4.6
        luminosity = constants.SOLAR_LUMINOSITY

    class Planet:
        star = Star()
        radius, distance, body_type = 6371.0, 1.0, "t"

    assert dust_supply.influx_of(Planet(), Planet(), [], []) == 0.0
