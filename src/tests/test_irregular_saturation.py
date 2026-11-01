# tests/test_irregular_saturation.py

"""
GEN.1079: the collisional-saturation check of the irregular zones (`planetgen.physics.irregular_population`; docs/design/systems/captured-asteroid-moons.md section 2.3 and item 5 of section 7).

The design expected the four giants to end with about the same population whatever they captured. Two ways to check it were named: reuse `belt_mass.belt_mass_earth` for the zone, or fit a scalar per planet.
The result (report gen1079-irregular-saturation-2026-10-30.md): the repo's belt law does not describe an irregular zone (it leaves 1e15 to 1e17 kg, 400 to 4400 times under what the giants hold, and varies a factor 125 across them),
and a fitted scalar fixes the counts above 5 km only by breaking the counts above the storage floor, which are the ones that are stored. The module's own grinding (intrinsic collision probability of the zone) is partly saturated and
keeps the four giants within a factor 2.5 of each other in counts and mass, so these tests pin that, with no factor added.
"""

import math

import pytest

from planetgen.physics import belt_mass, constants, irregular_population as ip, irregular_zone

SUN = constants.SOLAR_MASS_TO_KG
EARTH = constants.EARTH_MASS_TO_KG
AGE_YR = 4.5e9
FLOOR_RADIUS_KM = ip.radius_km(ip.STORE_FLOOR_KG)
GIANTS = {"Jupiter": (5.203, 1.898e27, True), "Saturn": (9.537, 5.683e26, True), "Uranus": (19.19, 8.681e25, False), "Neptune": (30.07, 1.024e26, False)}
FLOOR_COUNTS = {"Jupiter": 8, "Saturn": 5, "Uranus": 6, "Neptune": 8}         # CM6: irregulars above 1e16 kg
OBSERVED_MASS_KG = {"Jupiter": 3.9e18, "Saturn": 1.0e19, "Uranus": 7.5e18, "Neptune": 4.8e19}


def _giant(name, encounters=ip.REFERENCE_ENCOUNTERS):
    distance, mass, gas = GIANTS[name]
    hill = irregular_zone.hill_radius_au(distance, mass, SUN)
    pool = ip.make_pool(hill, encounters, ip.REFERENCE_DISC_MASS_EARTH, gas, 3.0 if gas else 0.0)
    return pool, ip.intrinsic_probability(mass, ip.zone_radius_au(hill, distance))


def _count(name, radius_km, encounters=ip.REFERENCE_ENCOUNTERS):
    pool, probability = _giant(name, encounters)
    return ip.ground_count_above(pool, radius_km, probability, AGE_YR)


def _mass(name):
    pool, probability = _giant(name)
    return ip.mass_below_kg(pool, probability, AGE_YR, ip.MODEL_MAX_RADIUS_KM)


def _spread(values):
    return max(values) / min(values)


def test_the_four_giants_end_within_a_factor_of_two_and_a_half_though_they_capture_differently():
    for radius in (5.0, FLOOR_RADIUS_KM, 40.0):
        counts = [_count(name, radius) for name in GIANTS]
        assert _spread(counts) < 2.5
    captured = [_giant(name)[0].count_above(5.0) for name in GIANTS]
    assert _spread(captured) > 1.0                           # the gas giants capture more (the gas window), the ice giants do not
    assert _spread([_mass(name) for name in GIANTS]) < 3.0


def test_the_ground_count_follows_what_was_captured_by_less_than_its_square_root_for_the_planets_that_grind():
    low, high = 49, 280                                      # the encounters of the three histories read (CM1)
    elasticity = {}
    for name in GIANTS:
        elasticity[name] = math.log(_count(name, 5.0, high) / _count(name, 5.0, low)) / math.log(high / low)
    assert all(0.2 < e < 0.9 for e in elasticity.values())   # 0 would be full saturation, 1 none
    assert max(elasticity, key=elasticity.get) == "Neptune"  # the big, slow zone grinds least
    assert sum(elasticity[name] for name in ("Jupiter", "Saturn", "Uranus")) / 3.0 < 0.5


def test_the_counts_above_the_storage_floor_are_the_observed_ones_within_a_factor_of_one_and_a_half():
    for name, observed in FLOOR_COUNTS.items():
        assert observed / 1.6 < _count(name, FLOOR_RADIUS_KM) < observed * 1.6


def test_a_single_scalar_on_the_collision_probability_cannot_fix_the_counts_above_5_km_without_breaking_the_floor_counts():
    observed_5km = {"Jupiter": 8, "Saturn": 11, "Uranus": 9}      # Neptune's survey is complete only above 16 km

    def ratios(scale):
        high = [ip.ground_count_above(_giant(n)[0], 5.0, _giant(n)[1] * scale, AGE_YR) / observed_5km[n] for n in observed_5km]
        floor = [ip.ground_count_above(_giant(n)[0], FLOOR_RADIUS_KM, _giant(n)[1] * scale, AGE_YR) / FLOOR_COUNTS[n] for n in observed_5km]
        return high, floor

    high, floor = ratios(1.0)                                     # as built: the 5 km counts are 2 to 3 times the observed, the floor counts within 1.5
    assert max(high) > 2.5 and all(1.0 / 1.6 < r < 1.6 for r in floor)
    high, floor = ratios(5.0)                                     # the scalar that fits Jupiter's 5 km count: the floor counts fall to a fifth
    assert high[0] == pytest.approx(1.0, abs=0.1) and min(floor) < 0.25


def test_the_repos_belt_law_does_not_describe_an_irregular_zone():
    """`belt_mass_earth` [R4785] for a zone round a planet, the planet as the central mass: far under the mass the giants hold, and a factor 100 across them."""
    left = {}
    for name, (distance, mass, _) in GIANTS.items():
        zone = irregular_zone.zone(distance, mass, SUN)
        lower, upper = 0.1 * zone.hill_radius_au, zone.retrograde_limit_au
        a = belt_mass.collisional_constant(lower, upper, mass / SUN, 100.0)
        left[name] = a / (AGE_YR / 1.0e6) * EARTH                  # the late-time limit M_max, kg
    assert all(left[name] < OBSERVED_MASS_KG[name] / 40.0 for name in GIANTS)
    assert _spread(list(left.values())) > 50.0
    assert _spread([_mass(name) for name in GIANTS]) < 3.0         # the module's own grinding is the near-invariant one
