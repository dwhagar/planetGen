# tests/test_cosmic_dust.py

"""GEN.1542: the cosmic dust that reaches a body, and the V8 gate it sets."""

import math

import pytest

from planetgen.physics import constants, cosmic_dust as dust, ring_origin, sulfur_cloud_set as sc

SUN = constants.SOLAR_MASS_TO_KG


def tonnes_a_day(mass_kg, radius_m, distance_au, star_mass_kg=SUN):
    return dust.influx_kg_yr(mass_kg, radius_m, distance_au, star_mass_kg) / 1000.0 / 365.25


def test_the_earth_gets_the_measured_28_tonnes_a_day():
    earth = dust.EARTH
    assert tonnes_a_day(earth["mass_kg"], earth["radius_m"], earth["distance_au"]) == pytest.approx(28.0, rel=0.03)


def test_the_influx_falls_with_the_distance_and_rises_with_size_and_the_star_mass():
    base = dust.influx_kg_yr(5.0e24, 6.0e6, 1.0, SUN)
    assert dust.influx_kg_yr(5.0e24, 6.0e6, 2.0, SUN) < base
    light = dust.influx_kg_yr(1.0e10, 6.0e6, 1.0, SUN)         # no gravity worth the name: no focusing, so the law alone
    assert dust.influx_kg_yr(1.0e10, 6.0e6, 2.0, SUN) == pytest.approx(light * 2.0 ** -dust.DISTANCE_EXPONENT, rel=1e-3)
    assert dust.influx_kg_yr(5.0e24, 7.0e6, 1.0, SUN) > base
    assert dust.influx_kg_yr(1.0e10, 6.0e6, 1.0, 2.0 * SUN) == pytest.approx(light * 2.0 ** dust.STAR_MASS_EXPONENT, rel=1e-3)


def test_a_body_with_no_size_or_place_gets_none():
    assert dust.influx_kg_yr(0.0, 6.0e6, 1.0, SUN) == 0.0
    assert dust.influx_kg_yr(5.0e24, 0.0, 1.0, SUN) == 0.0
    assert dust.influx_kg_yr(5.0e24, 6.0e6, None, SUN) == 0.0
    assert dust.influx_per_area_kg_m2_yr(5.0e24, 0.0, 1.0, SUN) == 0.0


def test_the_influx_per_area_is_the_influx_over_the_whole_surface():
    mass, radius = 5.0e24, 6.0e6
    assert dust.influx_per_area_kg_m2_yr(mass, radius, 1.0, SUN) == pytest.approx(dust.influx_kg_yr(mass, radius, 1.0, SUN) / (4.0 * math.pi * radius ** 2))


def test_venus_is_the_reference_of_the_supply_floor():
    assert dust.share_of_reference(**dust.VENUS) == pytest.approx(1.0)
    venus = dust.VENUS
    assert sc.dust_supply(venus["mass_kg"], venus["radius_m"], venus["distance_au"], venus["star_mass_kg"])


def test_the_floor_passes_the_inner_planets_and_closes_the_outer_belt():
    mars = dict(mass_kg=6.417e23, radius_m=3.3895e6, distance_au=1.524, star_mass_kg=SUN)
    belt_body = dict(mass_kg=1.0e21, radius_m=4.5e5, distance_au=2.8, star_mass_kg=SUN)
    cold = dict(mass_kg=1.0e22, radius_m=1.0e6, distance_au=8.0, star_mass_kg=SUN)
    assert sc.dust_supply(**mars)
    assert not sc.dust_supply(**belt_body)
    assert not sc.dust_supply(**cold)
    assert sc.number("dust_min_share_of_reference") == 0.1


def test_the_focusing_the_ring_law_uses_is_the_one_the_dust_uses():
    planet = ring_origin.Planet(mass_kg=5.0e24, radius_km=6000.0, distance_au=1.0, star_mass_kg=SUN)
    inner, outer = ring_origin.annulus_m(planet, True)
    focus = ring_origin.focusing(planet.mass_kg, 0.5 * (inner + outer), planet.distance_au, planet.star_mass_kg)
    v_inf = ring_origin.COMET_ENCOUNTER_FRACTION * math.sqrt(constants.G * SUN / constants.AU_TO_M)
    assert focus == pytest.approx(1.0 + 2.0 * constants.G * planet.mass_kg / (0.5 * (inner + outer)) / v_inf ** 2)
