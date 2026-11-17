# tests/test_dust_supply.py

"""GEN.1542 slice 0: the carbon a planet's cosmic dust supplies the V8 layer (`physics/dust_supply.py`), from the influx of GEN.1543."""

import math

import pytest

from planetgen.physics import dust_supply, influx

MAIN_BELT = influx.Reservoir("belt", 2.4e21, 2.2, 3.2, 4.6e9)     # a test input: a Sun with a main-belt sized reservoir
EARTH, VENUS, MARS = (6371.0, 5.972e24, 1.0), (6052.0, 4.867e24, 0.723), (3390.0, 6.4e23, 1.524)


def _flux(body, reservoirs=(MAIN_BELT,)):
    radius, mass, au = body
    influx_kg_yr = influx.planet_influx_kg_per_yr(radius, mass, au, 1.0, 1.0, list(reservoirs), [])[0]
    return dust_supply.carbon_flux_kg_m2_yr(influx_kg_yr, radius)


def test_the_carbon_is_the_dust_times_its_carbon_share_over_the_planets_area():
    area = 4.0 * math.pi * (6371.0e3) ** 2
    assert dust_supply.carbon_flux_kg_m2_yr(1.0e7, 6371.0) == pytest.approx(dust_supply.CARBON_SHARE_OF_DUST * 1.0e7 / area)
    assert dust_supply.carbon_flux_kg_m2_yr(2.0e7, 6371.0) == pytest.approx(2.0 * dust_supply.carbon_flux_kg_m2_yr(1.0e7, 6371.0))
    assert dust_supply.carbon_flux_kg_m2_yr(1.0e7, 12742.0) == pytest.approx(dust_supply.carbon_flux_kg_m2_yr(1.0e7, 6371.0) / 4.0), "the same dust over twice the radius is a quarter of the flux"


def test_no_dust_or_no_radius_gives_no_carbon():
    assert dust_supply.carbon_flux_kg_m2_yr(0.0, 6371.0) == 0.0
    assert dust_supply.carbon_flux_kg_m2_yr(1.0e7, 0.0) == 0.0
    assert dust_supply.carbon_flux_kg_m2_yr(-1.0, 6371.0) == 0.0
    assert not dust_supply.v8_supplied(0.0, 6371.0)


def test_the_inner_planets_of_a_belted_system_hold_the_floor_and_a_far_small_one_does_not():
    # the floor sits between a Mars-like body (about 9e-10) and an Earth- or Venus-like one (about 3e-9) of a Sun with a main-belt sized reservoir
    assert _flux(EARTH) >= dust_supply.V8_MIN_CARBON_FLUX_KG_M2_YR
    assert _flux(VENUS) >= dust_supply.V8_MIN_CARBON_FLUX_KG_M2_YR
    assert _flux(MARS) < dust_supply.V8_MIN_CARBON_FLUX_KG_M2_YR


def test_a_system_with_no_reservoir_gives_no_supply():
    assert _flux(VENUS, reservoirs=()) == 0.0


def test_the_gate_follows_the_dust_supplied():
    radius = 6052.0
    area = 4.0 * math.pi * (radius * 1000.0) ** 2
    needed = dust_supply.V8_MIN_CARBON_FLUX_KG_M2_YR * area / dust_supply.CARBON_SHARE_OF_DUST
    assert dust_supply.v8_supplied(needed * 1.01, radius) and not dust_supply.v8_supplied(needed * 0.99, radius)
