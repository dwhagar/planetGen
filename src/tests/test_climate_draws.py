# tests/test_climate_draws.py

"""
GEN.412: a planet or moon keeps the albedo and the greenhouse quantity it drew, and its surface temperature is exactly
what they give (`physics/climate_draws.py`); a stored row is inverted back to them with no draw, or takes a seeded class draw
where the inverse has no solution.
"""

import math

import pytest

from planetgen import tuning
from planetgen.generation import plausibility
from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.star import Star
from planetgen.physics import climate_draws, constants, orbit_elements, surface_air
from planetgen.util import draw

TOLERANCE = 1e-9
"""float: The relative error the inverse may have (it is exact up to floating-point rounding)."""


@pytest.fixture(scope="module")
def host_star():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    return Star(cfg)


def _flux(planet, star):
    # GEN.727: the light of an eccentric orbit, averaged over it, is that of a circle of radius a (1 - e^2)^(1/4)
    distance_km = orbit_elements.flux_distance_au(float(planet.distance), planet.eccentricity) * constants.AU_TO_KM
    return climate_draws.starlight_w_m2(star.luminosity, distance_km)


def _stored_row(planet):
    """The columns a saved body would hold, by their column names."""
    return {
        "planet_class": planet.planet_class, "body_type": planet.body_type, "atmosphere": planet.atmosphere,
        "surface_temperature_k": planet.surface_temperature, "scale_height_km": planet.scale_height,
        "atm_molar_density": planet.atm_molar_density, "gravity_g": planet.gravity, "atm_density": planet.atm_density,
        "atmospheric_pressure_pa": planet.atmospheric_pressure,
    }


def test_the_formulas_are_inverses():
    flux = 1361.0
    for albedo in (0.0, 0.12, 0.3, 0.8):
        for factor in (0.0, 0.5, 40.0, 100.0):
            temperature = climate_draws.surface_temperature_k(albedo, factor, flux)
            assert temperature >= climate_draws.equilibrium_temperature_k(albedo, flux)
            assert climate_draws.greenhouse_from_temperature(temperature, albedo, flux) == pytest.approx(
                factor, rel=1e-9, abs=1e-9)


def test_no_air_means_no_greenhouse():
    assert climate_draws.surface_temperature_k(0.3, 0.0, 1361.0) == climate_draws.equilibrium_temperature_k(0.3, 1361.0)


def test_nothing_is_colder_than_the_background():
    assert climate_draws.equilibrium_temperature_k(0.3, 1e-30) == constants.COSMIC_BACKGROUND_TEMPERATURE_K
    assert climate_draws.surface_temperature_k(0.3, 5.0, 1e-30) == constants.COSMIC_BACKGROUND_TEMPERATURE_K


def _temperature_of_draws(planet, flux):
    """What the stored albedo and greenhouse quantity give: the starlight alone, `T_e (1 + G)^(1/4)`, for a giant too (its own heat flux radiated with the
    starlight is inside the stored G, GEN.1030, so the pair reproduces the temperature, GEN.412); the air solve stops at `surface_air.MAX_SURFACE_TEMPERATURE_K`."""
    return min(climate_draws.surface_temperature_k(planet.bond_albedo, planet.greenhouse_factor, flux), surface_air.MAX_SURFACE_TEMPERATURE_K)


@pytest.mark.parametrize("zone", ["h", "e", "c"])
def test_every_body_stores_the_draws_its_temperature_is_made_of(host_star, zone):
    for seed in range(5):
        with draw.bound(seed):
            planet = Planet(SystemConfig(), host_star, host_star.habitable_zone,
                            plausibility.distance_for_zone(host_star, zone), zone_override=zone, moon_count=0)
        flux = _flux(planet, host_star)
        assert planet.albedo_basis == climate_draws.BASIS_GENERATED
        assert planet.bond_albedo is not None and 0.0 <= planet.bond_albedo < 1.0
        assert planet.greenhouse_factor is not None and planet.greenhouse_factor >= 0.0
        assert planet.surface_temperature == pytest.approx(_temperature_of_draws(planet, flux), rel=TOLERANCE)
        if planet.atmosphere == "None":
            assert planet.greenhouse_factor == 0.0
        if planet.body_type == "t":
            low, high = climate_draws.albedo_range(planet.body_type, planet.water_mass_fraction, flux, getattr(planet.star, "temperature", None))
            assert low <= planet.bond_albedo <= high


def test_a_giant_stores_the_pair_that_gives_its_temperature_with_its_own_heat_in_the_greenhouse(host_star):
    """GEN.1030: a giant radiates its own heat as well, which the air's factor leaves out of T_e; the stored G is against the starlight alone."""
    for seed in range(4):
        with draw.bound(seed):
            planet = Planet(SystemConfig(), host_star, host_star.habitable_zone, plausibility.distance_for_zone(host_star, "c"),
                            mass=300.0 * constants.EARTH_MASS_TO_KG, zone_override="c", moon_count=0)
        assert planet.body_type == "g"
        low, high = climate_draws.GIANT_BOND_ALBEDO_RANGE
        assert low <= planet.bond_albedo <= high
        assert planet.surface_temperature == pytest.approx(
            climate_draws.surface_temperature_k(planet.bond_albedo, planet.greenhouse_factor, _flux(planet, host_star)), rel=TOLERANCE)


def _old_convention_row(albedo, factor, flux, gravity_g, molar_kg_mol=0.029):
    """The columns of a body saved before the albedo and greenhouse columns existed: its scale height keeps the no-air temperature
    (`H = R T_e / (M g)`, the barometric law the generator then used), its temperature is `T_e (1 + G)^(1/4)`."""
    equilibrium = climate_draws.equilibrium_temperature_k(albedo, flux)
    gravity = gravity_g * constants.EARTH_GRAVITY
    height_m = constants.R * equilibrium / (molar_kg_mol * gravity)
    pressure = 101325.0
    return {
        "atmosphere": "a mix of nitrogen and oxygen", "surface_temperature_k": climate_draws.surface_temperature_k(albedo, factor, flux),
        "scale_height_km": height_m / 1000.0, "atm_molar_density": molar_kg_mol, "gravity_g": gravity_g,
        "atm_density": 2.0 * pressure / (gravity * height_m * gravity_g), "atmospheric_pressure_pa": pressure,
    }


@pytest.mark.parametrize("albedo, factor, gravity_g, flux", [
    (0.12, 0.0, 1.0, 1361.0), (0.30, 0.8, 1.0, 1361.0), (0.35, 40.0, 0.9, 1361.0), (0.05, 100.0, 1.3, 2400.0),
    (0.60, 0.05, 0.4, 600.0), (0.25, 3.0, 2.0, 100.0),
])
def test_the_inverse_returns_the_drawn_values_of_a_body_saved_before_the_columns(albedo, factor, gravity_g, flux):
    row = _old_convention_row(albedo, factor, flux, gravity_g)
    recovered_albedo, recovered_factor, reason = climate_draws.from_stored(row, flux)
    assert reason is None
    assert recovered_albedo == pytest.approx(albedo, rel=1e-9)
    assert recovered_factor == pytest.approx(factor, rel=1e-9, abs=1e-9)


def test_the_inverse_says_when_it_has_no_solution():
    base = {"surface_temperature_k": 280.0, "atmosphere": "N2", "scale_height_km": 8.0, "atm_molar_density": 0.029,
            "gravity_g": 1.0, "atm_density": 1.2, "atmospheric_pressure_pa": 101325.0}
    flux = 1361.0
    assert climate_draws.from_stored({**base, "surface_temperature_k": 2.725}, flux)[2] == "the surface sits on the background floor"
    assert "no starlight" in climate_draws.from_stored(base, 0.0)[2]
    assert "missing" in climate_draws.from_stored({**base, "scale_height_km": None}, flux)[2]
    assert "barometric law" in climate_draws.from_stored({**base, "atmospheric_pressure_pa": 1e9}, flux)[2]
    # a scale height that makes the no-air temperature hotter than the sun allows: the albedo would be negative
    assert "albedo" in climate_draws.from_stored({**base, "scale_height_km": 40.0}, flux)[2]


def test_the_backfill_draw_keeps_the_stored_temperature_and_is_seeded():
    row = {"body_type": "t", "water_mass_fraction": 1e-6, "surface_temperature_k": 288.0}
    flux = 1361.0
    first = climate_draws.backfill_draw(row, flux, "planets:7")
    assert first == climate_draws.backfill_draw(row, flux, "planets:7")
    assert first != climate_draws.backfill_draw(row, flux, "planets:8")
    albedo, factor = first
    low, high = climate_draws.albedo_range("t", 1e-6, flux)
    assert low <= albedo <= high
    assert factor >= 0.0
    assert climate_draws.surface_temperature_k(albedo, factor, flux) == pytest.approx(288.0, rel=1e-9)


def test_a_backfill_draw_never_gives_a_negative_greenhouse():
    # A cold body far from its star: the range's low albedo would need G < 0, so the albedo is lifted.
    row = {"body_type": "t", "water_mass_fraction": 1e-6, "surface_temperature_k": 40.0}
    flux = 0.5
    for seed in range(20):
        albedo, factor = climate_draws.backfill_draw(row, flux, f"moons:{seed}")
        assert factor >= -1e-9 and math.isfinite(albedo)
        assert climate_draws.surface_temperature_k(albedo, max(factor, 0.0), flux) == pytest.approx(40.0, rel=1e-6)


def test_the_albedo_range_follows_the_water_and_the_cold():
    sun = 1361.0
    assert climate_draws.albedo_range("g", 0.3, sun) == climate_draws.GIANT_BOND_ALBEDO_RANGE == (0.29, 0.503)
    assert climate_draws.albedo_range("t", 1e-6, sun) == climate_draws.DEFAULT_ALBEDO_RANGE
    assert climate_draws.albedo_range("t", 1e-3, 0.0) == climate_draws.DEFAULT_ALBEDO_RANGE
    assert climate_draws.albedo_range("t", 1e-3, sun) == tuning.ALBEDO_WET_RANGE
    assert climate_draws.albedo_range("t", 1e-3, 0.3 * 1361.0) == tuning.ALBEDO_FROZEN_RANGE


def test_a_giant_backfill_takes_the_open_bracket():
    assert climate_draws.backfill_albedo_range({"body_type": "g"}, 50.0) == climate_draws.GIANT_ALBEDO_RANGE
