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
from planetgen.physics import climate_draws, constants
from planetgen.util import draw

TOLERANCE = 1e-9
"""float: The relative error the inverse may have (it is exact up to floating-point rounding)."""


@pytest.fixture(scope="module")
def host_star():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    return Star(cfg)


def _flux(planet, star):
    distance_km = float(planet.distance) * constants.AU_TO_KM
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


@pytest.mark.parametrize("planet_class", sorted(tuning.PLANET_CLASSES))
def test_every_class_stores_the_draws_its_temperature_is_made_of(host_star, planet_class):
    data = tuning.PLANET_CLASSES[planet_class]
    zone = next((z for z in ("e", "h", "c") if data.get(z)), None)
    if zone is None:
        pytest.skip("a class no zone offers")
    for seed in range(5):
        with draw.bound(seed):
            planet = Planet(SystemConfig(), host_star, host_star.habitable_zone,
                            plausibility.distance_for_zone(host_star, zone),
                            planet_class=planet_class, zone_override=zone, moon_count=0)
        flux = _flux(planet, host_star)
        assert planet.albedo_basis == climate_draws.BASIS_GENERATED
        assert planet.bond_albedo is not None and 0.0 <= planet.bond_albedo < 1.0
        assert planet.greenhouse_factor is not None and planet.greenhouse_factor >= 0.0
        assert planet.surface_temperature == pytest.approx(
            climate_draws.surface_temperature_k(planet.bond_albedo, planet.greenhouse_factor, flux), rel=TOLERANCE)
        if planet.atmosphere == "None":
            assert planet.greenhouse_factor == 0.0
        low, high = climate_draws.class_albedo_range(planet.planet_class, planet.body_type)
        if planet_class in tuning.PLANET_CLASSES and "albedo_range" in data:
            assert low <= planet.bond_albedo <= high


@pytest.mark.parametrize("planet_class", ["F", "G", "N", "D", "P", "K", "S"])
def test_the_inverse_returns_the_drawn_values(host_star, planet_class):
    data = tuning.PLANET_CLASSES[planet_class]
    zone = next((z for z in ("e", "h", "c") if data.get(z)), None)
    recovered = 0
    for seed in range(12):
        with draw.bound(seed):
            planet = Planet(SystemConfig(), host_star, host_star.habitable_zone,
                            plausibility.distance_for_zone(host_star, zone),
                            planet_class=planet_class, zone_override=zone, moon_count=0)
        flux = _flux(planet, host_star)
        albedo, factor, reason = climate_draws.from_stored(_stored_row(planet), flux)
        if albedo is None:
            # The solution is refused only where the temperature sits on the floor or the row is not a plain one.
            assert reason
            continue
        recovered += 1
        assert albedo == pytest.approx(planet.bond_albedo, rel=1e-6, abs=1e-9)
        assert factor == pytest.approx(planet.greenhouse_factor, rel=1e-6, abs=1e-9)
    assert recovered, f"class {planet_class}: no body of 12 could be inverted"


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


def test_the_class_draw_keeps_the_stored_temperature_and_is_seeded():
    row = {"planet_class": "G", "body_type": "t", "surface_temperature_k": 288.0}
    flux = 1361.0
    first = climate_draws.class_draw(row, flux, "planets:7")
    assert first == climate_draws.class_draw(row, flux, "planets:7")
    assert first != climate_draws.class_draw(row, flux, "planets:8")
    albedo, factor = first
    low, high = climate_draws.class_albedo_range("G", "t")
    assert low <= albedo <= high
    assert factor >= 0.0
    assert climate_draws.surface_temperature_k(albedo, factor, flux) == pytest.approx(288.0, rel=1e-9)


def test_a_class_draw_never_gives_a_negative_greenhouse():
    # A cold body far from its star: the class range's low albedo would need G < 0, so the albedo is lifted.
    row = {"planet_class": "G", "body_type": "t", "surface_temperature_k": 40.0}
    flux = 0.5
    for seed in range(20):
        albedo, factor = climate_draws.class_draw(row, flux, f"moons:{seed}")
        assert factor >= -1e-9 and math.isfinite(albedo)
        assert climate_draws.surface_temperature_k(albedo, max(factor, 0.0), flux) == pytest.approx(40.0, rel=1e-6)


def test_a_giant_with_no_range_takes_the_open_bracket():
    assert climate_draws.class_albedo_range("no such class", "g") == climate_draws.GIANT_ALBEDO_RANGE
    assert climate_draws.class_albedo_range("no such class", "t") == climate_draws.DEFAULT_ALBEDO_RANGE
