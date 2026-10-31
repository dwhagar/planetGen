# tests/test_rogue_habitability.py

"""
A rogue planet's habitability score (`planetgen.physics.rogue_habitability`, schema v153): worked out from its stored surface
conditions with nothing drawn, so a saved row, a rebuilt object and a migrated old row all give the same numbers.
"""

import types

import pytest

from planetgen import tuning
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.physics import atmosphere, constants, habitability_world, rogue_habitability, rogue_surface
from planetgen.util import draw

EARTH = constants.EARTH_MASS_TO_KG


def _planet(**fields):
    base = dict(planet_type="t", planet_class="M", mass_kg=EARTH, radius_km=6371.0, age_gy=4.5,
                internal_heat_flux_w_m2=0.087, surface_regime="ice-shell-ocean", surface_temperature_k=35.0,
                surface_pressure_pa=0.3, ice_shell_thickness_km=13.0, ocean_depth_km=120.0, has_liquid_water=True,
                hp_ice_km=None)
    base.update(fields)
    return types.SimpleNamespace(**base)


def test_an_ocean_under_ice_scores_on_its_own_heat_and_the_ground_dose():
    planet = _planet()
    rogue_habitability.generate(planet)
    given = rogue_habitability.values(planet)
    assert given["hydrosphere"] == "ice-covered ocean" and given["ocean_class"] == "neutral"
    assert given["energy_flux_w_m2"] == pytest.approx(habitability_world.CHEMICAL_SHARE * 0.087)
    assert given["surface_dose_msv_yr"] > given["dose_ground_msv_yr"] > 0
    assert planet.phi4 is not None and 0.0 <= planet.phi4 <= 1.0
    assert planet.tier_temperature == 3 and planet.tier_pressure == 3  # 35 K and 0.3 Pa are Red
    assert planet.eq_body >= 3 and planet.eq_rad_suit > 0  # a hard suit and shielding at best
    assert not hasattr(planet, "band1_mean_k")  # no star, so no latitude bands
    assert (planet.life_microbial, planet.life_flora, planet.life_fauna, planet.life_sentient) == (0, 0, 0, 0)
    assert set(rogue_habitability.ROGUE_FIELDS) <= set(vars(planet))


def test_the_water_state_follows_the_regime_and_the_ice():
    sealed = _planet(hp_ice_km=30.0)
    assert rogue_habitability.values(sealed)["ocean_class"] == "ice-sealed"
    hycean = _planet(surface_regime="hycean", surface_temperature_k=300.0, surface_pressure_pa=5e6, ice_shell_thickness_km=None)
    given = rogue_habitability.values(hycean)
    assert given["hydrosphere"] == "hycean" and atmosphere.partial_kpa(given, "h2") == pytest.approx(5000.0)
    ice = _planet(has_liquid_water=False, ocean_depth_km=None)
    assert rogue_habitability.values(ice)["hydrosphere"] == "ice"
    dry = _planet(has_liquid_water=False, ocean_depth_km=None, ice_shell_thickness_km=None, surface_regime="bare-rock",
                  surface_pressure_pa=0.0)
    given = rogue_habitability.values(dry)
    assert given["hydrosphere"] == "dry" and given["water_mass_fraction"] is None
    frozen_air = _planet(has_liquid_water=False, ocean_depth_km=None, ice_shell_thickness_km=None,
                         surface_regime="frozen-atmosphere", surface_pressure_pa=2.0)
    assert atmosphere.partial_kpa(rogue_habitability.values(frozen_air), "n2") == pytest.approx(0.002)


def test_the_water_fraction_is_worked_back_from_the_column():
    planet = _planet()
    state, fraction = rogue_habitability._hydrosphere(planet)
    from planetgen.physics import hydrosphere
    column = 120.0 + 13.0 * hydrosphere.ICE_DENSITY_KG_M3 / hydrosphere.LIQUID_WATER_DENSITY_KG_M3
    assert state == "ice-covered ocean"
    assert hydrosphere.water_layer_depth_km(planet.mass_kg, planet.radius_km, fraction) == pytest.approx(column)


def test_a_gas_giant_or_brown_dwarf_scores_nothing():
    giant = _planet(planet_type="g", planet_class=None, surface_temperature_k=120.0)
    rogue_habitability.generate(giant)
    assert all(getattr(giant, name) is None for name in rogue_habitability.ROGUE_FIELDS)


def test_the_score_draws_nothing_and_repeats():
    draw.set_run_seed(5)
    first = _planet()
    rogue_habitability.generate(first)
    again = _planet()
    rogue_habitability.generate(again)
    assert {n: getattr(first, n) for n in rogue_habitability.ROGUE_FIELDS} == {n: getattr(again, n) for n in rogue_habitability.ROGUE_FIELDS}


def test_a_hycean_rogue_scores_above_a_bare_ice_rogue():
    hycean = _planet(surface_regime="hycean", surface_temperature_k=300.0, surface_pressure_pa=5e6, ice_shell_thickness_km=None)
    ice = _planet(has_liquid_water=False, ocean_depth_km=None, surface_temperature_k=35.0, surface_regime="bare-rock",
                  surface_pressure_pa=0.0)
    rogue_habitability.generate(hycean)
    rogue_habitability.generate(ice)
    assert hycean.phi4_water > ice.phi4_water
    assert hycean.phi4 >= ice.phi4


def test_every_generated_rogue_planet_is_scored_unless_it_is_a_giant():
    seen = set()
    for mass_bin in tuning.ROGUE_PLANET_MASS_BIN_CHOICES:
        for seed in range(6):
            with draw.bound(seed):
                planet = RoguePlanet(SystemConfig(), mass_bin=mass_bin)
            fields = [getattr(planet, name) for name in rogue_habitability.ROGUE_FIELDS]
            if planet.planet_type == "g":
                assert all(value is None for value in fields)
            else:
                optional = ("shelter_column_g_cm2",)  # NULL when no shelter applies
                assert all(getattr(planet, name) is not None for name in rogue_habitability.ROGUE_FIELDS if name not in optional)
                assert planet.energy_flux_w_m2 == pytest.approx(habitability_world.CHEMICAL_SHARE * planet.internal_heat_flux_w_m2)
                seen.add(planet.surface_regime)
    assert seen  # at least one terrestrial regime came up


def test_a_rebuilt_rogue_planet_keeps_its_score():
    with draw.bound(3):
        planet = RoguePlanet(SystemConfig(), mass_bin="terrestrial")
    data = planet.to_dict()
    assert data["phi4"] == planet.phi4
    rebuilt = RoguePlanet.from_dict(data, SystemConfig())
    assert all(getattr(rebuilt, name) == getattr(planet, name) for name in rogue_habitability.ROGUE_FIELDS)
