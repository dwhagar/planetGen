"""
Hydrosphere and ocean chemistry (GEN.88, `planetgen.physics.hydrosphere`):
the design doc's section 5 numbers, the ocean-class rules and what every
generated body stores.
"""

import types

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants, hydrosphere as hy, water_layer as wl
from planetgen.util import draw
from planetgen.util.random import log_uniform

from tests.air_support import air, fold_air
from tests.fuzz_support import deterministic_entropy

EARTH = constants.EARTH_MASS_TO_KG


@pytest.mark.parametrize("bottom_k, pressure_mpa", [
    (251.2, 210.0), (256.2, 350.0), (272.1, 620.0), (286.8, 800.0), (300.4, 1000.0), (327.0, 1500.0), (353.8, 2200.0),
])
def test_the_liquidus_follows_the_recomputed_melting_curves(bottom_k, pressure_mpa):
    assert wl.liquidus_pressure_pa(bottom_k) == pytest.approx(pressure_mpa * 1e6, rel=1e-9)


def test_the_deepest_liquid_at_earth_gravity_is_about_seventy_kilometres():
    """The design doc 5.2 had 73 km at a fixed 1050 kg/m3; the water is lighter than that over the first 700 MPa."""
    deepest = wl.column(EARTH, 6371.0, 0.02 * EARTH, bottom_temperature_k=280.0)
    assert deepest.liquid_km == pytest.approx(65.5, rel=0.03)
    assert deepest.hp_ice_km > 0 and deepest.base_pressure_pa > wl.liquidus_pressure_pa(280.0)


def test_earths_water_covers_seven_tenths():
    layer = hy.water_layer_depth_km(EARTH, 6371.0, hy.EARTH_OCEAN_FRACTION)
    assert layer == pytest.approx(2.74, rel=0.03)
    assert hy.ocean_fraction(layer, constants.EARTH_GRAVITY) == pytest.approx(0.72, abs=0.02)
    assert hy.ocean_fraction(hy.water_layer_depth_km(EARTH, 6371.0, 2e-3), constants.EARTH_GRAVITY) == 1.0


def test_europas_lid_and_tides():
    assert hy.ice_shell_thickness_km(0.05, 100.0) == pytest.approx(11.4, abs=0.2)
    assert hy.pressure_melted_shell_km(0.05, 100.0, 1.3) < hy.ice_shell_thickness_km(0.05, 100.0)
    europa = hy.tidal_heating_w(3e-4, constants.JUPITER_MASS_TO_KG, 1560.8, 6.709e8, 0.009)
    assert europa == pytest.approx(1.4e11, rel=0.15)


def test_the_lid_counts_by_mass_and_high_pressure_ice_seals_the_ocean():
    # A lid of 10 km of ice is 9.17e6 kg/m2: a column lighter than that freezes solid.
    area = 4 * 3.141592653589793 * (6371.0e3) ** 2
    light = wl.column(EARTH, 6371.0, 9.0e6 * area, shell_km=10.0, bottom_temperature_k=273.15)
    assert light.liquid_km == 0.0 and light.hp_ice_km == 0.0 and light.ice_shell_km == pytest.approx(9.0e6 / 917.0 / 1000.0)
    open_ocean = wl.column(EARTH, 6371.0, 2.0e7 * area, shell_km=10.0, bottom_temperature_k=273.15)
    assert open_ocean.ice_shell_km == 10.0 and open_ocean.liquid_km == pytest.approx(10.34, rel=0.02) and open_ocean.hp_ice_km == 0.0
    deep = wl.column(EARTH, 6371.0, 0.05 * EARTH, bottom_temperature_k=300.0)
    assert deep.liquid_km == pytest.approx(88.1, rel=0.03) and deep.hp_ice_km > 0


def _ocean(**gases):
    return types.SimpleNamespace(mantle_redox=gases.pop("redox", "intermediate"), **air(**gases))


def test_ocean_classes_follow_the_rules():
    assert hy.ocean_class(_ocean(), 5.0, 1e-3, 1.0) == "ice-sealed"
    assert hy.ocean_class(_ocean(), 0.0, 1e-5, 0.1) == "chloride brine"
    acid = _ocean(redox="oxidized", p_so2_kpa=1.0, p_co2_kpa=10.0)
    assert hy.ocean_class(acid, 0.0, 1e-3, 0.5) == "acid sulfate"
    soda = _ocean(redox="reduced", p_so2_kpa=0.0, p_co2_kpa=50.0)
    assert hy.ocean_class(soda, 0.0, 1e-3, 0.95) == "soda"
    assert hy.ocean_class(soda, 0.0, 1e-3, 0.5) == "neutral"


def _body(planet_class, temperature, pressure_pa, mass_earth=1.0, radius_km=6371.0, **extra):
    fields = dict(body_type="t", planet_class=planet_class, mass=mass_earth * EARTH, radius=radius_km,
                  water_mass_fraction=log_uniform(*hy.WATER_MASS_FRACTION[planet_class]),
                  gravity=mass_earth / (radius_km / 6371.0) ** 2, surface_temperature=temperature,
                  atmospheric_pressure=pressure_pa, is_moon=False, mantle_redox="intermediate",
                  p_co2_kpa=0.04, p_so2_kpa=0.0, p_h2_kpa=0.0)
    fields.update(extra)
    return types.SimpleNamespace(**fold_air(fields))


def test_states_by_temperature_and_pressure():
    draw.set_run_seed(88)
    earth = _body("M", 288.0, 101325.0)
    hy.generate_hydrosphere(earth, 0.087)
    assert earth.hydrosphere == "surface ocean" and earth.ocean_class in hy.OCEAN_CLASSES
    assert earth.ocean_fraction + earth.land_fraction == pytest.approx(1.0)
    assert 6.5 <= earth.ocean_ph <= 8.5 or earth.ocean_class != "neutral"
    venus = _body("M", 740.0, 9.2e6)
    hy.generate_hydrosphere(venus, 0.087)
    assert venus.hydrosphere == "vapour" and venus.ocean_depth_km is None
    moon = _body("C", 250.0, 0.0)
    hy.generate_hydrosphere(moon, 0.02)
    assert moon.hydrosphere == "dry" and moon.land_fraction == 1.0
    snowball = _body("P", 200.0, 50000.0)
    hy.generate_hydrosphere(snowball, 0.087)
    assert snowball.hydrosphere in ("ice", "ice-covered ocean") and snowball.ice_shell_km > 0
    hycean = _body("O", 300.0, 1e6, p_h2_kpa=900.0)
    hy.generate_hydrosphere(hycean, 0.087)
    assert hycean.hydrosphere == "hycean"
    giant = types.SimpleNamespace(body_type="g", water_mass_fraction=0.3)
    hy.generate_hydrosphere(giant, 1.0)
    assert giant.water_mass_fraction == 0.3
    assert all(getattr(giant, field) is None for field in hy.HYDROSPHERE_FIELDS if field != "water_mass_fraction")


def test_ocean_worlds_are_mostly_ice_sealed():
    draw.set_run_seed(9)
    worlds = [_body("O", 290.0, 101325.0) for _ in range(100)]
    for world in worlds:
        hy.generate_hydrosphere(world, 0.087)
    assert all(world.ocean_fraction == 1.0 for world in worlds)
    sealed = [world for world in worlds if world.ocean_class == "ice-sealed"]
    assert len(sealed) > 50
    assert all(world.hp_ice_km > 0 and world.phosphorus == "starved" for world in sealed)


def test_generated_bodies_store_a_consistent_hydrosphere():
    seen = set()
    for seed in range(15):
        with deterministic_entropy(seed):
            system = StarSystem(system_config=SystemConfig())
        for planet in system.planets:
            if not hasattr(planet, "moons"):
                continue
            for body in [planet, *planet.moons]:
                if body.body_type == "g":
                    assert all(getattr(body, field) is None for field in hy.HYDROSPHERE_FIELDS
                               if field != "water_mass_fraction")
                    continue
                seen.add(body.hydrosphere)
                assert body.hydrosphere in hy.HYDROSPHERE_STATES
                assert body.water_mass_fraction >= 0  # 0: its early runaway greenhouse took it all (GEN.639)
                if body.water_mass_fraction == 0:
                    assert body.hydrosphere == "dry"
                assert body.ocean_fraction + body.land_fraction == pytest.approx(1.0)
                liquid = body.hydrosphere in ("ice-covered ocean", "surface ocean", "hycean")
                assert (body.ocean_depth_km is not None) == liquid == (body.ocean_class is not None)
                if liquid:
                    (ph_lo, ph_hi), (aw_lo, aw_hi) = hy.OCEAN_CHEMISTRY[body.ocean_class]
                    assert ph_lo <= body.ocean_ph <= ph_hi and aw_lo <= body.water_activity <= aw_hi
                    assert body.phosphorus in ("high", "limited", "starved")
    assert {"ice", "dry", "vapour"} <= seen


def _moon(planet_class, mass_earth, radius_km, distance_au, **extra):
    return _body(planet_class, 100.0, 0.0, mass_earth=mass_earth, radius_km=radius_km, is_moon=True,
                 primary_mass_kg=constants.JUPITER_MASS_TO_KG, distance=distance_au, **extra)


def test_an_icy_moon_under_ice_keeps_its_tide_and_a_planet_or_a_warm_moon_keeps_none():
    for seed in range(8):
        draw.set_run_seed(seed)
        europa = _moon("D", 8.0e-3, 1560.8, 6.709e8 / constants.AU_TO_M)
        hy.generate_hydrosphere(europa, 0.002)
        # Europa's own 1.4e11 W over its surface is 4.6e-3 W/m2; the forced eccentricity is 0.001 to 0.01, so a tenth to 2 times
        assert 5e-5 < europa.tidal_flux_w_m2 < 2e-2
        draw.set_run_seed(seed)
        enceladus = _moon("D", 1.8e-5, 252.1, 2.38e8 / constants.AU_TO_M)
        hy.generate_hydrosphere(enceladus, 0.002)
        assert enceladus.tidal_flux_w_m2 > 0.0  # a small moon is boosted by its resonance
    planet = _body("D", 100.0, 0.0)
    hy.generate_hydrosphere(planet, 0.002)
    assert planet.tidal_flux_w_m2 is None
    warm = _moon("D", 8.0e-3, 1560.8, 6.709e8 / constants.AU_TO_M)
    warm.surface_temperature = 290.0
    warm.atmospheric_pressure = 101325.0
    hy.generate_hydrosphere(warm, 0.002)
    assert warm.tidal_flux_w_m2 is None  # a liquid surface needs no ice lid, so no tide is drawn


def test_the_tide_is_in_the_chemical_energy_of_an_ocean_under_ice():
    from planetgen.physics import habitability_world as hw
    body = types.SimpleNamespace(body_type="t", surface_temperature=100.0, tidal_flux_w_m2=0.004)
    values = {"hydrosphere": "ice-covered ocean"}
    heat = 0.002 + body.tidal_flux_w_m2
    assert hw.energy_flux_w_m2(values, None, None, heat) == pytest.approx(hw.CHEMICAL_SHARE * 0.006)
    assert hw.energy_flux_w_m2(values, None, None, heat) > 2.5 * hw.energy_flux_w_m2(values, None, None, 0.002)
