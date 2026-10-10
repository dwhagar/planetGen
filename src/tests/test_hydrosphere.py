"""
Hydrosphere and ocean chemistry (GEN.88, `planetgen.physics.hydrosphere`):
the design doc's section 5 numbers, the ocean-class rules and what every
generated body stores.
"""

import types

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants, hydrosphere as hy
from planetgen.util import draw

from tests.fuzz_support import deterministic_entropy

EARTH = constants.EARTH_MASS_TO_KG


@pytest.mark.parametrize("gravity, bottom_k, depth_km", [
    (9.8, 280.0, 73.0), (9.8, 300.0, 111.0), (9.8, 350.0, 204.0), (1.3, 280.0, 553.0), (20.0, 300.0, 54.0),
])
def test_the_deepest_liquid_matches_the_design_table(gravity, bottom_k, depth_km):
    assert hy.max_liquid_depth_km(bottom_k, gravity) == pytest.approx(depth_km, rel=0.03)


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


def test_the_lid_counts_as_water_and_high_pressure_ice_seals_the_ocean():
    # 10 km of ice is 9.17 km of water: a 9 km column freezes solid.
    assert hy.split_column(9.0, 10.0, 273.15, 9.8) == (pytest.approx(9.0 * hy.ICE_TO_WATER), None, 0.0)
    shell, liquid, hp = hy.split_column(20.0, 10.0, 273.15, 9.8)
    assert liquid == pytest.approx(20.0 - 10.0 / hy.ICE_TO_WATER) and hp == 0.0
    shell, liquid, hp = hy.split_column(500.0, 0.0, 300.0, 9.8)
    assert liquid == pytest.approx(111.0, rel=0.03)
    assert hp == pytest.approx((500.0 - liquid) / 1.3, rel=1e-6)


def _ocean(**gases):
    return types.SimpleNamespace(mantle_redox=gases.pop("redox", "intermediate"), **gases)


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
                  gravity=mass_earth / (radius_km / 6371.0) ** 2, surface_temperature=temperature,
                  atmospheric_pressure=pressure_pa, is_moon=False, mantle_redox="intermediate",
                  p_co2_kpa=0.04, p_so2_kpa=0.0, p_h2_kpa=0.0)
    fields.update(extra)
    return types.SimpleNamespace(**fields)


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
    giant = types.SimpleNamespace(body_type="g")
    hy.generate_hydrosphere(giant, 1.0)
    assert all(getattr(giant, field) is None for field in hy.HYDROSPHERE_FIELDS)


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
                    assert all(getattr(body, field) is None for field in hy.HYDROSPHERE_FIELDS)
                    continue
                seen.add(body.hydrosphere)
                assert body.hydrosphere in hy.HYDROSPHERE_STATES
                assert body.water_mass_fraction > 0
                assert body.ocean_fraction + body.land_fraction == pytest.approx(1.0)
                liquid = body.hydrosphere in ("ice-covered ocean", "surface ocean", "hycean")
                assert (body.ocean_depth_km is not None) == liquid == (body.ocean_class is not None)
                if liquid:
                    (ph_lo, ph_hi), (aw_lo, aw_hi) = hy.OCEAN_CHEMISTRY[body.ocean_class]
                    assert ph_lo <= body.ocean_ph <= ph_hi and aw_lo <= body.water_activity <= aw_hi
                    assert body.phosphorus in ("high", "limited", "starved")
    assert {"ice", "dry", "vapour"} <= seen
