"""Tests for physics/surface_composition.py (GEN.441): shares sum to one, agree with the bulk fractions, the temperature and the air."""

import math

import pytest

from planetgen.generation import possibility_table as pt
from planetgen.physics import surface_composition as sc
from planetgen.util import draw


@pytest.fixture(autouse=True)
def seeded():
    draw.set_run_seed(441)


def earth(**overrides):
    values = dict(planet_class="M", surface_temperature_k=288.0, pressure_pa=101325.0, o2_partial_pa=21000.0,
                  water_mass_fraction=2.3e-4, rock_mass_fraction=0.67, metal_mass_fraction=0.33, rock_feo_fraction=0.08,
                  hydrosphere="surface ocean", ocean_fraction=0.71, internal_heat_flux_w_m2=0.09, age_gy=4.5,
                  mobile_lid_chance=0.5)
    values.update(overrides)
    return sc.SurfaceInputs(**values)


def check_sums(result):
    assert math.fsum(result.mix.values()) == pytest.approx(1.0, abs=1e-9)
    assert math.fsum(result.components.values()) == pytest.approx(1.0, abs=1e-9)
    assert all(share >= 0.0 for share in result.mix.values())
    assert all(share >= 0.0 for share in result.components.values())


def test_tables_validate():
    pt.REGISTRY.validate()
    for row in sc.GROUND_TYPES.rows:
        assert math.fsum(row.normalised_fractions().values()) == pytest.approx(1.0)


def test_giants_have_no_ground():
    for letter in sc.GIANT_CLASSES:
        assert sc.draw_surface(earth(planet_class=letter)) is None


@pytest.mark.parametrize("planet_class", sorted(set("ABCDEFGHKLMNOPQSV")))
def test_every_class_sums_to_one(planet_class):
    for seed in range(20):
        draw.set_run_seed(seed)
        result = sc.draw_surface(earth(planet_class=planet_class))
        check_sums(result)
        assert result.grain in sc.GRAINS


def test_ocean_world_has_ocean_and_salts():
    result = sc.draw_surface(earth())
    assert result.mix[sc.OCEAN] == pytest.approx(0.71)
    assert result.components["water_liquid"] == pytest.approx(0.71 * (1 - sc.RANGES["ocean_salt"]), rel=1e-6)
    assert result.components["halite"] >= 0.71 * sc.RANGES["ocean_salt"]


def test_ice_needs_cold_enough_surface():
    warm = sc.draw_surface(earth(surface_temperature_k=300.0, hydrosphere="ice", ocean_fraction=0.9))
    assert sc.ICE not in warm.mix and "water_ice" not in warm.components
    cold = sc.draw_surface(earth(surface_temperature_k=250.0, hydrosphere="ice", ocean_fraction=0.9))
    assert cold.mix[sc.ICE] == pytest.approx(0.9)
    assert "water_ice" in cold.components


def test_ice_rich_frozen_body_is_all_ice():
    result = sc.draw_surface(earth(planet_class="D", surface_temperature_k=100.0, pressure_pa=0.0, water_mass_fraction=0.4,
                                   rock_mass_fraction=0.6, hydrosphere="ice", ocean_fraction=0.2))
    assert result.mix == {sc.ICE: 1.0}
    assert set(result.components) == {"water_ice"}


def test_other_ices_stand_only_where_cold_enough():
    shares = {"n2_ice": 0.6, "ch4_ice": 0.4}
    cold = sc.draw_surface(earth(planet_class="D", surface_temperature_k=30.0, pressure_pa=1.0, water_mass_fraction=0.3,
                                 ice_mass_fraction=0.2, rock_mass_fraction=0.5, volatile_shares=shares, hydrosphere="ice",
                                 ocean_fraction=0.1))
    assert {"water_ice", "n2_ice", "ch4_ice"} <= set(cold.components)
    warm = sc.draw_surface(earth(planet_class="D", surface_temperature_k=120.0, pressure_pa=1.0, water_mass_fraction=0.3,
                                 ice_mass_fraction=0.2, rock_mass_fraction=0.5, volatile_shares=shares, hydrosphere="ice",
                                 ocean_fraction=0.1))
    assert "n2_ice" not in warm.components and "ch4_ice" not in warm.components
    assert "water_ice" in warm.components


def test_water_ice_vacuum_and_air_limits():
    assert sc.ice_stable("water_ice", 200.0, 0.0)
    assert not sc.ice_stable("water_ice", 230.0, 0.0)
    assert sc.ice_stable("water_ice", 260.0, 101325.0)
    assert not sc.ice_stable("water_ice", 280.0, 101325.0)
    assert sc.ice_stable("co2_ice", 140.0, 600.0, partial_pa=570.0)
    assert not sc.ice_stable("co2_ice", 170.0, 600.0, partial_pa=570.0)


def test_modifiers_follow_the_air_and_the_water():
    dry = sc.draw_surface(earth(hydrosphere="dry", ocean_fraction=0.0, o2_partial_pa=0.0, mobile_lid_chance=0.0))
    assert not {"felsic", "clay", "oxidised", "evaporite", "sulfur"} & set(dry.mix)
    for seed in range(30):
        draw.set_run_seed(seed)
        wet = sc.draw_surface(earth(mobile_lid_chance=1.0))
        assert {"felsic", "clay", "oxidised", "evaporite"} <= set(wet.mix)


def test_tidal_sulfur():
    result = sc.draw_surface(earth(planet_class="D", hydrosphere="dry", ocean_fraction=0.0, tidal_heat_flux_w_m2=2.0,
                                   internal_heat_flux_w_m2=0.01, o2_partial_pa=0.0, mobile_lid_chance=0.0,
                                   surface_temperature_k=110.0, pressure_pa=0.0))
    assert "sulfur" in result.mix


def test_dominant_ground_rules():
    base = dict(hydrosphere="dry", ocean_fraction=0.0, o2_partial_pa=0.0, mobile_lid_chance=0.0)
    seen = {}
    for seed in range(60):
        draw.set_run_seed(seed)
        seen.setdefault("basaltic", sc.draw_surface(earth(internal_heat_flux_w_m2=1.0, **base)).ground_type)
        old = sc.draw_surface(earth(internal_heat_flux_w_m2=0.0, rock_feo_fraction=0.18, **base)).ground_type
        seen.setdefault("old", set()).add(old)
        assert sc.draw_surface(earth(metal_mass_fraction=0.8, rock_mass_fraction=0.2, **base)).ground_type == "metal"
        assert sc.draw_surface(earth(carbon_rock=True, **base)).ground_type == "carbon"
    assert seen["basaltic"] == "basaltic"
    assert seen["old"] == {"ultramafic", "anorthositic"}


def test_overturn_follows_feo():
    base = dict(hydrosphere="dry", ocean_fraction=0.0, o2_partial_pa=0.0, mobile_lid_chance=0.0, internal_heat_flux_w_m2=0.0)
    ultramafic = {}
    for feo in (0.0, 0.2):
        count = 0
        for seed in range(80):
            draw.set_run_seed(seed)
            count += sc.draw_surface(earth(rock_feo_fraction=feo, **base)).ground_type == "ultramafic"
        ultramafic[feo] = count
    assert ultramafic[0.0] == 0 and ultramafic[0.2] == 80


def test_grain_state():
    base = dict(hydrosphere="dry", ocean_fraction=0.0)
    assert sc.draw_surface(earth(surface_temperature_k=1800.0, **base)).grain == sc.SLAB
    assert sc.draw_surface(earth(age_gy=0.05, **base)).grain in (sc.SLAB, sc.CRUSHED)
    assert sc.draw_surface(earth(**base)).grain == sc.POWDER


def test_determinism():
    draw.set_run_seed(7)
    first = sc.draw_surface(earth())
    draw.set_run_seed(7)
    assert sc.draw_surface(earth()) == first


def test_albedo_is_between_the_ground_types():
    base = dict(hydrosphere="dry", ocean_fraction=0.0, o2_partial_pa=0.0, mobile_lid_chance=0.0, internal_heat_flux_w_m2=1.0)
    result = sc.draw_surface(earth(**base))
    bands, covered = sc.band_albedo(result)
    assert 0.0 < covered <= 1.0
    for value in bands.values():
        assert 0.0 <= value <= 1.0
    snow = sc.draw_surface(earth(surface_temperature_k=200.0, hydrosphere="ice", ocean_fraction=1.0))
    snow_bands, _ = sc.band_albedo(snow)
    assert snow_bands["alb_050_060"] > bands["alb_050_060"]


def test_reservoir_orders_heat_producers():
    def reservoir(name):
        mix = {name: 1.0}
        return sc.crust_reservoir(sc.SurfaceComposition(name, mix, {"olivine": 1.0}, sc.POWDER))[0]

    ultramafic, basaltic, felsic = (reservoir(n) for n in ("ultramafic", "basaltic", "felsic"))
    for key in ("k_ppm", "th_ppm", "u_ppm"):
        assert ultramafic[key] < basaltic[key] < felsic[key]
    assert sc.crust_reservoir(sc.SurfaceComposition("ice", {sc.ICE: 1.0}, {"water_ice": 1.0}, sc.POWDER)) == (None, 0.0)
