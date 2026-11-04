"""Tests for physics/surface_composition.py (GEN.441): shares sum to one, agree with the bulk fractions, the temperature and the air."""

import math
import types

import pytest

from planetgen.generation import possibility_table as pt
from planetgen.physics import rock_types, surface_composition as sc
from planetgen.util import draw


@pytest.fixture(autouse=True)
def seeded():
    draw.set_run_seed(441)


def earth(**overrides):
    values = dict(surface_temperature_k=288.0, pressure_pa=101325.0, o2_partial_pa=21000.0,
                  water_mass_fraction=2.3e-4, rock_mass_fraction=0.67, metal_mass_fraction=0.33, rock_feo_fraction=0.08,
                  hydrosphere="surface ocean", ocean_fraction=0.71, internal_heat_flux_w_m2=0.09, age_gy=4.5,
                  mobile_lid=True)
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
    assert sc.draw_surface(earth(giant=True)) is None


def test_the_ground_never_reads_the_class_letter():
    """GEN.917, GEN.915: the inputs hold no class letter, and the dominant ground is drawn on conditions alone (every row of its table
    is open to any class)."""
    assert "planet_class" not in sc.SurfaceInputs.__dataclass_fields__
    assert all(row.classes == (pt.ANY_CLASS,) for row in sc.DOMINANT_TABLE.rows)
    for seed in range(20):
        draw.set_run_seed(seed)
        result = sc.draw_surface(earth())
        check_sums(result)
        assert result.grain in sc.GRAINS


@pytest.mark.parametrize("kwargs", [dict(), dict(surface_temperature_k=150.0, pressure_pa=0.0, water_mass_fraction=0.0, hydrosphere="dry",
                                                 ocean_fraction=0.0),
                                    dict(carbon_rock=True, hydrosphere="dry", ocean_fraction=0.0),
                                    dict(metal_mass_fraction=0.9, rock_mass_fraction=0.1, hydrosphere="dry", ocean_fraction=0.0)])
def test_the_dominant_ground_follows_the_conditions(kwargs):
    """A body of these conditions has the ground its conditions say, the same on every seed's stream position."""
    draw.set_run_seed(5)
    first = sc.draw_surface(earth(**kwargs))
    draw.set_run_seed(5)
    assert sc.draw_surface(earth(**kwargs)).ground_type == first.ground_type


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
    result = sc.draw_surface(earth(surface_temperature_k=100.0, pressure_pa=0.0, water_mass_fraction=0.4,
                                   rock_mass_fraction=0.6, hydrosphere="ice", ocean_fraction=0.2))
    assert result.mix == {sc.ICE: 1.0}
    assert set(result.components) == {"water_ice"}


def test_other_ices_stand_only_where_cold_enough():
    shares = {"n2_ice": 0.6, "ch4_ice": 0.4}
    cold = sc.draw_surface(earth(surface_temperature_k=30.0, pressure_pa=1.0, water_mass_fraction=0.3,
                                 ice_mass_fraction=0.2, rock_mass_fraction=0.5, volatile_shares=shares, hydrosphere="ice",
                                 ocean_fraction=0.1))
    assert {"water_ice", "n2_ice", "ch4_ice"} <= set(cold.components)
    warm = sc.draw_surface(earth(surface_temperature_k=120.0, pressure_pa=1.0, water_mass_fraction=0.3,
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
    dry = sc.draw_surface(earth(hydrosphere="dry", ocean_fraction=0.0, o2_partial_pa=0.0, mobile_lid=False))
    assert not {"felsic", "clay", "oxidised", "evaporite", "sulfur"} & set(dry.mix)
    for seed in range(30):
        draw.set_run_seed(seed)
        wet = sc.draw_surface(earth(mobile_lid=True))
        assert {"felsic", "clay", "oxidised", "evaporite"} <= set(wet.mix)


def test_tidal_sulfur():
    result = sc.draw_surface(earth(hydrosphere="dry", ocean_fraction=0.0, tidal_heat_flux_w_m2=2.0,
                                   internal_heat_flux_w_m2=0.01, o2_partial_pa=0.0, mobile_lid=False,
                                   surface_temperature_k=110.0, pressure_pa=0.0))
    assert "sulfur" in result.mix


def test_dominant_ground_rules():
    base = dict(hydrosphere="dry", ocean_fraction=0.0, o2_partial_pa=0.0, mobile_lid=False)
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


def _oxides(mg_si, feo=0.08):
    return rock_types.rock_oxides(feo, 0.045, mg_si, 0.067)


def test_the_melt_share_rises_with_the_heat_flow_and_is_a_mid_ocean_ridge_s_at_earths():
    low, high = sc.MELT_FRACTION_RANGE
    shares = [sc.melt_fraction(flux) for flux in (1e-6, 3e-4, 1e-2, 0.09, 2.0, 1e6)]
    assert shares == sorted(shares) and shares[0] == low and shares[-1] == high
    assert shares[3] == pytest.approx(0.09, abs=0.02), "8 to 10 percent melts make the oceanic crust"


def test_the_volcanic_crust_is_the_melt_of_the_stored_mantle():
    """GEN.474: an Earth-like mantle makes basalt, a silica-richer one andesite or granite, and the types are the rows of `GROUND_TYPES`."""
    base = dict(hydrosphere="dry", ocean_fraction=0.0, o2_partial_pa=0.0, mobile_lid=False, internal_heat_flux_w_m2=0.09)
    seen = {}
    for mg_si in (1.2, 0.7, 0.6):
        seen[mg_si] = sc.draw_surface(earth(rock_oxides=_oxides(mg_si), **base))
        check_sums(seen[mg_si])
    assert seen[1.2].ground_type == "basaltic" and seen[0.7].ground_type == "intermediate" and seen[0.6].ground_type == "felsic"
    assert set(sc.GROUND_TYPES.row("intermediate").normalised_fractions()) <= set(sc.GROUND_TYPES.row("felsic").normalised_fractions()) | {"enstatite", "diopside"}
    assert sc.volcanic_ground_type(earth(rock_oxides=_oxides(0.9), internal_heat_flux_w_m2=0.09)) == "ultramafic"


def test_a_plagioclase_only_melt_reads_basaltic_never_anorthositic():
    """The stage melt of a calcium-poor stored rock is pure feldspar; an anorthosite is a flotation cumulate, not a melt of a mantle."""
    poor = rock_types.rock_oxides(0.08, 0.045, 1.2, 0.02)
    melt, _ = rock_types.partial_melt(poor, sc.melt_fraction(3e-4))
    assert rock_types.ground_type_of(melt) == "anorthositic"
    assert sc.volcanic_ground_type(earth(rock_oxides=poor, internal_heat_flux_w_m2=3e-4)) == "basaltic"


def test_the_primary_crust_and_the_bodies_without_a_stored_rock_keep_the_table_rules():
    base = dict(hydrosphere="dry", ocean_fraction=0.0, o2_partial_pa=0.0, mobile_lid=False)
    for seed in range(20):
        draw.set_run_seed(seed)
        assert sc.draw_surface(earth(internal_heat_flux_w_m2=1.0, **base)).ground_type == "basaltic"
        cold = sc.draw_surface(earth(internal_heat_flux_w_m2=0.0, rock_oxides=_oxides(0.7), rock_feo_fraction=0.18, **base)).ground_type
        assert cold in ("ultramafic", "anorthositic"), "no volcanism, no melt: the overturn rule decides"


def test_the_volcanic_crust_draws_nothing_the_table_rules_did_not():
    base = dict(hydrosphere="dry", ocean_fraction=0.0, o2_partial_pa=0.0, mobile_lid=True, internal_heat_flux_w_m2=0.09)
    after = []
    for oxides in (None, _oxides(1.2), _oxides(0.6)):
        draw.set_run_seed(5)
        sc.draw_surface(earth(rock_oxides=oxides, **base))
        after.append(draw.random())
    assert after[0] == after[1] == after[2]


def test_overturn_follows_feo():
    base = dict(hydrosphere="dry", ocean_fraction=0.0, o2_partial_pa=0.0, mobile_lid=False, internal_heat_flux_w_m2=0.0)
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
    base = dict(hydrosphere="dry", ocean_fraction=0.0, o2_partial_pa=0.0, mobile_lid=False, internal_heat_flux_w_m2=1.0)
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
        return sc.land_rock_ppm(types.SimpleNamespace(**{f"ground_{name}_fraction": 1.0}))

    ultramafic, basaltic, felsic = (reservoir(n) for n in ("ultramafic", "basaltic", "felsic"))
    for index in range(3):      # K, Th, U
        assert ultramafic[index] < basaltic[index] < felsic[index]


def _bodies(count=25):
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem

    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = False
    bodies = []
    for seed in range(count):
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        for planet in system.planets:
            if hasattr(planet, "moons"):
                bodies += [planet, *planet.moons]
    return bodies


def test_every_generated_body_has_a_consistent_ground():
    bodies = _bodies()
    solid = [b for b in bodies if b.body_type != "g"]
    assert solid and any(b.body_type == "g" for b in bodies)
    for body in bodies:
        if body.body_type == "g":
            assert all(getattr(body, field) is None for field in sc.COLUMN_FIELDS)
            continue
        mix = sc.mix_of(body)
        assert math.fsum(mix.values()) == pytest.approx(1.0, abs=1e-9)
        assert body.ground_type in sc.GROUND_FRACTION_TYPES and body.ground_grain in sc.GRAINS
        # agrees with the bulk fractions, the temperature and the water
        if mix.get(sc.ICE):
            assert body.surface_temperature < 273.15
        if mix.get(sc.OCEAN):
            assert body.hydrosphere in ("surface ocean", "brine seep", "hycean")
        if mix.get("carbon"):
            assert body.carbon_rock
        if mix.get("metal"):
            assert body.metal_mass_fraction >= sc.RANGES["metal_stripped"]
        assert "felsic" not in mix or body.hydrosphere is not None


def test_the_ground_dose_follows_the_stored_ground():
    felsic = types.SimpleNamespace(**{f"ground_{name}_fraction": 1.0 if name == "felsic" else 0.0 for name in sc.GROUND_FRACTION_TYPES})
    ultramafic = types.SimpleNamespace(**{f"ground_{name}_fraction": 1.0 if name == "ultramafic" else 0.0 for name in sc.GROUND_FRACTION_TYPES})
    none = types.SimpleNamespace()
    assert sc.land_rock_ppm(felsic)[1] > sc.land_rock_ppm(none)[1] > sc.land_rock_ppm(ultramafic)[1]
    assert sc.land_rock_ppm(none) == (sc.RESERVOIRS["basaltic"]["k_ppm"], sc.RESERVOIRS["basaltic"]["th_ppm"],
                                      sc.RESERVOIRS["basaltic"]["u_ppm"])
    ocean_only = types.SimpleNamespace(**{f"ground_{name}_fraction": 1.0 if name == "ocean" else 0.0 for name in sc.GROUND_FRACTION_TYPES})
    assert sc.land_rock_ppm(ocean_only) == sc.land_rock_ppm(none)
