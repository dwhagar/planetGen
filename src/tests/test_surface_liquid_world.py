# tests/test_surface_liquid_world.py

"""The liquid a world holds on its surface (GEN.1035 wiring, `physics/surface_liquid_world.py`): water of the hydrosphere, lava of a hot surface, a Titan-like
methane sea from the ice, a sea condensed from the air, and the worlds that hold none."""

import types

import pytest

from planetgen.physics import constants, surface_composition, surface_liquid_world as slw

from tests.air_support import fold_air

EARTH = constants.EARTH_MASS_TO_KG


def _body(**extra):
    """A Titan-like body (a small moon with a thick nitrogen air at 94 K) unless `extra` changes it."""
    fields = dict(body_type="t", mass=0.0225 * EARTH, radius=2575.0, gravity=0.138, surface_temperature=94.0, atmospheric_pressure=1.47e5,
                  hydrosphere="dry", p_n2_kpa=140.0, p_ch4_kpa=5.6)
    fields.update(extra)
    return types.SimpleNamespace(**fold_air(fields))


def test_a_titan_like_body_holds_a_methane_sea_from_its_ice():
    titan = _body(ice_ch4_mass_fraction=0.01)
    slw.generate(titan)
    assert titan.surface_liquid == "methane"
    assert titan.surface_liquid_fraction == pytest.approx(0.017, abs=0.004)
    assert 0.1 < titan.surface_liquid_depth_km < 3.0


def test_an_air_over_saturation_condenses_into_a_sea():
    sea = _body(surface_temperature=100.0, atmospheric_pressure=2.0e5, p_ch4_kpa=60.0)
    slw.generate(sea)
    assert sea.surface_liquid == "methane" and 0.0 < sea.surface_liquid_fraction <= 1.0


def test_a_water_world_has_water_over_the_cover_the_hydrosphere_stored():
    earth = _body(mass=EARTH, radius=6371.0, gravity=1.0, surface_temperature=288.0, atmospheric_pressure=1.013e5, p_n2_kpa=78.0, p_o2_kpa=21.0,
                  hydrosphere="surface ocean", ocean_fraction=0.7, ocean_depth_km=3.7)
    slw.generate(earth)
    assert (earth.surface_liquid, earth.surface_liquid_fraction, earth.surface_liquid_depth_km) == ("water", 0.7, 3.7)


def test_a_frozen_or_dry_water_world_has_no_water_surface_liquid():
    dry = _body(mass=EARTH, radius=6371.0, gravity=1.0, surface_temperature=288.0, atmospheric_pressure=1.013e5, p_n2_kpa=100.0, hydrosphere="dry")
    slw.generate(dry)
    assert dry.surface_liquid != "water"
    ice = _body(mass=EARTH, radius=6371.0, gravity=1.0, surface_temperature=200.0, hydrosphere="ice-covered ocean", ocean_fraction=0.9)
    slw.generate(ice)
    assert ice.surface_liquid != "water"


@pytest.mark.parametrize("ground,lava,kelvin", [("felsic", "rhyolite_melt", 1400.0), ("ultramafic", "komatiite_melt", 2000.0), ("basaltic", "basalt_melt", 1600.0)])
def test_a_hot_surface_is_a_lava_of_its_ground(ground, lava, kelvin):
    hot = _body(mass=EARTH, radius=6371.0, gravity=1.0, surface_temperature=kelvin, atmospheric_pressure=0.0, p_n2_kpa=0.0, p_ch4_kpa=0.0,
                **{f"ground_{ground}_fraction": 0.8})
    slw.generate(hot)
    assert hot.surface_liquid == lava and hot.surface_liquid_fraction == 1.0 and hot.surface_liquid_depth_km is None


def test_a_surface_hotter_than_the_peridotite_liquidus_is_a_magma_ocean():
    hot = _body(mass=EARTH, radius=6371.0, gravity=1.0, surface_temperature=2400.0, atmospheric_pressure=0.0, p_n2_kpa=0.0, p_ch4_kpa=0.0)
    slw.generate(hot)
    assert hot.surface_liquid == "peridotite_magma"


def test_airless_cold_and_giant_bodies_hold_none_and_clear_what_was_set():
    cold = _body(surface_temperature=40.0, atmospheric_pressure=0.0, p_n2_kpa=0.0, p_ch4_kpa=0.0, surface_liquid="water", surface_liquid_fraction=1.0)
    slw.generate(cold)
    assert all(getattr(cold, field) is None for field in slw.COLUMN_FIELDS)
    giant = _body(body_type="g", surface_liquid="water")
    slw.generate(giant)
    assert all(getattr(giant, field) is None for field in slw.COLUMN_FIELDS)


def test_the_columns_are_the_surface_compositions_so_the_store_saves_them():
    assert set(slw.COLUMN_FIELDS) <= set(surface_composition.COLUMN_FIELDS)


def test_a_cooler_surface_than_its_grounds_liquidus_is_solid():
    cool = _body(mass=EARTH, radius=6371.0, gravity=1.0, surface_temperature=1500.0, atmospheric_pressure=0.0, p_n2_kpa=0.0, p_ch4_kpa=0.0,
                 ground_ultramafic_fraction=0.8)
    slw.generate(cool)
    assert cool.surface_liquid is None


def test_an_ammonia_and_water_ice_body_holds_the_ammonia_water_eutectic_where_pure_ammonia_would_be_frozen():
    """GEN.1218a: at 200 K pure ammonia is a solid (195.5 K melting is under it only just, and it boils at 240 K), but the 32 percent solution is liquid from 176 K."""
    body = _body(surface_temperature=190.0, atmospheric_pressure=1.0e4, p_n2_kpa=10.0, p_ch4_kpa=0.0, ice_nh3_mass_fraction=0.01, water_mass_fraction=0.05)
    slw.generate(body)
    assert body.surface_liquid == "ammonia_water"
    assert body.surface_liquid_ammonia_wt == slw.AMMONIA_WATER_MASS_PERCENT
    dry = _body(surface_temperature=190.0, atmospheric_pressure=1.0e4, p_n2_kpa=10.0, p_ch4_kpa=0.0, ice_nh3_mass_fraction=0.01, water_mass_fraction=0.0)
    slw.generate(dry)
    assert dry.surface_liquid is None            # ammonia with no water makes no solution, and pure ammonia is frozen at 190 K


def test_methane_is_the_sea_where_methane_and_carbon_monoxide_are_both_liquid_and_supplied():
    """GEN.1218a, Boss's answer 6: the larger carbon monoxide ice does not take the sea from methane where both are liquid."""
    from planetgen.physics import surface_liquids
    both = surface_liquids.candidates(105.0, 5.0e6, {}, {"carbon_monoxide": 0.05, "methane": 0.01})
    assert [c.liquid_id for c in both][:2] == ["methane", "carbon_monoxide"]
    only_co = surface_liquids.candidates(75.0, 1.0e5, {}, {"carbon_monoxide": 0.05, "methane": 0.01})
    assert [c.liquid_id for c in only_co] == ["carbon_monoxide"]          # methane is a solid at 75 K
