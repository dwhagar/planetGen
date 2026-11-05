# tests/test_surface_liquid_mix.py

"""The makeup, the salt and the second liquid of the surface liquids (GEN.1218b, `physics/surface_liquid_world.py`)."""

import types

import pytest

from planetgen.physics import constants, surface_liquid_world as slw

from tests.air_support import fold_air
from tests.test_surface_liquid_world import _body

EARTH = constants.EARTH_MASS_TO_KG


def test_a_titan_like_sea_has_a_makeup_that_adds_up():
    titan = _body(ice_ch4_mass_fraction=0.01)
    slw.generate(titan)
    total = titan.surface_liquid_methane_mol + titan.surface_liquid_ethane_mol + titan.surface_liquid_nitrogen_mol
    assert titan.surface_liquid == "methane" and total == pytest.approx(1.0, abs=0.05)
    assert titan.surface_liquid_ethane_mol > 0.2 and 0.1 < titan.surface_liquid_methane_mol < 0.9 and 0.0 <= titan.surface_liquid_nitrogen_mol < 0.25


def test_a_sea_supplied_by_ice_with_none_of_its_gas_in_the_air_is_that_liquid_alone():
    sea = _body(ice_ch4_mass_fraction=0.01, p_ch4_kpa=0.0)
    slw.generate(sea)
    assert sea.surface_liquid == "methane" and sea.surface_liquid_methane_mol == 1.0 and sea.surface_liquid_ethane_mol == 0.0


def test_a_water_sea_has_the_salt_of_its_water_activity_and_a_dry_body_none():
    ocean = _body(mass=EARTH, radius=6371.0, gravity=1.0, surface_temperature=288.0, atmospheric_pressure=1.013e5, p_n2_kpa=78.0, p_o2_kpa=21.0,
                  hydrosphere="surface ocean", ocean_fraction=0.7, ocean_depth_km=3.7, water_activity=0.981)
    slw.generate(ocean)
    assert ocean.surface_liquid_salt_wt == pytest.approx(3.5, rel=0.2)       # seawater: 35 g/kg (3.5 percent) at a_w 0.981
    assert ocean.surface_liquid_methane_mol is None
    ocean.water_activity = 1.0
    slw.generate(ocean)
    assert ocean.surface_liquid_salt_wt == 0.0
    ocean.water_activity = 0.5
    slw.generate(ocean)
    assert ocean.surface_liquid_salt_wt is None       # a brine of other salts has no sodium chloride equivalent
    dry = _body(hydrosphere="dry")
    slw.generate(dry)
    assert dry.surface_liquid_salt_wt is None


def test_a_body_with_two_stable_supplied_liquids_keeps_the_second():
    body = _body(ice_ch4_mass_fraction=0.01, ice_other_mass_fraction=0.01, surface_temperature=120.0, atmospheric_pressure=3.0e5, p_n2_kpa=250.0, p_ch4_kpa=1.0)
    slw.generate(body)
    if body.surface_liquid_2 is not None:
        assert body.surface_liquid_2 != body.surface_liquid
        assert 0.0 < body.surface_liquid_2_fraction <= 1.0 - body.surface_liquid_fraction + 1e-12


def test_a_lava_surface_and_a_single_liquid_have_no_second():
    hot = _body(mass=EARTH, radius=6371.0, gravity=1.0, surface_temperature=1600.0, atmospheric_pressure=0.0, p_n2_kpa=0.0, p_ch4_kpa=0.0)
    slw.generate(hot)
    assert hot.surface_liquid_2 is None and hot.surface_liquid_2_fraction is None
    titan = _body(ice_ch4_mass_fraction=0.01)
    slw.generate(titan)
    assert titan.surface_liquid_2 is None


def test_every_new_field_is_cleared_for_a_giant():
    giant = types.SimpleNamespace(body_type="g", **{field: 1.0 for field in slw.COLUMN_FIELDS})
    slw.generate(giant)
    assert all(getattr(giant, field) is None for field in slw.COLUMN_FIELDS)
