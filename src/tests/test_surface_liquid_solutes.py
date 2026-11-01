# tests/test_surface_liquid_solutes.py

"""What the surface liquids hold (GEN.1215, `physics/surface_liquid_world.py`): the strength and place of the liquid sulfuric acid the air makes, the
ammonia of a liquid water ocean, and the life classes reading the stored values in place of working them out again."""

import types

import pytest

from planetgen.physics import acid_medium, constants, strange_solvents, surface_liquid_world as slw

from tests.air_support import fold_air

EARTH = constants.EARTH_MASS_TO_KG


def _body(**extra):
    fields = dict(body_type="t", mass=EARTH, radius=6371.0, gravity=1.0, surface_temperature=330.0, atmospheric_pressure=1.0e5, hydrosphere="dry",
                  p_co2_kpa=80.0, p_n2_kpa=19.0, p_so2_kpa=1.0, p_h2o_kpa=0.002)
    fields.update(extra)
    return types.SimpleNamespace(**fold_air(fields))


def test_a_warm_dry_sulfurous_ground_stores_a_strong_acid_pool():
    body = _body()
    slw.generate(body)
    assert body.surface_liquid_acid_layer == "pool"
    assert body.surface_liquid_acid_wt >= acid_medium.number("concentrated_wt_percent")


def test_air_without_the_supply_stores_no_acid():
    body = _body(p_so2_kpa=0.0)
    slw.generate(body)
    assert body.surface_liquid_acid_wt is None and body.surface_liquid_acid_layer is None


def test_an_ocean_with_other_ice_stores_its_ammonia_as_a_percent_of_the_liquid():
    body = _body(p_so2_kpa=0.0, hydrosphere="ice-covered ocean", ocean_fraction=0.9, water_mass_fraction=0.1, ice_mass_fraction=0.02)
    slw.generate(body)
    share = strange_solvents.number("nh3_share_of_ice") * 0.02 / 0.1
    assert body.surface_liquid_ammonia_wt == pytest.approx(100.0 * share / (1.0 + share))


@pytest.mark.parametrize("state", ["dry", "ice", "vapour"])
def test_a_body_without_a_liquid_ocean_stores_no_ammonia(state):
    body = _body(p_so2_kpa=0.0, hydrosphere=state, water_mass_fraction=0.1, ice_mass_fraction=0.02)
    slw.generate(body)
    assert body.surface_liquid_ammonia_wt is None


def test_a_giant_stores_neither():
    body = _body(body_type="g")
    slw.generate(body)
    assert all(getattr(body, name) is None for name in slw.SOLUTE_FIELDS)


def test_the_life_classes_read_the_stored_values_and_work_them_out_when_there_are_none():
    assert slw.acid_level({"surface_liquid_acid_layer": "pool"}) == acid_medium.POOL
    assert slw.acid_level({"surface_liquid_acid_layer": "cloud"}) == acid_medium.CLOUD
    values = {"surface_temperature_k": None}
    assert slw.acid_level(values) == acid_medium.level(values)
    share = 0.37
    assert slw.ammonia_water_mass_fraction({"surface_liquid_ammonia_wt": 100.0 * share / (1.0 + share)}) == pytest.approx(share)
    values = {"ice_mass_fraction": 0.02, "water_mass_fraction": 0.1}
    assert slw.ammonia_water_mass_fraction(values) == strange_solvents.ammonia_water_mass_fraction(values)
