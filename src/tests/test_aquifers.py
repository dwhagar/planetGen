# tests/test_aquifers.py

"""
GEN.654: a world with no liquid at its surface can hold water in the pores of its crust. The layer follows the geotherm and the pore
closure, the water fills it from the bottom, and it changes the habitability score in four places (solvent, dose, energy, Water rung).
"""

import pytest

from planetgen import tuning
from planetgen.generation import life as planet_life
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import aquifers, habitability_world

from tests.fuzz_support import deterministic_entropy

MARS = dict(surface_temperature_k=210.0, surface_pressure_pa=600.0, gravity_m_s2=3.71, flux_w_m2=0.03)


def test_the_top_of_the_layer_is_the_273_k_level_when_there_is_water_to_fill_the_pores():
    depth, thickness = aquifers.layer_km(**MARS, water_depth_m=5000.0)
    gradient = 0.03 / tuning.AQUIFER_CONDUCTIVITY_W_M_K * 1000.0
    assert depth == pytest.approx((273.15 - 210.0) / gradient, rel=1e-6)
    assert 1.0 < thickness < 10.0


def test_a_smaller_inventory_leaves_a_deeper_and_thinner_layer():
    layers = [aquifers.layer_km(**MARS, water_depth_m=water) for water in (5000.0, 50.0, 5.0)]
    assert all(layer[0] is not None for layer in layers)
    depths, thicknesses = zip(*layers)
    assert list(depths) == sorted(depths) and list(thicknesses) == sorted(thicknesses, reverse=True)
    assert depths[0] < depths[1] < depths[2]


def test_no_water_no_heat_or_a_surface_too_hot_for_life_is_no_layer():
    assert aquifers.layer_km(**MARS, water_depth_m=0.0) == (None, None)
    assert aquifers.layer_km(**{**MARS, "flux_w_m2": 0.0}, water_depth_m=100.0) == (None, None)
    assert aquifers.layer_km(**{**MARS, "surface_temperature_k": tuning.AQUIFER_LIFE_MAX_K + 1.0}, water_depth_m=100.0) == (None, None)
    assert aquifers.layer_km(**MARS, water_depth_m=1e-4) == (None, None)  # under the minimum thickness


def test_a_warm_dry_surface_boils_the_first_metres_away_and_a_lighter_world_has_deeper_pores():
    depth, _thickness = aquifers.layer_km(380.0, 1.0e5, 9.8, 0.08, 5000.0)
    assert 0.0 < depth < 0.05
    earth = aquifers.layer_km(240.0, 600.0, 9.8, 0.06, 5000.0)
    moon = aquifers.layer_km(240.0, 600.0, 1.6, 0.06, 5000.0)
    assert moon[1] > earth[1]
    assert aquifers.layer_km(210.0, 600.0, 9.8, 0.03, 5000.0) == (None, None)  # the pores close above the 273 K level


def _bodies(count):
    for seed in range(count):
        config = SystemConfig()
        config.MOONS = True
        with deterministic_entropy(9100 + seed):
            system = StarSystem(system_config=config)
        for planet in system.planets + system.secondary_planets:
            if planet.body_type == "a":
                continue
            yield planet
            yield from planet.moons


def test_only_a_world_with_no_liquid_at_the_surface_has_an_aquifer():
    with_layer = 0
    for body in _bodies(30):
        if body.body_type == "g":
            assert body.aquifer_depth_km is None and body.aquifer_thickness_km is None
            continue
        if body.aquifer_thickness_km is None:
            assert body.aquifer_depth_km is None
            continue
        assert body.hydrosphere in aquifers.SURFACE_STATES
        assert body.aquifer_thickness_km >= tuning.AQUIFER_MIN_THICKNESS_KM and body.aquifer_depth_km >= 0.0
        with_layer += 1
    assert with_layer > 5


def _a_world_with_a_reachable_layer():
    for body in _bodies(80):
        if body.body_type != "g" and body.aquifer_thickness_km and body.aquifer_depth_km <= tuning.AQUIFER_REACH_KM:
            return body
    pytest.fail("no generated world with a reachable aquifer in 80 seeds")


def test_a_reachable_aquifer_is_a_yellow_water_rung_and_lets_microbes_score_but_not_flora():
    body = _a_world_with_a_reachable_layer()
    values = habitability_world.values_from_body(body)
    assert habitability_world.groundwater(values) == "reachable"
    assert habitability_world.solvent(values) == "water"
    assert habitability_world.water_source(values) in ("aquifer", "ice")
    assert body.tier_water == 2  # Yellow
    assert body.life_flora == body.life_sentient == 0
    assert body.life_fauna <= body.life_microbial   # fauna comes from the class's own producer base, not from flora (docs/design/habitability/life-classes.md section 7): a Class 3 aquifer world can hold animal grade A1
    without = habitability_world.compute({**values, "aquifer_depth_km": None, "aquifer_thickness_km": None})
    assert without["l_solv"] == 0.0 and without["tier_water"] == 3  # Red: no liquid at all
    with_layer = habitability_world.compute(values)
    assert with_layer["l_solv"] == 1.0 and with_layer["phi4_water"] == pytest.approx(0.6)


def test_a_deep_aquifer_shelters_microbes_but_the_water_is_out_of_reach():
    body = _a_world_with_a_reachable_layer()
    values = {**habitability_world.values_from_body(body), "aquifer_depth_km": tuning.AQUIFER_REACH_KM + 1.0}
    assert habitability_world.groundwater(values) == "deep"
    scores = habitability_world.compute(values)
    assert scores["l_solv"] == 1.0 and scores["tier_water"] == 3
