# tests/test_photochemical_loss.py

"""GEN.1732: the star's ultraviolet takes the ammonia and methane out of an air (`physics/photochemical_loss.py`)."""

import math

import pytest

from planetgen.physics import aerosols, gases, photochemical_loss as pl, photolysis, unit_registry

SUN_AT_1_AU = aerosols.FarUv(1.0, 0.0)
SUN_T_K = 5772.0
GYR_S = 1.0e9 * unit_registry.unit("yr").factor
EARTH_G = 9.80665


def test_a_source_against_a_loss_stands_at_the_source_over_the_rate_and_never_over_the_inventory():
    rate = 1.0e-4
    assert pl.standing(2.0e-6, rate) == pytest.approx(2.0e-2)
    assert pl.standing(0.0, rate) == 0.0 and pl.standing(1.0, 0.0) == math.inf
    for inventory, age_s in ((10.0, 1.0e9), (1.0e3, 1.0e15), (5.0, 1.0)):
        for loss in (1.0e-12, 1.0e-9, 1.0e-4):
            left = pl.standing_from_inventory(inventory, age_s, loss)
            assert 0.0 <= left <= inventory * (1.0 + 1e-12)
    assert pl.standing_from_inventory(10.0, 1.0e12, 1.0e-3) == pytest.approx(10.0 / 1.0e12 / 1.0e-3, rel=1e-9)       # a short life: the source over the rate
    assert pl.standing_from_inventory(10.0, 1.0e3, 1.0e-12) == pytest.approx(10.0, rel=1e-6)                         # a long life: nearly all of it
    assert pl.standing_from_inventory(10.0, 0.0, 1.0) == 10.0 and pl.standing_from_inventory(10.0, 5.0, 0.0) == 10.0


def test_the_exact_step_decays_with_no_source_fills_to_the_standing_amount_and_does_not_depend_on_the_step():
    rate = 3.0e-9
    assert pl.advance(8.0, 0.0, rate, 1.0 / rate) == pytest.approx(8.0 / math.e, rel=1e-12)
    assert pl.advance(0.0, 5.0e-9, rate, 100.0 / rate) == pytest.approx(pl.standing(5.0e-9, rate), rel=1e-9)
    once = pl.advance(3.0, 2.0e-9, rate, 4.0e8)
    twice = pl.advance(pl.advance(3.0, 2.0e-9, rate, 2.0e8), 2.0e-9, rate, 2.0e8)
    assert once == pytest.approx(twice, rel=1e-12)
    assert pl.advance(3.0, 2.0, 0.0, 10.0) == pytest.approx(23.0) and pl.advance(3.0, 2.0, rate, 0.0) == 3.0


def test_a_source_over_what_the_light_can_break_gains_the_difference_and_the_law_is_continuous_there():
    rate, age = 1.0e-3, 1.0e6
    limit = 2.0e-6
    below = pl.standing_from_inventory(limit * age * (1.0 - 1e-9), age, rate, limit)
    above = pl.standing_from_inventory(limit * age * (1.0 + 1e-9), age, rate, limit)
    assert above == pytest.approx(below, rel=1e-5)
    more = pl.standing_from_inventory(4.0 * limit * age, age, rate, limit)
    assert more > limit / rate and more < 4.0 * limit * age
    assert more == pytest.approx(limit / rate + (4.0 * limit - limit) * (age + math.log1p(-0.25) / rate), rel=1e-9)
    assert pl.advance(0.0, 4.0 * limit, rate, age, limit) == pytest.approx(3.0 * limit * age, rel=1e-9)


def test_the_ammonia_of_a_warm_air_is_gone_under_the_quiet_sun_and_the_nitrogen_it_leaves_is_half_a_mole_a_mole():
    assert photolysis.lifetime_s("NH3") == pytest.approx(6.7e3, rel=0.01)
    air = {"N2": 78.0e3, "NH3": 1000.0, "CH4": 200.0, "H2O": 1500.0, "CO2": 30.0}
    after = pl.lose(air, EARTH_G, 4.5 * GYR_S, SUN_AT_1_AU, SUN_T_K)
    assert after["NH3"] < 1e-6 * air["NH3"] and after["CH4"] < 1e-3 * air["CH4"]
    assert after["H2O"] == air["H2O"] and after["N2"] > air["N2"] and after["CO2"] > air["CO2"]
    moles = lambda gas, pa: pa / gases.gas(gas).molar_mass
    nitrogen_before = 2.0 * moles("N2", air["N2"]) + moles("NH3", air["NH3"])
    nitrogen_after = 2.0 * moles("N2", after["N2"]) + moles("NH3", after["NH3"])
    assert nitrogen_after == pytest.approx(nitrogen_before, rel=1e-9)
    carbon_before = moles("CO2", air["CO2"]) + moles("CH4", air["CH4"])
    carbon_after = moles("CO2", after["CO2"]) + moles("CH4", after["CH4"])
    assert carbon_after == pytest.approx(carbon_before, rel=1e-9)


def test_a_young_body_or_a_dim_star_keeps_more_and_no_far_ultraviolet_keeps_all():
    air = {"N2": 1.0e5, "NH3": 500.0}
    old = pl.lose(air, EARTH_G, 4.5 * GYR_S, SUN_AT_1_AU, SUN_T_K)["NH3"]
    young = pl.lose(air, EARTH_G, 1.0e-6 * GYR_S, SUN_AT_1_AU, SUN_T_K)["NH3"]
    far = pl.lose(air, EARTH_G, 4.5 * GYR_S, aerosols.FarUv(1.0e-4, 0.0), SUN_T_K)["NH3"]
    assert young > far > old
    assert pl.lose(air, EARTH_G, 4.5 * GYR_S, None, SUN_T_K) == air
    assert pl.lose({"N2": 1e5}, EARTH_G, 4.5 * GYR_S, SUN_AT_1_AU, SUN_T_K) == {"N2": 1e5}
    assert pl.lose(air, EARTH_G, 4.5 * GYR_S, SUN_AT_1_AU, None) == air


def test_the_columns_and_pressures_go_both_ways():
    column = pl.column_mol_m2(1000.0, "NH3", EARTH_G)
    assert column == pytest.approx(1000.0 / (EARTH_G * gases.gas("NH3").molar_mass))
    assert pl.partial_pa_of(column, "NH3", EARTH_G) == pytest.approx(1000.0)
    assert pl.light_limit("CH4", None, SUN_T_K) == math.inf
    assert 0.0 < pl.light_limit("CH4", SUN_AT_1_AU, SUN_T_K) < math.inf
