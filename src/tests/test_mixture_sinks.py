"""
What an ocean takes out of a drawn mixture (GEN.757, `planetgen.physics.mixture_sinks`): the water vapour cap and the sulfur
uptake. The Earth is a test target, never an input of the code.
"""

import math

import pytest

from planetgen.physics import gases, mixture_sinks as sinks

G_EARTH = 9.81
EARTH_OCEAN_KG_M2 = 2.7e6
"""float: The global mean of the Earth's ocean, 1.4e21 kg over 5.1e14 m2."""
AIR = {"N2": 0.78, "O2": 0.2, "H2O": 0.03, "SO2": 1.0e-4, "H2S": 1.0e-4, "CO2": 0.0198}


def _air():
    total = sum(AIR.values())
    return {formula: value / total for formula, value in AIR.items()}


def test_the_water_over_an_earth_ocean_is_held_at_saturation():
    sunk = sinks.apply_ocean_sinks(_air(), 1.0e5, 288.0, G_EARTH, ocean_present=True)
    assert sunk.mole_fractions["H2O"] * sunk.pressure_pa == pytest.approx(sinks.water_vapour_cap_pa(288.0), rel=1e-9)
    assert sinks.water_vapour_cap_pa(288.0) == pytest.approx(1.69e3, rel=0.01)
    assert sunk.condensed_water_pa == pytest.approx(_air()["H2O"] * 1.0e5 - sinks.water_vapour_cap_pa(288.0), rel=1e-9)
    assert sum(sunk.mole_fractions.values()) == pytest.approx(1.0)
    assert sunk.pressure_pa == pytest.approx(1.0e5 - sunk.condensed_water_pa)


def test_a_drier_air_than_saturation_is_left_alone():
    air = {"N2": 0.9899, "H2O": 0.01, "CO2": 0.0001}
    sunk = sinks.apply_ocean_sinks(air, 1.0e5, 288.0, G_EARTH, ocean_present=True)
    assert sunk.condensed_water_pa == 0.0
    assert sunk.pressure_pa == pytest.approx(1.0e5)
    for formula, fraction in air.items():
        assert sunk.mole_fractions[formula] == pytest.approx(fraction)


def test_without_an_ocean_nothing_is_taken():
    sunk = sinks.apply_ocean_sinks(_air(), 1.0e5, 288.0, G_EARTH, ocean_present=False, ocean_column_kg_m2=EARTH_OCEAN_KG_M2)
    assert sunk.pressure_pa == pytest.approx(1.0e5)
    assert sunk.condensed_water_pa == 0.0 and sunk.dissolved_pa == {}
    assert sinks.water_vapour_cap_pa(288.0, ocean_present=False) == math.inf


def test_above_the_critical_point_water_is_not_capped():
    assert sinks.water_vapour_cap_pa(700.0) == math.inf
    air = {"H2O": 0.9, "CO2": 0.1}
    sunk = sinks.apply_ocean_sinks(air, 1.0e7, 700.0, 9.0, ocean_present=True)
    assert sunk.pressure_pa == pytest.approx(1.0e7) and sunk.condensed_water_pa == 0.0


def test_the_cap_rises_with_temperature():
    caps = [sinks.water_vapour_cap_pa(t) for t in (280.0, 300.0, 330.0, 400.0, 500.0)]
    assert caps == sorted(caps)


def test_an_earth_ocean_keeps_a_few_percent_of_the_sulfur_dioxide_and_about_half_the_hydrogen_sulfide():
    so2 = sinks.sulfur_remaining_fraction("SO2", 288.0, G_EARTH, EARTH_OCEAN_KG_M2)
    h2s = sinks.sulfur_remaining_fraction("H2S", 288.0, G_EARTH, EARTH_OCEAN_KG_M2)
    assert 0.01 < so2 < 0.1
    assert 0.3 < h2s < 0.7
    assert so2 < h2s


def test_the_uptake_follows_henrys_law_in_the_formula():
    kh = gases.gas("SO2").henry.kh(288.0)
    expected = 1.0 / (1.0 + kh * EARTH_OCEAN_KG_M2 * G_EARTH * gases.gas("SO2").molar_mass / 1.0e5)
    assert sinks.sulfur_remaining_fraction("SO2", 288.0, G_EARTH, EARTH_OCEAN_KG_M2) == pytest.approx(expected, rel=1e-12)


def test_more_ocean_or_more_gravity_or_a_colder_ocean_takes_more():
    base = sinks.sulfur_remaining_fraction("SO2", 288.0, G_EARTH, EARTH_OCEAN_KG_M2)
    assert sinks.sulfur_remaining_fraction("SO2", 288.0, G_EARTH, 2.0 * EARTH_OCEAN_KG_M2) < base
    assert sinks.sulfur_remaining_fraction("SO2", 288.0, 2.0 * G_EARTH, EARTH_OCEAN_KG_M2) < base
    assert sinks.sulfur_remaining_fraction("SO2", 275.0, G_EARTH, EARTH_OCEAN_KG_M2) < base
    assert sinks.sulfur_remaining_fraction("SO2", 288.0, G_EARTH, 0.0) == 1.0


def test_the_henry_law_range_clamps_instead_of_failing():
    assert sinks.sulfur_remaining_fraction("SO2", 400.0, G_EARTH, EARTH_OCEAN_KG_M2) == pytest.approx(
        sinks.sulfur_remaining_fraction("SO2", gases.HENRY_RANGE_K[1], G_EARTH, EARTH_OCEAN_KG_M2))
    assert sinks.sulfur_remaining_fraction("SO2", 250.0, G_EARTH, EARTH_OCEAN_KG_M2) == pytest.approx(
        sinks.sulfur_remaining_fraction("SO2", gases.HENRY_RANGE_K[0], G_EARTH, EARTH_OCEAN_KG_M2))


def test_the_sulfur_taken_and_the_sulfur_left_add_up_and_the_other_gases_keep_their_pressure():
    air, pressure = _air(), 1.0e5
    sunk = sinks.apply_ocean_sinks(air, pressure, 288.0, G_EARTH, ocean_present=True, ocean_column_kg_m2=EARTH_OCEAN_KG_M2)
    for formula in sinks.SULFUR_GASES:
        left = sunk.mole_fractions[formula] * sunk.pressure_pa
        assert left + sunk.dissolved_pa[formula] == pytest.approx(air[formula] * pressure, rel=1e-9)
        assert sunk.dissolved_pa[formula] > 0.0
    for formula in ("N2", "O2", "CO2"):
        assert sunk.mole_fractions[formula] * sunk.pressure_pa == pytest.approx(air[formula] * pressure, rel=1e-9)


def test_an_air_of_one_gas_with_nothing_to_take_is_returned_whole():
    sunk = sinks.apply_ocean_sinks({"N2": 1.0}, 5.0e4, 288.0, G_EARTH, ocean_present=True, ocean_column_kg_m2=EARTH_OCEAN_KG_M2)
    assert sunk.mole_fractions == {"N2": pytest.approx(1.0)} and sunk.pressure_pa == pytest.approx(5.0e4)


def test_an_air_that_the_ocean_takes_entirely_is_empty():
    sunk = sinks.apply_ocean_sinks({"H2O": 1.0}, 1.0e5, 280.0, G_EARTH, ocean_present=True)
    assert sunk.pressure_pa == pytest.approx(sinks.water_vapour_cap_pa(280.0))
    assert sunk.mole_fractions == {"H2O": pytest.approx(1.0)}


def test_bad_input_is_refused():
    with pytest.raises(ValueError):
        sinks.sulfur_remaining_fraction("CO2", 288.0, G_EARTH, 1.0e6)
    with pytest.raises(ValueError):
        sinks.sulfur_remaining_fraction("SO2", 288.0, 0.0, 1.0e6)
    for bad in ((0.0, 288.0, G_EARTH), (1.0e5, 0.0, G_EARTH), (1.0e5, 288.0, -1.0)):
        with pytest.raises(ValueError):
            sinks.apply_ocean_sinks({"N2": 1.0}, *bad, ocean_present=True)
