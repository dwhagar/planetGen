"""
What an ocean takes out of a drawn mixture (GEN.757, `planetgen.physics.mixture_sinks`): the sulfur uptake and the nitrogen uptake (GEN.1834). The Earth is a test target, never an input of the code.
"""

import pytest

from planetgen.physics import gases, mixture_sinks as sinks

G_EARTH = 9.81
EARTH_OCEAN_KG_M2 = 2.7e6
"""float: The global mean of the Earth's ocean, 1.4e21 kg over 5.1e14 m2."""
AIR = {"N2": 0.78, "O2": 0.2, "H2O": 0.03, "SO2": 1.0e-4, "H2S": 1.0e-4, "CO2": 0.0198}


def _air():
    total = sum(AIR.values())
    return {formula: value / total for formula, value in AIR.items()}


def test_an_earth_ocean_keeps_a_few_percent_of_the_sulfur_dioxide_and_about_half_the_hydrogen_sulfide():
    so2 = sinks.remaining_fraction("SO2", 288.0, G_EARTH, EARTH_OCEAN_KG_M2)
    h2s = sinks.remaining_fraction("H2S", 288.0, G_EARTH, EARTH_OCEAN_KG_M2)
    assert 0.01 < so2 < 0.1
    assert 0.3 < h2s < 0.7
    assert so2 < h2s


def test_an_earth_ocean_keeps_the_nitrogen_it_takes_under_a_percent_dissolved_and_the_air_the_rest():
    """GEN.1834: 1.4e21 kg of water at 288 K and 0.78 bar of N2 hold 2.2e16 kg of N2, 0.56 percent of the air's 3.9e18 kg."""
    kept = sinks.remaining_fraction("N2", 288.0, G_EARTH, EARTH_OCEAN_KG_M2)
    assert 0.99 < kept < 0.997
    assert sinks.remaining_fraction("N2", 288.0, G_EARTH, 10.0 * EARTH_OCEAN_KG_M2) == pytest.approx(1.0 / (1.0 + 10.0 * (1.0 / kept - 1.0)), rel=1e-12)
    assert kept > sinks.remaining_fraction("SO2", 288.0, G_EARTH, EARTH_OCEAN_KG_M2) > 0.0


def test_the_uptake_follows_henrys_law_in_the_formula():
    kh = gases.gas("SO2").henry.kh(288.0)
    expected = 1.0 / (1.0 + kh * EARTH_OCEAN_KG_M2 * G_EARTH * gases.gas("SO2").molar_mass / 1.0e5)
    assert sinks.remaining_fraction("SO2", 288.0, G_EARTH, EARTH_OCEAN_KG_M2) == pytest.approx(expected, rel=1e-12)


def test_more_ocean_or_more_gravity_or_a_colder_ocean_takes_more():
    base = sinks.remaining_fraction("SO2", 288.0, G_EARTH, EARTH_OCEAN_KG_M2)
    assert sinks.remaining_fraction("SO2", 288.0, G_EARTH, 2.0 * EARTH_OCEAN_KG_M2) < base
    assert sinks.remaining_fraction("SO2", 288.0, 2.0 * G_EARTH, EARTH_OCEAN_KG_M2) < base
    assert sinks.remaining_fraction("SO2", 275.0, G_EARTH, EARTH_OCEAN_KG_M2) < base
    assert sinks.remaining_fraction("SO2", 288.0, G_EARTH, 0.0) == 1.0


def test_the_henry_law_range_clamps_instead_of_failing():
    assert sinks.remaining_fraction("SO2", 400.0, G_EARTH, EARTH_OCEAN_KG_M2) == pytest.approx(
        sinks.remaining_fraction("SO2", gases.HENRY_RANGE_K[1], G_EARTH, EARTH_OCEAN_KG_M2))
    assert sinks.remaining_fraction("SO2", 250.0, G_EARTH, EARTH_OCEAN_KG_M2) == pytest.approx(
        sinks.remaining_fraction("SO2", gases.HENRY_RANGE_K[0], G_EARTH, EARTH_OCEAN_KG_M2))


def test_nitrogen_stays_whole_outside_the_henry_range_where_the_sulfur_gases_take_the_nearest_end():
    """GEN.1834: a frozen sea is cut off from the air and a hot one is past the data; the sulfur gases keep the clamp of GEN.757."""
    low, high = gases.HENRY_RANGE_K
    deep = 1000.0 * EARTH_OCEAN_KG_M2
    assert sinks.remaining_fraction("N2", low - 1.0, G_EARTH, deep) == 1.0
    assert sinks.remaining_fraction("N2", high + 1.0, G_EARTH, deep) == 1.0
    assert sinks.remaining_fraction("N2", low, G_EARTH, deep) < 0.5
    assert sinks.remaining_fraction("N2", high, G_EARTH, deep) < 0.5
    assert sinks.remaining_fraction("SO2", low - 1.0, G_EARTH, deep) < 1.0


def test_bad_input_is_refused():
    with pytest.raises(ValueError):
        sinks.remaining_fraction("CO2", 288.0, G_EARTH, 1.0e6)
    with pytest.raises(ValueError):
        sinks.remaining_fraction("SO2", 288.0, 0.0, 1.0e6)
