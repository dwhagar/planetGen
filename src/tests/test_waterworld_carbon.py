"""
The carbon sink of a world with no land (GEN.814, `planetgen.physics.waterworld_carbon`): ocean carbonate chemistry under a high-pressure-ice
seal, and seafloor weathering. [R2381, R2489] are the test targets, never inputs of the code.
"""

import math
import random

import pytest

from planetgen.physics import waterworld_carbon as wc
from planetgen.util import draw

EARTH_OCEAN_KG_M2 = 2.7e6
SEALED_OCEAN_KG_M2 = 2.0e8
"""float: An ocean 200 km deep: 2 GPa on the seafloor at 10 m/s2, 1.7 percent of an Earth mass of water."""
G = 10.0
CO2_KG_MOL = 0.04401


def _carbon(ratio, column=SEALED_OCEAN_KG_M2):
    """The carbon of a closed system whose dioxide is `ratio` of the ocean's mass (Kite and Ford's `cC / fW`), mol per square metre."""
    return ratio * column / CO2_KG_MOL


def _pco2(ratio, km=10.0, temperature=300.0, column=SEALED_OCEAN_KG_M2):
    sodium, calcium = wc.leached_cations_mol_m2(km)
    return wc.ocean_pco2(_carbon(ratio, column), temperature, column, G, sodium_mol_m2=sodium, calcium_mol_m2=calcium)


def test_the_carbonate_constants_reproduce_the_values_at_298_k():
    k0, k1, k2, ksp, kw = wc._log10_constants(298.15)
    assert 10.0 ** k0 == pytest.approx(0.0340, rel=0.01)
    assert -k1 == pytest.approx(6.35, abs=0.01) and -k2 == pytest.approx(10.33, abs=0.01)
    assert -ksp == pytest.approx(8.48, abs=0.01) and -kw == pytest.approx(14.0, abs=0.01)


def test_the_carbonate_constants_hold_at_the_ends_of_their_range_and_beyond_use_the_end_values():
    cold, hot = wc.TEMPERATURE_RANGE_K
    assert wc._log10_constants(cold - 50.0) == wc._log10_constants(cold)
    assert wc._log10_constants(hot + 300.0) == wc._log10_constants(hot)
    assert 10.0 ** wc._log10_constants(cold)[0] > 10.0 ** wc._log10_constants(hot)[0]       # carbon dioxide leaves a warmer ocean
    assert wc._log10_constants(hot)[3] < wc._log10_constants(cold)[3]                        # calcite is retrograde: less soluble when hot


def test_the_leached_cations_scale_with_the_crust_and_match_the_concentration_of_kite_and_ford():
    sodium, calcium = wc.leached_cations_mol_m2(1.0)
    assert sodium == pytest.approx(1000.0 * 3000.0 * 0.028 / 0.02299, rel=0.01)
    assert calcium == pytest.approx(1000.0 * 3000.0 * 0.085 / 0.040078, rel=0.01)
    assert wc.leached_cations_mol_m2(7.0) == pytest.approx((7.0 * sodium, 7.0 * calcium))
    assert wc.leached_cations_mol_m2(0.0) == (0.0, 0.0)
    # their cations total 0.02 to 9 Eq per kg of water over 0.1 to 50 km of rock: here the same span for a 200 km ocean
    for km, low, high in ((0.1, 0.005, 0.1), (50.0, 1.0, 9.0)):
        na, ca = wc.leached_cations_mol_m2(km)
        assert low <= (na + 2.0 * ca) / SEALED_OCEAN_KG_M2 <= high, km


def test_the_seafloor_pressure_seals_an_ocean_deeper_than_the_seal():
    assert not wc.is_sealed(wc.seafloor_pressure_pa(EARTH_OCEAN_KG_M2, 9.8), wc.SEAL_PRESSURE_PA[0])
    assert wc.seafloor_pressure_pa(SEALED_OCEAN_KG_M2, G) == pytest.approx(2.0e9)
    assert wc.is_sealed(2.0e9, 1.5e9) and not wc.is_sealed(2.0e9, 2.5e9)
    assert wc.is_sealed(1.0e9, 1.0e9)
    assert wc.seafloor_pressure_pa(-1.0, G) == 0.0


def test_the_dispatch_follows_the_seal_then_the_land():
    assert wc.sink_kind(True, 0.3) == wc.OCEAN_CHEMISTRY and wc.sink_kind(True, 0.0) == wc.OCEAN_CHEMISTRY
    assert wc.sink_kind(False, 0.3) == wc.CONTINENTAL_WEATHERING
    assert wc.sink_kind(False, 0.0) == wc.SEAFLOOR_WEATHERING


@pytest.mark.parametrize("ratio", [1e-5, 1e-3, 1e-2, 2e-2, 3e-2, 1e-1])
@pytest.mark.parametrize("km", [0.1, 3.0, 50.0])
@pytest.mark.parametrize("temperature", [273.15, 300.0, 400.0])
def test_the_carbon_is_conserved(ratio, km, temperature):
    result = _pco2(ratio, km, temperature)
    total = result.atmosphere_mol_m2 + result.dissolved_mol_m2 + result.precipitated_mol_m2
    assert total == pytest.approx(_carbon(ratio), rel=1e-9)
    assert min(result.atmosphere_mol_m2, result.dissolved_mol_m2, result.precipitated_mol_m2) >= 0.0


def test_an_ocean_with_the_alkalinity_of_the_earths_holds_an_earthlike_air_and_ph():
    """Fresh-water constants, no calcite, 2.3 mmol/kg of alkalinity and 2.4 of carbon: within a factor of three of 400 ppm and a pH of 8.1."""
    column = EARTH_OCEAN_KG_M2
    carbon = 0.0024 * column + 4.0e-4 * 1.0e5 / (9.8 * CO2_KG_MOL)
    result = wc.ocean_pco2(carbon, 288.15, column, 9.8, sodium_mol_m2=0.0023 * column, calcium_mol_m2=0.0)
    assert 4.0e-4 / 3.0 < result.pco2_bar < 4.0e-4 * 3.0
    assert result.ph == pytest.approx(8.1, abs=0.6)
    assert result.precipitated_mol_m2 == 0.0 and not result.saturated


def test_the_air_of_a_closed_system_rises_with_the_temperature_and_the_carbon_and_falls_with_the_alkalinity():
    base = _pco2(0.03).pco2_bar
    assert _pco2(0.03, temperature=340.0).pco2_bar > base > _pco2(0.03, temperature=280.0).pco2_bar
    assert _pco2(0.05).pco2_bar > base > _pco2(0.025).pco2_bar
    assert _pco2(0.03, km=30.0).pco2_bar < base < _pco2(0.03, km=3.0).pco2_bar


def test_a_world_with_more_cations_than_carbon_has_almost_no_carbon_dioxide_in_the_air():
    result = _pco2(1.0e-3, 10.0)
    assert result.pco2_bar == wc.MIN_PCO2_BAR and result.ph > 12.0
    assert result.precipitated_mol_m2 > 0.99 * _carbon(1.0e-3)
    assert wc.ocean_pco2(0.0, 300.0, SEALED_OCEAN_KG_M2, G, sodium_mol_m2=1.0e6, calcium_mol_m2=1.0e6).pco2_bar == 0.0


def test_calcite_holds_the_dissolved_carbon_of_a_calcium_rich_ocean_below_a_mole_per_kilogram():
    """Kite and Ford: calcite buffers the dissolved carbon to about 0.1 mol/kg or less."""
    for ratio in (2e-2, 3e-2, 5e-2):
        result = _pco2(ratio)
        assert result.saturated
        assert result.dissolved_mol_m2 / SEALED_OCEAN_KG_M2 < 1.0


def test_the_ocean_takes_almost_all_the_carbon_a_bare_atmosphere_would_hold():
    """Without the ocean the carbon of a ratio of 0.03 would be a pressure of 600 bar; with it, a few bar."""
    result = _pco2(0.03)
    in_air_alone = _carbon(0.03) * CO2_KG_MOL * G / 1.0e5
    assert in_air_alone == pytest.approx(600.0)
    assert result.pco2_bar < in_air_alone / 50.0
    assert result.atmosphere_mol_m2 < 0.02 * _carbon(0.03)


def test_the_prior_of_kite_and_ford_puts_a_share_of_the_worlds_at_a_fraction_of_a_bar_to_twenty_bar():
    """Their carbon (cC / fW 1e-6 to 1e-1) and rock layer (0.1 to 50 km) at 300 K: pCO2 'frequently 0.2 to 20 bar'."""
    rng = random.Random(814)
    pressures = []
    for _ in range(300):
        km = 10.0 ** rng.uniform(-1.0, math.log10(50.0))
        ratio = 10.0 ** rng.uniform(-6.0, -1.0)
        pressures.append(_pco2(ratio, km).pco2_bar)
    sweet = sum(0.2 <= p <= 20.0 for p in pressures) / len(pressures)
    assert 0.1 < sweet < 0.5
    assert sum(p < 1.0e-6 for p in pressures) / len(pressures) > 0.4                # the cations hold the carbon whole
    assert max(pressures) < 1000.0


def test_the_ocean_must_be_there_and_the_air_must_have_mass():
    with pytest.raises(ValueError):
        wc.ocean_pco2(1.0e6, 300.0, 0.0, G, sodium_mol_m2=0.0, calcium_mol_m2=0.0)
    with pytest.raises(ValueError):
        wc.ocean_pco2(1.0e6, 300.0, 1.0e6, 0.0, sodium_mol_m2=0.0, calcium_mol_m2=0.0)


# --- seafloor weathering ----------------------------------------------------------------------------------------------------

LAW = wc.SeafloorLaw(activation_j_mol=75.0e3, ph_exponent=0.27, spreading_exponent=0.47)
"""The posterior medians of [R2489: Table 1]."""


def _surface(p_bar):
    return 255.0 + 25.0 * math.log10(1.0 + p_bar / 3.0e-4) / 3.0


def test_the_seafloor_rate_is_the_modern_one_at_the_modern_state():
    assert wc.seafloor_rate_over_modern(wc.CO2_MODERN_BAR, wc.TS_MODERN_K, wc.HEAT_FLUX_EARTH_W_M2, LAW) == pytest.approx(1.0)
    assert wc.pore_temperature_k(285.0) == pytest.approx(274.037 + 9.0)


def test_the_seafloor_rate_rises_with_the_temperature_the_carbon_dioxide_and_the_heat_flow():
    base = wc.seafloor_rate_over_modern(1.0e-3, 289.0, 0.0614, LAW)
    assert wc.seafloor_rate_over_modern(1.0e-3, 300.0, 0.0614, LAW) > base
    assert wc.seafloor_rate_over_modern(1.0e-2, 289.0, 0.0614, LAW) > base
    assert wc.seafloor_rate_over_modern(1.0e-3, 289.0, 0.12, LAW) > base
    assert wc.seafloor_rate_over_modern(1.0e-2, 289.0, 0.0614, LAW) / base == pytest.approx(10.0 ** (0.27 / 1.34), rel=1e-6)
    assert wc.seafloor_rate_over_modern(1.0e-3, 289.0, 0.0, LAW) == 0.0 or wc.seafloor_rate_over_modern(1.0e-3, 289.0, 0.0, LAW) < base


def test_the_modern_seafloor_sink_is_the_small_share_of_the_outgassing_of_the_source():
    assert wc.SEAFLOOR_SHARE_OF_OUTGASSING == pytest.approx(0.45 / 6.6)
    assert 0.06 < wc.SEAFLOOR_SHARE_OF_OUTGASSING < 0.07


def test_the_seafloor_balance_is_set_by_the_temperature_more_than_by_the_pressure_and_takes_the_supply_when_that_is_less():
    """The pCO2 enters to the power 0.2 and the temperature through an Arrhenius law: a 300 K world balances at about half a bar, a 288 K one needs
    hundreds of bar, a 310 K one a few millibar; and a supply under the balance stays in the air whole."""
    at = lambda t: wc.seafloor_pco2_bar(1.0e5, 0.0614, lambda p: t, LAW, 1.5)
    assert 0.1 < at(300.0) < 2.0
    assert at(288.0) > 100.0 and 1.0e-3 < at(310.0) < 1.0e-2
    assert at(270.0) > at(288.0) > at(300.0) > at(310.0)
    pco2 = at(300.0)
    assert wc.seafloor_pco2_bar(pco2 / 10.0, 0.0614, lambda p: 300.0, LAW, 1.5) == pytest.approx(pco2 / 10.0)
    assert wc.seafloor_pco2_bar(0.0, 0.0614, lambda p: 300.0, LAW, 1.5) == 0.0


def test_the_seafloor_balance_sits_where_the_two_rates_are_equal():
    pco2 = wc.seafloor_pco2_bar(100.0, 0.0614, _surface, LAW, 1.5)
    sink = wc.SEAFLOOR_SHARE_OF_OUTGASSING * wc.seafloor_rate_over_modern(pco2, _surface(pco2), 0.0614, LAW)
    assert sink == pytest.approx(1.0, rel=1e-6)


def test_a_higher_activation_energy_lowers_the_balance_of_a_warm_world_and_a_smaller_outgassing_lowers_it_too():
    warm = lambda p: 320.0
    cold_law = wc.SeafloorLaw(40.0e3, 0.27, 0.47)
    hot_law = wc.SeafloorLaw(110.0e3, 0.27, 0.47)
    assert wc.seafloor_pco2_bar(1.0e4, 0.0614, warm, hot_law, 1.5) < wc.seafloor_pco2_bar(1.0e4, 0.0614, warm, cold_law, 1.5)
    assert wc.seafloor_pco2_bar(1.0e4, 0.03, warm, LAW, 1.5) < wc.seafloor_pco2_bar(1.0e4, 0.12, warm, LAW, 1.5)


def test_with_no_heat_flow_there_is_no_outgassing_and_the_air_has_none():
    assert wc.seafloor_pco2_bar(10.0, 0.0, _surface, LAW, 1.5) == wc.MIN_PCO2_BAR


# --- the draws --------------------------------------------------------------------------------------------------------------

def test_the_drawn_brackets_stay_inside_their_ranges_and_are_recorded():
    recorded = {}
    with draw.bound(814):
        worlds = [wc.draw_waterworld(recorded) for _ in range(30)]
        law = wc.draw_seafloor_law(recorded)
    assert all(wc.SEAL_PRESSURE_PA[0] <= w.seal_pa <= wc.SEAL_PRESSURE_PA[1] for w in worlds)
    assert all(wc.ALTERED_CRUST_KM[0] <= w.altered_km <= wc.ALTERED_CRUST_KM[1] for w in worlds)
    assert len({w.seal_pa for w in worlds}) > 1
    assert wc.SEAFLOOR_ACTIVATION_KJ_MOL[0] * 1e3 <= law.activation_j_mol <= wc.SEAFLOOR_ACTIVATION_KJ_MOL[1] * 1e3
    assert 0.0 <= law.ph_exponent <= 0.5 and 0.0 <= law.spreading_exponent <= 1.0
    assert set(recorded) == {"waterworld_seal_pa", "waterworld_altered_km", "seafloor_activation_kj_mol", "seafloor_ph_exponent",
                             "seafloor_spreading_exponent"}


def test_the_seafloor_rate_is_zero_on_a_body_too_cold_for_pore_water():
    law = wc.SeafloorLaw(activation_j_mol=55000.0, ph_exponent=0.46, spreading_exponent=0.22)
    assert wc.seafloor_rate_over_modern(5.0, 5.0, 0.001, law) == 0.0
