"""
The carbon balance of a world with land and an ocean (GEN.1713, `planetgen.physics.continental_weathering`): the weathering law of Graham and Pierrehumbert 2020 with its limits,
normalised to the Earth, and the seafloor sink beside it.
"""

import pytest

from planetgen.physics import continental_weathering as cw, surface_pressure as sp, waterworld_carbon as wc

EARTH_LAW = wc.SeafloorLaw(activation_j_mol=75.0e3, ph_exponent=0.27, spreading_exponent=0.47)
LAND = sp.LAND_FRACTION_EARTH
Q = sp.HEAT_FLUX_EARTH_W_M2


def _balance(supply_bar=10.0, flux=Q, land=LAND, temperature=sp.TS_MODERN_K, absorbed=cw.EARTH_ABSORBED_FLUX_W_M2, cycle=(1.0, 1.5), law=EARTH_LAW):
    return cw.pco2_bar(supply_bar, flux, land, temperature, absorbed, cycle, law)


def test_the_cycle_draws_stay_in_their_brackets_and_are_recorded():
    recorded = {}
    lambda_scale, exponent = cw.draw_cycle(recorded)
    assert cw.LAMBDA_SCALE[0] <= lambda_scale <= cw.LAMBDA_SCALE[1]
    assert sp.OUTGASSING_EXPONENT[0] <= exponent <= sp.OUTGASSING_EXPONENT[1]
    assert recorded == {"weathering_lambda_scale": lambda_scale, "outgassing_exponent": exponent}


def test_the_law_at_the_modern_reference_is_the_paper_s_form_and_its_level_is_the_factor_of_0_36_the_paper_misses_by():
    per_land = cw.rate_mol_m2_yr(289.0, sp.CO2_MODERN_BAR, cw.EARTH_ABSORBED_FLUX_W_M2)
    assert per_land == pytest.approx(0.0184, rel=0.01)                                 # Eq. 40 with Table 1, per square metre of land
    assert per_land * 0.29 / 0.0147 == pytest.approx(0.36, abs=0.02)                  # GP20's own reference flux (V_ref over the Earth's area) is 2.8 times that


@pytest.mark.parametrize("lambda_scale, m, law", [(0.1, 1.0, None), (1.0, 1.5, EARTH_LAW), (10.0, 2.0, EARTH_LAW)])
def test_the_earth_is_the_balance_at_its_modern_flux_land_and_temperature(lambda_scale, m, law):
    assert cw.rate_over_modern(sp.TS_MODERN_K, sp.CO2_MODERN_BAR, cw.EARTH_ABSORBED_FLUX_W_M2, lambda_scale) == pytest.approx(1.0)
    assert _balance(cycle=(lambda_scale, m), law=law) == pytest.approx(sp.CO2_MODERN_BAR, rel=1e-6)


def test_the_runoff_is_held_under_the_energy_limit_and_the_rate_under_the_supply_limit():
    free = cw.runoff_m_yr(300.0, 1.0e6)                                                # light to spare: the precipitation law alone, Gamma times 0.99 (1 + 0.03 x 12)
    assert free == pytest.approx(0.2 * 0.99 * 1.36)
    starved = cw.runoff_m_yr(300.0, 50.0)                                              # 50 W per square metre cannot evaporate that much water
    assert starved < free / 1.5 and starved == pytest.approx(0.2 * 50.0 / (1.918e9 * (300.0 / (300.0 - 33.91)) ** 2) * 3.15576e7)
    assert cw.runoff_m_yr(240.0, 238.0) == 0.0 and cw.runoff_m_yr(20.0, 238.0) == 0.0     # the linear rise has run out; the fit is undefined under its offset
    ceiling = cw.ALPHA / (cw.MOLAR_MASS_KG_MOL * cw.SPECIFIC_AREA_M2_KG * cw.SOIL_AGE_YR)
    assert cw.rate_mol_m2_yr(340.0, 1.0, 1.0e6, 10.0) < ceiling                        # the supply limit is a ceiling that no heat or CO2 passes
    assert cw.rate_mol_m2_yr(340.0, 1.0, 1.0e6, 10.0) > 0.2 * ceiling


def test_more_volcanism_or_less_land_raises_the_balance_and_a_warmer_air_lowers_it():
    base = _balance()
    assert _balance(flux=2.0 * Q) > base
    assert _balance(land=0.5 * LAND) > base
    assert _balance(temperature=sp.TS_MODERN_K + 10.0) < base


def test_the_seafloor_sink_lowers_the_carbon_dioxide_of_a_warm_world_and_the_energy_limit_lifts_it():
    hot = 315.0
    assert _balance(temperature=hot, law=EARTH_LAW) < _balance(temperature=hot, law=None)                  # two sinks take more than one
    assert _balance(temperature=hot, law=None, absorbed=60.0) > _balance(temperature=hot, law=None)        # a dimmer sun evaporates less water, so less runoff and less weathering


def test_a_frozen_world_with_no_rain_keeps_its_supply_and_the_balance_never_exceeds_it():
    assert _balance(supply_bar=0.02, temperature=200.0) == 0.02                         # no runoff, no weathering: the CO2 builds up until the ice melts
    assert _balance(supply_bar=1e-5) == pytest.approx(1e-5)
    assert _balance(supply_bar=0.0) == 0.0
    assert _balance(flux=0.0) == sp.MIN_PCO2_BAR


def test_no_ocean_floor_means_no_seafloor_sink_and_the_land_alone_sets_the_balance():
    with_sea = _balance(temperature=310.0)
    alone = _balance(temperature=310.0, law=None)
    assert alone > with_sea
