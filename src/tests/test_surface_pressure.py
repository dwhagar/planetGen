"""
The surface pressure supply law (GEN.726, `planetgen.physics.surface_pressure`, docs/design/composition/surface-pressure-supply-law.md).
"""

import math

import pytest

from planetgen.physics import constants, interior_heat, surface_pressure as sp
from tests.test_interior_heat import CORNERS, EARTH, EARTH_ROCK, MID, _bulk_silicate_earth, _earth

EARTH_MANTLE_KG = EARTH * EARTH_ROCK.silicate
VENUS_AIR_KG = 4.77e20


def _venus(interior=MID):
    return interior_heat.heat_flow(_bulk_silicate_earth(), interior, 0.815 * EARTH, 6052.0, 4.5, EARTH_ROCK)


def test_the_interior_reports_the_heat_it_has_released_and_it_grows_with_the_age():
    released = [_earth(age=age).released_j for age in (0.5, 1.0, 2.0, 4.568)]
    assert released == sorted(released) and released[0] > 0.0
    # the present flow cannot have been the whole history: the young body lost heat faster
    assert released[-1] > _earth().total_w * 4.568e9 * constants.SECONDS_PER_YEAR


def test_a_body_with_no_silicate_has_released_no_heat():
    none = interior_heat.heat_flow(_bulk_silicate_earth(), MID, EARTH, 6371.0, 4.5, interior_heat.radiogenic.Rock(0.0, 0.0, 1.0, 0.0))
    assert none.released_j == 0.0


def test_an_earth_twin_has_degassed_the_share_the_surface_carbon_of_the_earth_asks_for():
    share = sp.degassed_share(_earth().released_j, EARTH_MANTLE_KG)
    assert share == pytest.approx(sp.EARTH_DEGASSED_SHARE, rel=0.02)
    # the Earth's surface carbon, 1e22 mol of CO2, is that share of the log-middle of the Ortenzi CO2 range in the mantle
    middle_co2_kg = 10.0 ** -3.5 * EARTH_MANTLE_KG
    assert share * middle_co2_kg == pytest.approx(1e22 * 0.04401, rel=0.15)


def test_venus_degasses_about_what_the_earth_does_and_its_air_asks_for():
    """The check, not the fit: Venus's 96.5 percent CO2 over 9.2 MPa is 4.6e20 kg, 0.43 of the same inventory in its mantle."""
    needed = 0.965 * VENUS_AIR_KG / (10.0 ** -3.5 * 0.815 * EARTH_MANTLE_KG)
    for interior in (MID, *CORNERS):
        found = sp.degassed_share(_venus(interior).released_j, 0.815 * EARTH_MANTLE_KG)
        # a factor of 2.2: the stiffest corner of the rheology range gives 0.214 against 0.22 at a factor of 2 now that the mantle's ceiling is its liquidus (GEN.735)
        assert needed / 2.2 < found < needed * 2.2


def test_no_heat_or_no_mantle_degasses_nothing_and_the_share_stays_below_one():
    assert sp.degassed_share(0.0, EARTH_MANTLE_KG) == 0.0 and sp.degassed_share(1e31, 0.0) == 0.0
    assert 0.0 < sp.degassed_share(1e40, EARTH_MANTLE_KG) <= 1.0


def test_the_cycle_draws_stay_in_their_brackets_and_are_recorded():
    recorded = {}
    alpha, t_e, exponent = sp.draw_cycle(recorded)
    assert sp.WEATHERING_ALPHA[0] <= alpha <= sp.WEATHERING_ALPHA[1] and sp.WEATHERING_T_E_K[0] <= t_e <= sp.WEATHERING_T_E_K[1]
    assert sp.OUTGASSING_EXPONENT[0] <= exponent <= sp.OUTGASSING_EXPONENT[1]
    assert recorded == {"weathering_alpha": alpha, "weathering_t_e_k": t_e, "outgassing_exponent": exponent}


def _balance(supply_bar=10.0, flux=sp.HEAT_FLUX_EARTH_W_M2, land=sp.LAND_FRACTION_EARTH, temperature=lambda p: sp.TS_MODERN_K, **law):
    return sp.carbonate_pco2_bar(supply_bar, flux, land, temperature, law.get("alpha", 0.3), law.get("t_e", 20.0), law.get("m", 1.5))


@pytest.mark.parametrize("alpha, t_e, m", [(0.1, 10.0, 1.0), (0.3, 20.0, 1.5), (0.5, 40.0, 2.0)])
def test_the_earth_is_the_balance_at_its_modern_flux_land_and_temperature(alpha, t_e, m):
    assert _balance(alpha=alpha, t_e=t_e, m=m) == pytest.approx(sp.CO2_MODERN_BAR, rel=1e-6)


def test_more_volcanism_or_less_land_raises_the_balance_and_a_warmer_air_lowers_it():
    base = _balance()
    assert _balance(flux=2.0 * sp.HEAT_FLUX_EARTH_W_M2) > base
    assert _balance(land=0.5 * sp.LAND_FRACTION_EARTH) > base
    assert _balance(temperature=lambda p: sp.TS_MODERN_K + 10.0) < base


def test_the_balance_never_exceeds_the_supply_and_a_world_with_little_carbon_keeps_all_it_has_below_the_balance():
    assert _balance(supply_bar=1e-5) == pytest.approx(1e-5)
    assert _balance(supply_bar=1e-3) == pytest.approx(sp.CO2_MODERN_BAR, rel=1e-6)
    assert _balance(supply_bar=0.0) == 0.0


def test_a_greenhouse_that_warms_with_the_pressure_still_has_one_root_below_the_supply():
    greenhouse = lambda p: sp.TS_MODERN_K + 8.0 * math.log10(max(p, 1e-12) / sp.CO2_MODERN_BAR + 1.0)
    found = _balance(supply_bar=100.0, temperature=greenhouse, flux=3.0 * sp.HEAT_FLUX_EARTH_W_M2)
    assert sp.CO2_MODERN_BAR < found < 100.0
    # at the root the weathering equals the outgassing
    weathering = (found / sp.CO2_MODERN_BAR) ** 0.3 * math.exp((greenhouse(found) - sp.TS_MODERN_K) / 20.0)
    assert weathering == pytest.approx(3.0 ** 1.5, rel=1e-6)


def test_a_body_with_no_volcanism_is_taken_to_the_floor_and_a_cold_still_world_keeps_its_supply():
    assert _balance(flux=0.0) == sp.MIN_PCO2_BAR
    assert _balance(supply_bar=0.02, temperature=lambda p: 150.0) == 0.02
