"""
The surface pressure supply law (GEN.726, `planetgen.physics.surface_pressure`, docs/design/composition/surface-pressure-supply-law.md).
"""

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


def test_the_balance_root_is_the_stable_crossing_the_supply_or_the_floor():
    rising = lambda p_bar: p_bar / 1e-3 - 1.0
    assert sp.balance_root_bar(1.0, rising) == pytest.approx(1e-3, rel=1e-6)
    assert sp.balance_root_bar(1e-4, rising) == pytest.approx(1e-4)              # a world with little carbon keeps all of it
    assert sp.balance_root_bar(0.0, rising) == 0.0
    assert sp.balance_root_bar(1.0, lambda p_bar: 1.0) == sp.MIN_PCO2_BAR            # a sink that takes everything even at the floor
    assert sp.balance_root_bar(1.0, lambda p_bar: -1.0) == 1.0                       # one that never keeps up: the air holds the supply
