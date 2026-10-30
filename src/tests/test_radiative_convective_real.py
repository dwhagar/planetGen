"""
The radiative-convective model on the real downloaded tables (GEN.715) against numbers measured in the sky and the checks GEN.717 and GEN.718
asked for: the response of Earth to a doubling of carbon dioxide with the vapour held, Titan's surface with and without its haze, and the
surface temperatures of Earth, Mars and Venus. Earth, Venus, Mars and Titan are test targets, never inputs of the code. Skipped on a machine
whose cache folder does not hold the two infrared datasets (`python3 -m planetgen.cli.datasets status`); the offline fixtures of
`test_line_opacity.py` and `test_radiative_convective.py` cover the code.
"""

import numpy as np
import pytest

from planetgen.datasets import fetch, hitran_cia, ir_ktables
from planetgen.physics import constants
from planetgen.physics import radiative_convective as rc

pytestmark = pytest.mark.skipif(
    any(fetch.locate(dataset).missing for dataset in (ir_ktables.DATASET, hitran_cia.DATASET)),
    reason="the infrared datasets are not in this machine's cache folder")

SIGMA = constants.STEFAN_BOLTZMANN_CONSTANT
EARTH = {"N2": 78084.0, "O2": 20946.0, "Ar": 934.0, "CO2": 41.0, "H2O": 1000.0}
EARTH_FLUX = 240.0
"""float: Earth's absorbed sunlight with the clouds left in (the model has none, so its surface is a clear-sky one, cooler than the measured 288 K)."""
MARS = {"CO2": 600.0 * 0.953, "N2": 600.0 * 0.027, "Ar": 600.0 * 0.016}
VENUS = {"CO2": 0.965 * 92.0e5, "N2": 0.035 * 92.0e5}
TITAN = {"N2": 1.47e5 * 0.985, "CH4": 1.47e5 * 0.015}


def _co2(ppm):
    return {**EARTH, "CO2": ppm * 1.0e-6 * 1.0e5}


@pytest.fixture(scope="module")
def earth():
    return rc.solve(_co2(280.0), 9.81, EARTH_FLUX)


def test_a_doubling_of_carbon_dioxide_with_the_vapour_held_warms_by_about_a_kelvin(earth):
    """GEN.717: the model was off against this (0.3 to 0.6 K). Held vapour is the no-feedback response (about 1.2 K from the Planck response alone
    with the stratosphere adjusted, Myhre et al. 1998 forcing 3.7 W/m2 over a Planck feedback of about 3.2 W/m2/K)."""
    held = earth.water_fraction
    base = rc.solve(_co2(280.0), 9.81, EARTH_FLUX, water_fraction=held)
    doubled = rc.solve(_co2(560.0), 9.81, EARTH_FLUX, water_fraction=held)
    assert base.converged and doubled.converged
    assert 0.8 < doubled.surface_temperature_k - base.surface_temperature_k < 1.4


def test_with_the_vapour_free_the_doubling_warms_more_than_with_it_held(earth):
    held = rc.solve(_co2(560.0), 9.81, EARTH_FLUX, water_fraction=earth.water_fraction)
    free = rc.solve(_co2(560.0), 9.81, EARTH_FLUX)
    assert free.surface_temperature_k - earth.surface_temperature_k > held.surface_temperature_k - earth.surface_temperature_k


def test_earth_sends_out_what_it_absorbs_and_sits_where_a_clear_sky_earth_does(earth):
    assert earth.converged and earth.outgoing_w_m2 == pytest.approx(EARTH_FLUX, rel=0.01)
    assert 270.0 < earth.surface_temperature_k < 282.0                 # 288 K measured, with clouds that warm and cool; clear-sky models give 275 to 280
    assert 170.0 < earth.tropopause_temperature_k < 230.0
    assert 5.0e3 < earth.tropopause_pressure_pa < 4.0e4                # the tropopause: 100 to 200 hPa


def test_without_the_vapour_or_the_carbon_dioxide_the_earth_is_colder():
    dry = rc.solve({"N2": 78084.0, "O2": 20946.0, "Ar": 934.0, "CO2": 41.0}, 9.81, EARTH_FLUX).surface_temperature_k
    bare = rc.solve({"N2": 1.0e5}, 9.81, EARTH_FLUX).surface_temperature_k
    wet = rc.solve(EARTH, 9.81, EARTH_FLUX).surface_temperature_k
    assert bare < dry < wet
    assert wet - dry > 8.0                                              # water vapour is the larger part of the natural greenhouse


def test_mars_sits_near_its_measured_surface_temperature():
    mars = rc.solve(MARS, 3.71, 110.0)
    assert mars.converged and 205.0 < mars.surface_temperature_k < 225.0          # 210 to 215 K measured


def test_venus_is_hot_when_its_far_wings_are_included():
    """The k-tables cut every line at 25 cm-1; the labelled far-wing stand-in (IR3) fills what Koll 2019 says a dense carbon dioxide air holds."""
    venus = rc.solve(VENUS, 8.87, 160.0)
    assert venus.converged and 680.0 < venus.surface_temperature_k < 770.0        # 735 K measured


def test_titan_with_its_haze_takes_its_measured_surface_and_without_it_is_too_warm():
    """GEN.718: Titan came out twice too warm with a clear sky. Its haze absorbs about 60 percent of the sunlight high up (McKay et al. 1991)."""
    clear = rc.solve(TITAN, 1.35, 2.61)
    hazy = rc.solve(TITAN, 1.35, 2.61, aloft_share=0.6, aloft_pressure_pa=3.0e3)
    assert clear.converged and hazy.converged
    assert 85.0 < hazy.surface_temperature_k < 102.0                     # 94 K measured
    assert clear.surface_temperature_k - hazy.surface_temperature_k > 8.0


def test_the_warming_factor_is_the_model_in_the_repos_variable(earth):
    warm = rc.warming(_co2(280.0), 9.81, EARTH_FLUX)
    assert warm.surface_temperature_k == pytest.approx(earth.surface_temperature_k)
    assert warm.factor == pytest.approx((warm.surface_temperature_k / warm.effective_temperature_k) ** 4 - 1.0)
    assert warm.effective_temperature_k == pytest.approx((EARTH_FLUX / SIGMA) ** 0.25)


# ---- the tropopause is a model level (GEN.852 part 3, row IR10) ----

def test_the_surface_temperature_does_not_jump_when_the_tropopause_moves_to_another_level():
    air = {"N2": 0.78e5, "O2": 0.21e5, "CO2": 40.0, "H2O": 1500.0}
    fluxes = (237.0, 238.0, 239.0, 240.0, 241.0)
    solutions = [rc.solve(air, 9.81, f) for f in fluxes]
    pressures = [s.tropopause_pressure_pa for s in solutions]
    assert len(set(pressures)) > 1                    # the sweep crosses a change of level (the sweep of 237 to 241 W/m2 moves it between 12.7 and 20 kPa)
    steps = np.diff([s.surface_temperature_k for s in solutions])
    assert (steps > 0).all()
    assert steps.max() - steps.min() < 0.15           # about 0.4 K per W/m2 throughout; the move of a level costs under a tenth of a kelvin
