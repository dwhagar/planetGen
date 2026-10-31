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
"""float: Earth's absorbed sunlight with the clouds left in (the model's cloud is the share `CLOUD_LONGWAVE_SHARE` of its outgoing flux, IR12)."""
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


def test_with_the_vapour_free_the_doubling_warms_more_than_with_it_held(monkeypatch):
    """Clear sky: the cloud falls with the surface humidity (a fixed amount of vapour is a drier air as it warms), which takes back what the vapour gives."""
    monkeypatch.setattr(rc, "CLOUD_LONGWAVE_SHARE", 0.0)
    base = rc.solve(_co2(280.0), 9.81, EARTH_FLUX)
    held = rc.solve(_co2(560.0), 9.81, EARTH_FLUX, water_fraction=base.water_fraction)
    free = rc.solve(_co2(560.0), 9.81, EARTH_FLUX)
    assert free.surface_temperature_k - base.surface_temperature_k > held.surface_temperature_k - base.surface_temperature_k + 0.3


def test_earth_sends_out_what_it_absorbs_and_sits_where_the_cloudy_earth_does(earth):
    assert earth.converged and earth.outgoing_w_m2 == pytest.approx(EARTH_FLUX, rel=0.01)
    assert 280.0 < earth.surface_temperature_k < 290.0                 # 288 K measured; with its vapour fixed at 1 kPa the surface is 284 K
    # The measured cold-point and lapse-rate tropopauses are 190 to 217 K (tropics to mid-latitudes), but they sit under an ozone-heated stratosphere; the model's stratosphere is
    # clear-sky radiative equilibrium without ozone (a limit it states: 140 to 160 K here), so the adiabat meets it colder, and the level is one of 40 (a step of about 12 K, GEN.852).
    assert 155.0 < earth.tropopause_temperature_k < 230.0
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


def test_venus_with_the_sourced_far_wing_alone_is_warmer_than_without_and_cooler_than_measured():
    """GEN.1019: the Perrin and Hartmann wing takes the clear-sky 92 bar air from 489 K to about 540 K; the rest of Venus's greenhouse is the cloud."""
    venus = rc.solve(VENUS, 8.87, 160.0)
    assert venus.converged and 500.0 < venus.surface_temperature_k < 600.0


def test_venus_with_the_acid_cloud_is_near_its_measured_surface_and_the_droplet_size_moves_it():
    """GEN.1020: a cloud 25 deep at 550 nm of 1 micron droplets (the middle of the published modes) brings the clear-sky air to about 757 K
    (735 K measured); droplets of 0.3 micron, the smallest published mode, leave it near 620 K."""
    cloudy = rc.solve(VENUS, 8.87, 160.0, cloud=(25.0, 1.0))
    small = rc.solve(VENUS, 8.87, 160.0, cloud=(25.0, 0.3))
    assert cloudy.converged and small.converged
    assert 700.0 < cloudy.surface_temperature_k < 800.0 and small.surface_temperature_k < cloudy.surface_temperature_k - 50.0


def test_titan_with_its_haze_takes_its_measured_surface_and_without_it_is_too_warm():
    """GEN.718: Titan came out twice too warm with a clear sky. Its haze absorbs about 60 percent of the sunlight high up (McKay et al. 1991)."""
    clear = rc.solve(TITAN, 1.35, 2.61)
    hazy = rc.solve(TITAN, 1.35, 2.61, aloft_share=0.6, aloft_pressure_pa=3.0e3)
    assert clear.converged and hazy.converged
    assert 85.0 < hazy.surface_temperature_k < 102.0                     # 94 K measured
    assert clear.surface_temperature_k - hazy.surface_temperature_k > 8.0


@pytest.mark.parametrize("name, partials, gravity, effective_k, measured_k, tolerance", [
    ("Jupiter", {"H2": 0.864e5, "He": 0.136e5, "CH4": 0.002e5}, 24.79, 125.0, 166.0, 0.03),
    ("Uranus", {"H2": 0.825e5, "He": 0.152e5, "CH4": 0.023e5}, 8.69, 59.1, 76.0, 0.06),
    ("Neptune", {"H2": 0.80e5, "He": 0.19e5, "CH4": 0.015e5}, 11.15, 59.3, 72.0, 0.10),
])
def test_the_giants_sit_near_their_measured_temperature_at_one_bar(name, partials, gravity, effective_k, measured_k, tolerance):
    """GEN.1030: the methane pseudo-adiabat and the ortho-para equilibrium of cold hydrogen bring the ice giants' 1 bar temperature down from 20 percent
    too warm (the effective temperatures are the measured ones)."""
    giant = rc.warming(partials, gravity, SIGMA * effective_k ** 4)
    assert giant.converged and giant.surface_temperature_k == pytest.approx(measured_k, rel=tolerance), name


def test_the_warming_factor_is_the_model_in_the_repos_variable(earth):
    warm = rc.warming(_co2(280.0), 9.81, EARTH_FLUX)
    assert warm.surface_temperature_k == pytest.approx(earth.surface_temperature_k)
    assert warm.factor == pytest.approx((warm.surface_temperature_k / warm.effective_temperature_k) ** 4 - 1.0)
    assert warm.effective_temperature_k == pytest.approx((EARTH_FLUX / SIGMA) ** 0.25)


# ---- the tropopause is a model level (GEN.852 part 3, row IR10) ----

def test_the_surface_temperature_does_not_jump_when_the_tropopause_moves_to_another_level():
    air = {"N2": 0.78e5, "O2": 0.21e5, "CO2": 40.0, "H2O": 1500.0}
    fluxes = (244.0, 245.0, 246.0, 247.0, 248.0)
    solutions = [rc.solve(air, 9.81, f) for f in fluxes]
    pressures = [s.tropopause_pressure_pa for s in solutions]
    assert len(set(pressures)) > 1                    # the sweep crosses a change of level (the sweep of 244 to 248 W/m2 moves it from 10.1 to 12.7 kPa; GEN.1032 moved the sweep when the far wing and the acid cloud changed the model)
    steps = np.diff([s.surface_temperature_k for s in solutions])
    assert (steps > 0).all()
    assert steps.max() - steps.min() < 0.15           # about 0.19 K per W/m2 throughout; the move of a level costs under a tenth of a kelvin
