"""
Climate sensitivity sweep (GEN.975, TEST.509): the surface-air fixed point of an Earth twin, solved again at a range of absorbed sunlight from a cooling
to a warming sweep. Every solve must converge, the temperature must rise with the light, and the response to a change of the absorbed flux must stay
inside the published range of Earth's climate sensitivity: an equilibrium warming of 2 to 5 K for a doubling of carbon dioxide (the very likely
range of IPCC AR6, 2021, chapter 7), a forcing of 3.7 W/m2 for the doubling, so 0.54 to 1.35 K per W/m2 of absorbed flux.

The twin is the Earth itself: its dry air (78.08 percent nitrogen, 20.95 oxygen, 0.93 argon, 420 ppm carbon dioxide, 1 bar), its gravity and radius, and its
ocean (1.4e21 kg over 71 percent of the surface), held in the fixed point through a stand-in for the mixture draw, so the air is the Earth's and the fixed point
still brings the vapour, the pressure and the temperature. (TEST.509: the twin of GEN.975 was a random drawn body at 313 K with a hydrogen and carbon dioxide air,
whose response said nothing of the Earth's.) The model stands about 6.5 K under the Earth's 288 K at the Earth's own 238 W/m2 absorbed, because it has no cloud
(GEN.974), so the absorbed flux of the twin is set, by bisection, to the flux at which the fixed point stands at the measured 288 K: the sweep then measures
the response of a 288 K Earth, the number the published range is about.

Where the fixed point stands: it converges and warms with the light, with a response of 0.40 K per W/m2 at 288 K, under the published lower bound of 0.54. Of the
feedbacks of IPCC AR6 (2021) table 7.10 only the Planck response and the water vapour with the lapse rate are in it (1 / (3.22 - 1.30) = 0.52 K per W/m2 by the
assessment's own numbers; the fixed point's outgoing flux rises faster with the temperature, 2.5 W/m2 per K, a gap GEN.973 follows). Adding the surface albedo
(+0.35 W/m2 per K, GEN.972) and the cloud (+0.42, GEN.974) feedbacks gives 1 / (2.5 - 0.77) = 0.58, inside the range. The lower bound is held by a strict
expected failure that names them; when GEN.972, GEN.973 and GEN.974 bring the response inside the range, the test passes, the expected failure fails the run, and
the mark comes off.

The twin is a test target, never an input of the code. The solve is replayed on one stream (`draw.bound`), so the sweep changes the light and nothing else;
the albedo is held in the absorbed flux until the feedbacks that move it (GEN.972 ice-albedo, GEN.974 cloud albedo) are in the fixed point.
"""

import math
import types

import pytest
from scipy.optimize import brentq

from planetgen.physics import atmosphere_state, gases, mixture_draw, surface_air as sa
from planetgen.util import draw

SENSITIVITY_K_PER_W_M2 = (2.0 / 3.7, 5.0 / 3.7)
"""tuple: The very likely range of the equilibrium climate sensitivity (2 to 5 K per doubling of CO2, AR6) over the forcing of a doubling (3.7 W/m2)."""
EARTH_DRY_AIR = {"N2": 0.7808, "O2": 0.2095, "Ar": 0.0093, "CO2": 4.2e-4}
"""dict: The mole fractions of the Earth's dry air (a test target)."""
EARTH_GRAVITY_M_S2, EARTH_RADIUS_M, EARTH_MANTLE_KG = 9.80665, 6.371e6, 4.0e24
EARTH_WATER_KG, EARTH_OCEAN_COVER = 1.4e21, 0.71
"""float: Earth's ocean (a test target)."""
EARTH_SURFACE_K = 288.0
"""float: The Earth's mean surface temperature, K, the temperature the absorbed flux of the twin is set to give."""
START_K = 255.0
"""float: The temperature the first step of the fixed point starts from (the airless Earth's)."""
FACTORS = (0.90, 0.95, 1.0, 1.05, 1.10)
"""tuple: The sweep, as a share of the absorbed flux of the twin: a cooling half and a warming half."""


def _earth_twin(flux_w_m2):
    return sa.Setting(gravity_m_s2=EARTH_GRAVITY_M_S2, radius_m=EARTH_RADIUS_M, mantle_kg=EARTH_MANTLE_KG, released_j=0.0, heat_flux_w_m2=0.087,
                      ocean_cover=EARTH_OCEAN_COVER, water_kg=EARTH_WATER_KG, absorbed_flux_w_m2=flux_w_m2, released_at_seal_j=0.0,
                      subject_of=lambda *_: None, age_gy=4.5)


def _solve(flux_w_m2):
    with draw.bound("climate-sweep"):
        return sa.solve(_earth_twin(flux_w_m2), START_K)


def _hold_the_air_of_the_earth(monkeypatch):
    """The mixture draw stands in with the Earth's dry air, in the moles one bar of it holds over the Earth's area."""
    area = 4.0 * math.pi * EARTH_RADIUS_M ** 2
    mean_molar_mass = sum(x * gases.gas(formula).molar_mass for formula, x in EARTH_DRY_AIR.items())
    moles = 1.0e5 * area / EARTH_GRAVITY_M_S2 / mean_molar_mass
    air = types.SimpleNamespace(family="rocky", mole_fractions=dict(EARTH_DRY_AIR), moles=moles / (EARTH_MANTLE_KG / 0.1), pressure_bar=1.0, recipe=None)
    monkeypatch.setattr(mixture_draw, "draw_mixture", lambda subject, **_: air)


@pytest.fixture(scope="module")
def sweep():
    with pytest.MonkeyPatch.context() as monkeypatch:
        _hold_the_air_of_the_earth(monkeypatch)
        flux = brentq(lambda f: _solve(f).temperature_k - EARTH_SURFACE_K, 200.0, 320.0, xtol=1.0e-3)
        return flux, {factor: _solve(flux * factor) for factor in FACTORS}


def test_every_step_of_the_cooling_and_the_warming_sweep_converges(sweep):
    _, airs = sweep
    assert [factor for factor, air in airs.items() if not air.converged] == []


def test_the_surface_warms_with_the_light(sweep):
    _, airs = sweep
    temperatures = [airs[factor].temperature_k for factor in FACTORS]
    assert temperatures == sorted(temperatures) and len(set(temperatures)) == len(temperatures)


def _responses(sweep):
    flux, airs = sweep
    return [(airs[second].temperature_k - airs[first].temperature_k) / (flux * (second - first)) for first, second in zip(FACTORS[:-1], FACTORS[1:])]


def test_the_twin_stands_at_the_earths_temperature_with_the_earths_air(sweep):
    flux, airs = sweep
    assert airs[1.0].temperature_k == pytest.approx(EARTH_SURFACE_K, abs=0.05)
    assert 220.0 < flux < 275.0                                    # about 7 percent over the Earth's 238 W/m2: the cloud the model lacks
    assert airs[1.0].partials_pa["N2"] / airs[1.0].pressure_pa > 0.75 and 0.9e5 < airs[1.0].pressure_pa < 1.1e5


def test_the_response_to_the_flux_is_not_above_the_published_climate_sensitivity(sweep):
    assert all(0.0 < response < SENSITIVITY_K_PER_W_M2[1] for response in _responses(sweep))


@pytest.mark.xfail(strict=True, reason="GEN.972 and GEN.974: no ice-albedo or cloud feedback in the fixed point yet, so the response (0.40 K per W/m2 at 288 K) is the Planck one "
                                       "with the vapour and the lapse rate, under the published lower bound of 0.54; with them (+0.35 and +0.42 W/m2 per K, AR6 table 7.10) "
                                       "it is 0.58. GEN.1194 (the outgoing flux rising too fast with the temperature) moves it too")
def test_the_response_to_the_flux_is_not_below_the_published_climate_sensitivity(sweep):
    assert all(response > SENSITIVITY_K_PER_W_M2[0] for response in _responses(sweep))


def test_the_vapour_of_the_ocean_is_a_feedback_in_the_fixed_point(sweep):
    """GEN.973: the fixed point already lets the vapour follow the temperature (the ocean holds it at a share of saturation at each step's own temperature, and the
    column's humidity follows Manabe and Wetherald). Hold the air of the unchanged twin and ask the table for each flux of the sweep: the response is smaller than the
    fixed point's, so the vapour feedback is in it and the held-air response is the Planck one with the lapse rate."""
    flux, airs = sweep
    held = airs[1.0].partials_pa
    held_response = (atmosphere_state.greenhouse_warming(held, EARTH_GRAVITY_M_S2, flux * FACTORS[-1]).surface_temperature_k
                     - atmosphere_state.greenhouse_warming(held, EARTH_GRAVITY_M_S2, flux * FACTORS[0]).surface_temperature_k) / (flux * (FACTORS[-1] - FACTORS[0]))
    free_response = (airs[FACTORS[-1]].temperature_k - airs[FACTORS[0]].temperature_k) / (flux * (FACTORS[-1] - FACTORS[0]))
    assert 0.0 < held_response < free_response
    assert airs[FACTORS[-1]].partials_pa["H2O"] > 1.5 * airs[FACTORS[0]].partials_pa["H2O"]          # the vapour that carried the extra warming
