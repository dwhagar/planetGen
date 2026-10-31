"""
Climate sensitivity sweep (GEN.975): the surface-air fixed point of an Earth-like twin, solved again at a range of absorbed sunlight from a cooling
to a warming sweep. Every solve must converge, the temperature must rise with the light, and the response to a change of the absorbed flux must stay
inside the published range of Earth's climate sensitivity: an equilibrium warming of 2 to 5 K for a doubling of carbon dioxide (the very likely
range of IPCC AR6, 2021, chapter 7), a forcing of 3.7 W/m2 for the doubling, so 0.54 to 1.35 K per W/m2 of absorbed flux.

Where the fixed point stands: it converges and warms with the light, with a response (about 0.23 K per W/m2 with the greenhouse table of GEN.850; 0.4 with the old grey closure) under the published lower
bound, because the feedbacks that make a climate more sensitive than its Planck response are not in it yet. The lower bound is held by a strict expected
failure that names them; when GEN.972 (ice-albedo), GEN.973 (water vapour) and GEN.974 (cloud albedo) bring the response inside the range, the test passes, the
expected failure fails the run, and the mark comes off.

The body is Earth-like and is a test target, never an input of the code. The solve is replayed on one stream (`draw.bound`), so the sweep changes the
light and nothing else; the albedo is the body's own, held, until the feedbacks that move it (GEN.972 ice-albedo, GEN.974 cloud albedo) are in the fixed point.
"""

import dataclasses

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.physics import atmosphere_state, constants, surface_air as sa
from planetgen.util import draw

SENSITIVITY_K_PER_W_M2 = (2.0 / 3.7, 5.0 / 3.7)
"""tuple: The very likely range of the equilibrium climate sensitivity (2 to 5 K per doubling of CO2, AR6) over the forcing of a doubling (3.7 W/m2)."""
EARTH_WATER_KG, EARTH_OCEAN_COVER = 1.4e21, 0.71
"""float: Earth's ocean (a test target): the twin's own water is a dry planet's, and a sweep through it crosses the freezing of its small sea."""
TEMPERATE_K = (260.0, 340.0)
"""tuple: The surface temperatures, K, between which a twin counts as Earth-like (Earth 288 K) for the sweep."""
FACTORS = (0.90, 0.95, 1.0, 1.05, 1.10)
"""tuple: The sweep, as a share of the absorbed flux of the twin: a cooling half and a warming half."""


def _earth_twin_setting(monkeypatch):
    from planetgen.generation.planet import Planet
    from planetgen.generation.star import Star

    captured = {}
    real = sa.solve

    def spy(setting, start_temperature_k):
        captured["setting"], captured["start"] = setting, start_temperature_k
        return real(setting, start_temperature_k)

    monkeypatch.setattr(sa, "solve", spy)
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    with draw.bound("climate-sweep-star"):
        star = Star(config)
    inner, outer = star.habitable_zone
    for seed in range(3, 60):       # the first draw whose air with Earth's ocean stands at a temperate temperature: the planet a seed makes moves with the generator
        with draw.bound(seed):
            Planet(config, star, star.habitable_zone, 0.5 * (inner + outer), mass=constants.EARTH_MASS_TO_KG, moon_count=0)
        setting = dataclasses.replace(captured["setting"], water_kg=EARTH_WATER_KG, ocean_cover=EARTH_OCEAN_COVER)
        with draw.bound("climate-sweep"):
            air = real(setting, captured["start"])
        if air.converged and TEMPERATE_K[0] < air.temperature_k < TEMPERATE_K[1]:
            monkeypatch.undo()
            return setting, captured["start"]
    monkeypatch.undo()
    pytest.fail("no Earth-mass planet of 57 seeds stands at a temperate temperature with Earth's ocean")


@pytest.fixture(scope="module")
def sweep():
    with pytest.MonkeyPatch.context() as monkeypatch:
        setting, start = _earth_twin_setting(monkeypatch)
    airs = {}
    for factor in FACTORS:
        with draw.bound("climate-sweep"):
            airs[factor] = sa.solve(dataclasses.replace(setting, absorbed_flux_w_m2=setting.absorbed_flux_w_m2 * factor), start)
    return setting.absorbed_flux_w_m2, airs


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


def test_the_response_to_the_flux_is_not_above_the_published_climate_sensitivity(sweep):
    assert all(0.0 < response < SENSITIVITY_K_PER_W_M2[1] for response in _responses(sweep))


@pytest.mark.xfail(strict=True, reason="GEN.972, GEN.973 and GEN.974: no ice-albedo, water-vapour or cloud-albedo feedback in the fixed point yet, so the response "
                                       "is the Planck response with the vapour of the ocean's evaporation, under the published lower bound")
def test_the_response_to_the_flux_is_not_below_the_published_climate_sensitivity(sweep):
    assert all(response > SENSITIVITY_K_PER_W_M2[0] for response in _responses(sweep))


@pytest.fixture(scope="module")
def twin_gravity():
    with pytest.MonkeyPatch.context() as monkeypatch:
        setting, _start = _earth_twin_setting(monkeypatch)
    return setting.gravity_m_s2


def test_the_vapour_of_the_ocean_is_a_feedback_in_the_fixed_point(sweep, twin_gravity):
    """GEN.973: the fixed point already lets the vapour follow the temperature (the ocean holds it at a share of saturation at each step's own temperature, and the
    column's humidity follows Manabe and Wetherald). Hold the air of the unchanged twin and ask the table for each flux of the sweep: the response is smaller than the
    fixed point's, so the vapour feedback is in it and the held-air response is the Planck one with the lapse rate."""
    flux, airs = sweep
    held = airs[1.0].partials_pa
    held_response = (atmosphere_state.greenhouse_warming(held, twin_gravity, flux * FACTORS[-1]).surface_temperature_k
                     - atmosphere_state.greenhouse_warming(held, twin_gravity, flux * FACTORS[0]).surface_temperature_k) / (flux * (FACTORS[-1] - FACTORS[0]))
    free_response = (airs[FACTORS[-1]].temperature_k - airs[FACTORS[0]].temperature_k) / (flux * (FACTORS[-1] - FACTORS[0]))
    assert 0.0 < held_response < free_response
    assert airs[FACTORS[-1]].partials_pa["H2O"] > 1.5 * airs[FACTORS[0]].partials_pa["H2O"]          # the vapour that carried the extra warming
