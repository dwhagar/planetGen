"""
Climate sensitivity sweep (GEN.975): the surface-air fixed point of an Earth-like twin, solved again at a range of absorbed sunlight from a cooling
to a warming sweep. Every solve must converge, the temperature must rise with the light, and the response to a change of the absorbed flux must stay
inside the published range of Earth's climate sensitivity: an equilibrium warming of 2 to 5 K for a doubling of carbon dioxide (the very likely
range of IPCC AR6, 2021, chapter 7), a forcing of 3.7 W/m2 for the doubling, so 0.54 to 1.35 K per W/m2 of absorbed flux.

Where the fixed point stands: it converges and warms with the light, with a response (about 0.4 K per W/m2 with the old grey closure) under the published lower
bound, because the feedbacks that make a climate more sensitive than its Planck response are not in it yet. The lower bound is held by a strict expected
failure that names them; when GEN.972 (ice-albedo), GEN.973 (water vapour) and GEN.974 (cloud albedo) bring the response inside the range, the test passes, the
expected failure fails the run, and the mark comes off.

The body is Earth-like and is a test target, never an input of the code. The solve is replayed on one stream (`draw.bound`), so the sweep changes the
light and nothing else; the albedo is the body's own, held, until the feedbacks that move it (GEN.972 ice-albedo, GEN.974 cloud albedo) are in the fixed point.
"""

import dataclasses

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.physics import constants, surface_air as sa
from planetgen.util import draw

SENSITIVITY_K_PER_W_M2 = (2.0 / 3.7, 5.0 / 3.7)
"""tuple: The very likely range of the equilibrium climate sensitivity (2 to 5 K per doubling of CO2, AR6) over the forcing of a doubling (3.7 W/m2)."""
EARTH_WATER_KG, EARTH_OCEAN_COVER = 1.4e21, 0.71
"""float: Earth's ocean (a test target): the twin's own water is a dry planet's, and a sweep through it crosses the freezing of its small sea."""
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
    with draw.bound(3):
        Planet(config, star, star.habitable_zone, 0.5 * (inner + outer), mass=constants.EARTH_MASS_TO_KG, moon_count=0)
    monkeypatch.undo()
    return (dataclasses.replace(captured["setting"], water_kg=EARTH_WATER_KG, ocean_cover=EARTH_OCEAN_COVER), captured["start"])


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
