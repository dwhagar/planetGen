# tests/test_surface_air_knee.py

"""
GEN.1344: an air whose fixed point does not settle stores the knee its answers chatter across, not its last answer.

Where the step finds a temperature that jumps down as the temperature it starts from rises (a haze or an ocean that switches on: the hotter air is the cooler answer), the
map has no fixed point; the damped iterate stays at the jump and its answers swing to either side. On train 15 a body at 205 K answered 350 K and 205 K in turn (the pressure
10 bar and 1.9 bar), and 4 to 5 percent of all airs (93 of 1,841 on main, 71 of 1,614 on train 15) do not settle in `MAX_STEPS`. The stored value is the mean of the answers
over the last cycles (the time average of a switch, `surface_air.sliding_temperature`), with the greenhouse factor taken from it; the air stays flagged unsettled.
"""

import dataclasses
import math

import pytest

from planetgen.physics import atmosphere, surface_air as sa
from tests.test_surface_air import _earth_like_planet


def test_the_answers_of_a_chatter_across_a_jump_give_their_mean():
    window = [(250.0 + (0.5 if k % 2 else -0.5), 320.0 if k % 2 else 180.0) for k in range(sa.SLIDING_WINDOW)]
    assert sa.sliding_temperature(window, 180.0) == (250.0, True)


def test_a_cycle_of_three_gives_the_mean_of_its_answers():
    pattern = [(249.5, 300.0), (250.0, 300.0), (250.5, 150.0)]
    window = [pattern[k % 3] for k in range(sa.SLIDING_WINDOW)]
    assert sa.sliding_temperature(window, 150.0) == (pytest.approx(250.0), True)


def test_a_creep_a_drift_and_noise_keep_the_last_answer():
    creep = [(200.0 + 0.1 * k, 200.1 + 0.1 * k) for k in range(sa.SLIDING_WINDOW)]
    assert sa.sliding_temperature(creep, 201.3) == (201.3, False)
    drift = [(200.0 + 5.0 * k, 205.0 + 5.0 * k) for k in range(sa.SLIDING_WINDOW)]              # the iterate moves as far as the answers swing: not a chatter
    assert sa.sliding_temperature(drift, 260.0) == (260.0, False)
    noise = [(250.0, 250.0 + (0.3 if k % 2 else -0.3)) for k in range(sa.SLIDING_WINDOW)]        # a swing under 10 times the tolerance
    assert sa.sliding_temperature(noise, 249.7) == (249.7, False)
    assert sa.sliding_temperature(creep[:5], 201.3) == (201.3, False)


def test_a_planet_whose_greenhouse_answers_alternate_stores_the_knee_and_a_factor_that_makes_it(monkeypatch):
    real = sa.atmosphere_state.greenhouse_warming
    answers = [320.0, 180.0]
    calls = []

    def alternating(*args):
        calls.append(1)
        return dataclasses.replace(real(*args), surface_temperature_k=answers[len(calls) % 2])

    monkeypatch.setattr(sa.atmosphere_state, "greenhouse_warming", alternating)
    planet = _earth_like_planet(3)
    assert planet.atmospheric_pressure > 0.0
    assert planet.surface_temperature == pytest.approx(250.0, abs=10.0)          # the mean of 320 and 180, not the last answer
    # the stored factor and albedo make the stored temperature (GEN.412): T = T_e (1 + G)^(1/4)
    star_flux = planet.star.luminosity / (4.0 * math.pi * (planet._star_distance_au * 1.495978707e11) ** 2)
    equilibrium = (star_flux * (1.0 - planet.bond_albedo) / (4.0 * 5.670374419e-8)) ** 0.25
    assert planet.surface_temperature == pytest.approx(equilibrium * (1.0 + planet.greenhouse_factor) ** 0.25, rel=1e-6)
