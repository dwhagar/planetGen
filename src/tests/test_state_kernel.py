# tests/test_state_kernel.py

"""GEN.289: the compiled state step agrees with `kepler.universal_step`, and the helpers around it."""

import math

import numpy as np
import pytest

from planetgen.physics import body_state, constants as pc, kepler, state_kernel
from planetgen.physics.orbits import orbital_position_au

MU = 1.327e11
"""float: The Sun's GM, km^3/s^2."""


def _circular(distance_km, mu=MU):
    return (distance_km, 0.0, 0.0), (0.0, math.sqrt(mu / distance_km), 0.0)


def _states():
    cases = []
    for distance in (5.8e7, 1.5e8, 7.8e8, 4.5e9):
        position, velocity = _circular(distance)
        cases.append(position + velocity)
    cases.append((1.5e8, 0.0, 0.0, 0.0, 40.0, 3.0))  # eccentric and inclined
    cases.append((1.5e8, 2.0e7, -1.0e6, -5.0, 29.0, 1.0))
    cases.append((1.5e8, 0.0, 0.0, 0.0, 12.0, 0.0))  # a long thin ellipse
    cases.append((1.5e8, 0.0, 0.0, 0.0, 45.0, 0.0))  # open: a hyperbola
    return cases


@pytest.mark.parametrize("dt", [0.0, 1.0e3, 1.0e7, 3.1557600e7, -1.0e6, 3.1557600e10])
def test_the_compiled_step_matches_the_python_step(dt):
    states = np.array(_states())
    moved, status = state_kernel.step_many(states, np.full(len(states), MU), np.full(len(states), dt))
    assert (status == state_kernel.OK).all()
    for state, got in zip(states, moved):
        position, velocity = kepler.universal_step(tuple(state[:3]), tuple(state[3:]), MU, dt)
        expected = np.array(position + velocity)
        scale = max(1.0, float(np.abs(expected[:3]).max()))
        assert np.abs(got[:3] - expected[:3]).max() <= 1e-12 * scale
        assert np.abs(got[3:] - expected[3:]).max() <= 1e-12 * max(1.0, float(np.abs(expected[3:]).max()))


def test_a_bad_row_keeps_its_state_and_does_not_stop_the_others():
    good = _circular(1.5e8)
    states = np.array([good[0] + good[1], (0.0, 0.0, 0.0, 1.0, 0.0, 0.0), (float("nan"), 0.0, 0.0, 0.0, 1.0, 0.0),
                       good[0] + good[1]])
    mus = np.array([MU, MU, MU, -1.0])
    moved, status = state_kernel.step_many(states, mus, np.full(4, 1.0e6))
    assert list(status) == [state_kernel.OK, state_kernel.BAD_INPUT, state_kernel.BAD_INPUT, state_kernel.BAD_INPUT]
    assert moved[0][1] != states[0][1]
    assert np.array_equal(moved[1], states[1])
    assert np.array_equal(moved[3], states[3])


def test_a_circular_state_stays_on_its_circle_over_a_thousand_years():
    position, velocity = _circular(1.5e8)
    period = 2 * math.pi * math.sqrt(1.5e8 ** 3 / MU)
    states = np.array([position + velocity])
    moved, status = state_kernel.step_many(states, np.array([MU]), np.array([1000.0 * pc.SECONDS_PER_YEAR]))
    assert status[0] == state_kernel.OK
    assert math.hypot(*moved[0][:3]) == pytest.approx(1.5e8, rel=1e-9)
    assert math.hypot(*moved[0][3:]) == pytest.approx(math.sqrt(MU / 1.5e8), rel=1e-9)
    turn = (1000.0 * pc.SECONDS_PER_YEAR % period) / period * 2 * math.pi
    assert math.atan2(moved[0][1], moved[0][0]) % (2 * math.pi) == pytest.approx(turn % (2 * math.pi), abs=1e-6)


def test_mu_is_keplers_third_law():
    period_years = 2 * math.pi * math.sqrt(1.5e8 ** 3 / MU) / pc.SECONDS_PER_YEAR
    assert body_state.mu_km3_s2(1.5e8, period_years) == pytest.approx(MU, rel=1e-12)


@pytest.mark.parametrize("inclination,node,phase", [(0.0, 0.0, 10.0), (5.0, 80.0, 200.0), (30.0, 300.0, 359.0),
                                                   (89.0, 12.0, 91.0), (3.0, 180.0, 0.5)])
def test_the_argument_of_latitude_reads_back_the_phase_a_position_was_made_with(inclination, node, phase):
    position = tuple(a * pc.AU_TO_KM for a in orbital_position_au(1.0, inclination, node, phase))
    assert body_state.argument_of_latitude_deg(position, inclination, node) == pytest.approx(phase, abs=1e-9)


def test_advance_moves_a_body_a_whole_orbit_back_to_where_it_started():
    position, velocity = _circular(1.5e8)
    period_years = 2 * math.pi * math.sqrt(1.5e8 ** 3 / MU) / pc.SECONDS_PER_YEAR
    got, _ = body_state.advance(position, velocity, 1.5e8, period_years, 7 * period_years * pc.SECONDS_PER_YEAR)
    assert got == pytest.approx(position, abs=1e-3)
