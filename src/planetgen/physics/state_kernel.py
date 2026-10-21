# planetgen/physics/state_kernel.py

"""
The compiled state-vector step (GEN.289)
========================================

`kepler.universal_step` for a whole table of planets or moons in one compiled call. The arithmetic is the Python
function's, line for line (the Stumpff series near zero, Halley's method kept inside a bracket that always holds the
root, bisection when a step leaves it, a whole number of orbits dropped before the bracket), so the two agree to a few
parts in 1e15; `tests/test_state_kernel.py` holds them to 1e-12.

Units are the caller's, so long as they are consistent: `physics/body_state.py` passes kilometres and seconds, with
`mu` in km^3/s^2.

Like `field_kernel`, the module is **kept out of the web app's import path** (numba's cache needs a writable code
tree): the nightly orbit update imports it inside the function that calls it, and
`test_the_web_app_never_imports_numba` holds that.
"""

import math

import numpy as np
from numba import njit

from planetgen.physics import kepler

TOLERANCE = 1e-13
"""float: `kepler.UNIVERSAL_TOLERANCE`."""

MAX_ITERATIONS = 60
"""int: `kepler.UNIVERSAL_MAX_ITERATIONS`."""

OK = 0
"""int: The status of a row that was stepped."""

BAD_INPUT = 1
"""int: The status of a row with a non-positive `mu`, a non-finite input or a body on its primary."""

NO_CONVERGENCE = 2
"""int: The status of a row whose iteration reached `MAX_ITERATIONS`."""


_C2_COEFFICIENTS = tuple(1.0 / math.factorial(2 * k + 2) for k in range(kepler.STUMPFF_SERIES_TERMS))
_C3_COEFFICIENTS = tuple(1.0 / math.factorial(2 * k + 3) for k in range(kepler.STUMPFF_SERIES_TERMS))
_SERIES_LIMIT = kepler.STUMPFF_SERIES_LIMIT
_SERIES_TERMS = kepler.STUMPFF_SERIES_TERMS


@njit(cache=True)
def _stumpff(z):
    """`kepler.stumpff_c2_c3`, compiled: the series below `kepler.STUMPFF_SERIES_LIMIT`, the closed forms beyond."""
    if abs(z) < _SERIES_LIMIT:
        c2 = 0.0
        c3 = 0.0
        for k in range(_SERIES_TERMS - 1, -1, -1):
            c2 = c2 * -z + _C2_COEFFICIENTS[k]
            c3 = c3 * -z + _C3_COEFFICIENTS[k]
        return c2, c3
    if z > 0:
        root = math.sqrt(z)
        return (1 - math.cos(root)) / z, (root - math.sin(root)) / (root ** 3)
    root = math.sqrt(-z)
    return (math.cosh(root) - 1) / -z, (math.sinh(root) - root) / (root ** 3)


@njit(cache=True)
def _residual(chi, s0, alpha, r0, sqrt_mu, dt):
    z = alpha * chi * chi
    c2, c3 = _stumpff(z)
    f = s0 * chi * chi * c2 + (1 - alpha * r0) * chi ** 3 * c3 + r0 * chi - sqrt_mu * dt
    df = s0 * chi * (1 - z * c3) + (1 - alpha * r0) * chi * chi * c2 + r0
    d2f = s0 * (1 - z * c2) + (1 - alpha * r0) * chi * (1 - z * c3)
    return f, df, d2f


@njit(cache=True)
def step_one(state, mu, dt, out):
    """
    Moves one state `dt` forward on its conic into `out` (six floats: position, then velocity) and returns the status
    (`OK`, `BAD_INPUT` or `NO_CONVERGENCE`); `out` is left as `state` unless the status is `OK`.
    """
    for k in range(6):
        out[k] = state[k]
    if not (mu > 0.0) or not math.isfinite(mu):
        return BAD_INPUT
    if not math.isfinite(dt):
        return BAD_INPUT
    for k in range(6):
        if not math.isfinite(state[k]):
            return BAD_INPUT
    if dt == 0.0:
        return OK
    r0 = math.sqrt(state[0] * state[0] + state[1] * state[1] + state[2] * state[2])
    if r0 == 0.0:
        return BAD_INPUT
    v0_sq = state[3] * state[3] + state[4] * state[4] + state[5] * state[5]
    sqrt_mu = math.sqrt(mu)
    s0 = (state[0] * state[3] + state[1] * state[4] + state[2] * state[5]) / sqrt_mu
    alpha = 2 / r0 - v0_sq / mu

    if alpha > 1e-12 / r0:
        chi = sqrt_mu * dt * alpha
        period = 2 * math.pi / (sqrt_mu * alpha ** 1.5)
        if abs(dt) > period:
            dt = np.fmod(dt, period)
            chi = sqrt_mu * dt * alpha
    elif alpha < -1e-12 / r0:
        a = 1 / alpha
        sign = 1.0 if dt > 0 else -1.0
        inner = (-2 * mu * alpha * dt) / (s0 * sqrt_mu + sign * math.sqrt(-mu * a) * (1 - r0 * alpha))
        if inner > 0:
            chi = sign * math.sqrt(-a) * math.log(inner)
        else:
            chi = sign * math.sqrt(-a)
    else:
        chi = sqrt_mu * dt / r0

    guess = abs(chi) if abs(chi) != 0.0 else 1.0
    if dt > 0:
        low, high = 0.0, guess
    else:
        low, high = -guess, 0.0
    for _ in range(200):
        if dt > 0 and _residual(high, s0, alpha, r0, sqrt_mu, dt)[0] < 0:
            low, high = high, high * 2
        elif dt < 0 and _residual(low, s0, alpha, r0, sqrt_mu, dt)[0] > 0:
            low, high = low * 2, low
        else:
            break
    if not low <= chi <= high:
        chi = 0.5 * (low + high)

    converged = False
    for _ in range(MAX_ITERATIONS):
        f, df, d2f = _residual(chi, s0, alpha, r0, sqrt_mu, dt)
        if f > 0:
            high = chi
        else:
            low = chi
        denominator = 2 * df * df - f * d2f
        if denominator != 0:
            step = 2 * f * df / denominator
        else:
            step = f / df
        new = chi - step
        if not low < new < high:
            new = 0.5 * (low + high)
        if abs(new - chi) <= TOLERANCE * max(1.0, abs(new)) or high - low <= 4e-16 * max(1.0, abs(new)):
            chi = new
            converged = True
            break
        chi = new
    if not converged:
        return NO_CONVERGENCE

    z = alpha * chi * chi
    c2, c3 = _stumpff(z)
    f = 1 - chi * chi * c2 / r0
    g = dt - chi ** 3 * c3 / sqrt_mu
    rx = f * state[0] + g * state[3]
    ry = f * state[1] + g * state[4]
    rz = f * state[2] + g * state[5]
    r = math.sqrt(rx * rx + ry * ry + rz * rz)
    g_dot = 1 - chi * chi * c2 / r
    f_dot = sqrt_mu * chi * (z * c3 - 1) / (r * r0)
    out[0] = rx
    out[1] = ry
    out[2] = rz
    out[3] = f_dot * state[0] + g_dot * state[3]
    out[4] = f_dot * state[1] + g_dot * state[4]
    out[5] = f_dot * state[2] + g_dot * state[5]
    return OK


@njit(cache=True)
def _step_many(states, mus, dts, out, status):
    for i in range(states.shape[0]):
        status[i] = step_one(states[i], mus[i], dts[i], out[i])


def step_many(states, mus, dts):
    """
    Steps every row of `states` (an (n, 6) array: position, then velocity) by its own `dts` on its own `mus`.

    Returns:
        tuple: `(states, status)`, the stepped (n, 6) array (a row whose status is not `OK` is unchanged) and an (n,)
            array of `OK`, `BAD_INPUT` or `NO_CONVERGENCE`.
    """
    states = np.ascontiguousarray(states, dtype=np.float64)
    mus = np.ascontiguousarray(mus, dtype=np.float64)
    dts = np.ascontiguousarray(dts, dtype=np.float64)
    out = np.empty_like(states)
    status = np.empty(states.shape[0], dtype=np.int64)
    _step_many(states, mus, dts, out, status)
    return out, status
