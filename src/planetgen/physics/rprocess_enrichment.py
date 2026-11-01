# planetgen/physics/rprocess_enrichment.py

"""
The thorium-to-uranium enrichment of a star's gas from the disc's history (GEN.663, docs/design/composition/radioisotope-tuning.md section 2)
===========================================================================================================================================

Row o4 of the atmosphere recipes needs the factor `F` of the thorium-to-uranium ratio of the gas a star formed from:

    Th/U (mass) = 0.97483 * PR * F * exp[(l238 - l232) * age]

with `PR` the production ratio of the r-process (1.667 to 2.0, [R1447]) and `age` the time since the star formed. `F` is
what the decay of the earlier ejecta did to the ratio before the star formed: thorium-232 (half-life 14.05 Gyr) outlives
uranium-238 (4.468 Gyr), so gas enriched long before the star formed is thorium-rich.

The row used a constant r-process input over a duration `T` of 3 to 13 Gyr, calibrated on one Solar System datum. This module
takes the input from the disc's own history instead:

1. **The sources.** The r-process input at time `t` is the star formation `SFR(t - tau)` convolved with a delay-time law
   `1/tau` from `TAU_MIN_GY` to the age of the disc (neutron star mergers, [R2010, R4101]). The slope 1
   and the 30 Myr floor are RECALLED, a stand-in: the delay law is the one number that shapes `F` and is labelled in the gap
   plan registry. Setting `TAU_MIN_GY` to nothing larger than a time step gives prompt sources (collapsars).
2. **The star formation** is the disc's (`galaxy.disk_history`, GEN.469), so a galaxy that formed its stars early gives thorium-rich
   gas to the stars that form late and a galaxy that forms them now does not.
3. **The decay.** Ejecta of time `t'` have decayed by `exp(-l (t_b - t'))` when a star forms at `t_b`; `F` is the ratio of the
   thorium-232 to the uranium-238 that is left, summed over the input.

Left out: the gas before the disc began (the halo and thick disc ejecta), so the factor is a lower bound for old stars; and the
variation with galactocentric radius. The Solar System datum (CI chondrites at about 4.5 Gyr need an effective constant input of
6.5 Gyr) is a test point and not an input.
"""

import functools
import math

import numpy as np

from planetgen.physics import elements

TAU_MIN_GY = 0.03
"""float: The shortest merger delay, Gyr (RECALLED stand-in; see the module note)."""

DELAY_NODES = 160
"""int: Log-spaced delay nodes of the convolution."""

SECONDS_PER_GYR = 3.15576e16
"""float: Seconds in a Gyr (Julian year)."""

PRODUCTION_RATIO_RANGE = (1.667, 2.0)
"""tuple: The r-process production ratio of thorium-232 to uranium-238 by number, the range of the theoretical predictions
([R1447: Table 1], read in docs/design/composition/radioisotope-tuning.md section 2); a body draws one value from it."""


def decay_constants_per_gy():
    """`(l232, l238)`: the decay constants of thorium-232 and uranium-238 per Gyr, from the isotope table."""
    return tuple(math.log(2.0) / (elements.isotope(name).half_life_s / SECONDS_PER_GYR) for name in ("Th-232", "U-238"))


def constant_input_factor(duration_gy):
    """The factor `F` of a constant r-process input over `duration_gy` before the star formed (the row's former form)."""
    l232, l238 = decay_constants_per_gy()
    return ((1.0 - math.exp(-l232 * duration_gy)) / l232) / ((1.0 - math.exp(-l238 * duration_gy)) / l238)


def effective_duration_gy(factor):
    """The duration of a constant input that gives `factor`, Gyr, by bisection (`factor` is 1 at 0 and rises with duration)."""
    low, high = 0.0, 60.0
    for _ in range(80):
        middle = 0.5 * (low + high)
        if constant_input_factor(middle) < factor:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def source_rate(times_gy, star_formation, tau_min_gy=None):
    """
    The r-process input at each time: `star_formation` (on `times_gy`, increasing from 0) convolved with the delay-time law
    1/tau from `tau_min_gy` to the span of the times. Arbitrary units.
    """
    times = np.asarray(times_gy, dtype=float)
    sfr = np.asarray(star_formation, dtype=float)
    tau_min = TAU_MIN_GY if tau_min_gy is None else tau_min_gy
    span = float(times[-1])
    taus = np.geomspace(tau_min, span, DELAY_NODES)
    # d tau / tau = d ln tau: trapezoids on the log grid
    weights = np.full(DELAY_NODES, math.log(span / tau_min) / (DELAY_NODES - 1))
    weights[0] *= 0.5
    weights[-1] *= 0.5
    out = np.zeros(len(times))
    for tau, weight in zip(taus, weights):
        delayed = times - tau
        out += weight * np.where(delayed >= 0.0, np.interp(np.maximum(delayed, 0.0), times, sfr), 0.0)
    return out


def th_u_factor(times_gy, star_formation, birth_gy, tau_min_gy=None):
    """
    The factor `F` of a star born `birth_gy` after the disc began.

    Args:
        times_gy (sequence): Increasing times since the disc began, from 0 to its age.
        star_formation (sequence): The disc's star formation at those times.
        birth_gy (float): When the star formed, on the same clock.
        tau_min_gy (float, optional): The shortest merger delay.

    Returns:
        float: `F`, at least 1 and at most that of an input of the whole disc life. Zero birth time (no earlier input) gives 1.
    """
    times = np.asarray(times_gy, dtype=float)
    birth = min(max(float(birth_gy), 0.0), float(times[-1]))
    if birth <= 0.0:
        return 1.0
    rate = source_rate(times, star_formation, tau_min_gy)
    grid = np.linspace(0.0, birth, 2001)
    input_rate = np.interp(grid, times, rate)
    l232, l238 = decay_constants_per_gy()
    thorium = np.trapezoid(input_rate * np.exp(-l232 * (birth - grid)), grid)
    uranium = np.trapezoid(input_rate * np.exp(-l238 * (birth - grid)), grid)
    if not uranium > 0.0:
        return 1.0
    return float(thorium / uranium)


def th_u_factor_of_history(history, host_age_gy, tau_min_gy=None):
    """
    The factor `F` of a star of age `host_age_gy` in the disc of a `galaxy.disk_history.DiskHistory`: the disc began
    `history.age_gy[0]` ago, and the star formed `host_age_gy` ago (older than the disc is held at the disc's start).
    """
    span = float(history.age_gy[0])
    times = span - np.asarray(history.age_gy, dtype=float)
    star_formation = np.sum(history.radius_scale_lengths[:, None] * history.rate, axis=0)
    return th_u_factor(times, star_formation, span - host_age_gy, tau_min_gy)


@functools.lru_cache(maxsize=4096)
def th_u_factor_of_shape(shape, host_age_gy):
    """The factor `F` of a star of age `host_age_gy` in the thin disc of the galaxy of this `GalaxyShape` (`th_u_factor_of_history` of its
    `galaxy.disk_history.history_of`), cached because every planet and moon of a star asks for it."""
    from planetgen.galaxy import disk_history
    return th_u_factor_of_history(disk_history.history_of(shape), host_age_gy)


@functools.lru_cache(maxsize=64)
def typical_th_u_factor(shape):
    """The factor `F` of a typical star of the galaxy of this shape, for a body with no host star (a rogue planet): that of the star born when
    half of the disc's stars had formed (the median of its star formation, the rings summed as `th_u_factor_of_history` does)."""
    from planetgen.galaxy import disk_history
    history = disk_history.history_of(shape)
    span = float(history.age_gy[0])
    times = span - np.asarray(history.age_gy, dtype=float)
    rate = np.sum(history.radius_scale_lengths[:, None] * history.rate, axis=0)
    cumulative = np.concatenate(([0.0], np.cumsum(0.5 * (rate[1:] + rate[:-1]) * np.diff(times))))
    median = float(np.interp(0.5 * cumulative[-1], cumulative, times)) if cumulative[-1] > 0.0 else 0.5 * span
    return th_u_factor(times, rate, median)
