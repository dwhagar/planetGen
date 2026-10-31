# planetgen/galaxy/stepper.py

"""
The star stepper: kick-drift-kick in the galactic potential, a step per object, on arrays (GEN.340; docs/design/orbits/orbital-solvers-and-integrators.md section 5)
==================================================================================================================================================================

Replaces the turn about the galactic pole (`store._galactic_turn`), which is exact only for a planar circular orbit, with the second-order symplectic
kick-drift-kick step through `galaxy.acceleration`. Everything is written on `(n, 3)` arrays of positions in parsecs and velocities in km/s, so a million stars
are one call (the NumPy batching GEN.212 asks for), and the same inputs give the same numbers.

- **The step.** `crossing_time_myr` is a star's own timescale, the smaller of `r / v` (the time to cross its distance from the centre) and `sqrt(r / a)` (the
  free-fall time of the pull on it); `step_myr` is `STEP_ETA` (0.05) of it, no more than `STEP_CAP_MYR` (5 Myr). An orbit through its pericentre is the
  fastest part of the path, and the timescale at the star's present point is the cheap stand-in for the pericentre crossing time (registry row ST1).
- **Block steps.** A span of time is cut for each star into `2^k` equal steps, the smallest `k` for which a step is within the star's own, so the stars of one
  level go through the same number of steps together (`advance` loops over the levels, not the stars). `MAX_LEVEL` bounds `k`: a star that would need more steps
  than `2^MAX_LEVEL` is stepped at that cap and comes back flagged (`capped`), as the orbital update records a step cap.
- **An extra pull.** `advance(extra=...)` adds a callable's acceleration, `extra(positions, velocities, time_myr)` in (km/s)^2 per pc, to the same kick: the
  influence set of GEN.109 is added this way and the stepper knows nothing of it. A star with no encounter takes the potential's pull alone.
- **Time.** A Myr is `potential.KMS_TO_PC_PER_MYR` (1.0227) pc per (km/s).

Not here: reading and writing the stars of a sector, refiling a star that crosses a sector's edge, and the next-update value from the new speed (the one core
function of GEN.105 does each); the rotation about the pole stays for a facility with no vertical offset.
"""

import numpy as np

from planetgen.galaxy import potential as galactic

STEP_ETA = 0.05
"""float: A step is this share of the star's own crossing time (registry row ST1; the number is GEN.340's own)."""
STEP_CAP_MYR = 5.0
"""float: The longest step, Myr (registry row ST1; GEN.340's own)."""
MAX_LEVEL = 12
"""int: The most halvings of the cap one star may need, so at most 4096 steps a span; past it the star is stepped at the cap and flagged (registry row ST2)."""
TINY = 1.0e-300


def _norm(vectors):
    return np.sqrt(np.einsum("ij,ij->i", vectors, vectors))


def _pull(pot, positions):
    ax, ay, az = galactic.acceleration(pot, positions[:, 0], positions[:, 1], positions[:, 2])
    return np.stack([ax, ay, az], axis=1)


def crossing_time_myr(positions, velocities, pot):
    """
    Each star's own timescale, Myr: the smaller of `r / v` and `sqrt(r / a)` at its present point (`r` its distance from the centre, `v` its speed, `a` the pull on
    it). A star at rest at the centre, where both are 0 or infinite, takes the longest timescale the arithmetic can hold (the cap then limits it).
    """
    positions = np.atleast_2d(np.asarray(positions, dtype=float))
    velocities = np.atleast_2d(np.asarray(velocities, dtype=float))
    distance = _norm(positions)
    speed = _norm(velocities)
    pull = _norm(_pull(pot, positions))
    with np.errstate(divide="ignore", invalid="ignore"):
        by_speed = np.where(speed > 0.0, distance / np.maximum(speed, TINY), np.inf)
        by_pull = np.where(pull > 0.0, np.sqrt(distance / np.maximum(pull, TINY)), np.inf)
    crossing = np.minimum(by_speed, by_pull)
    # r / v is in pc per (km/s), sqrt(r / a) too: a Myr is KMS_TO_PC_PER_MYR of them.
    return np.where(np.isfinite(crossing), crossing / galactic.KMS_TO_PC_PER_MYR, np.inf)


def step_myr(positions, velocities, pot, eta=STEP_ETA, cap_myr=STEP_CAP_MYR):
    """Each star's step, Myr: `eta` of its crossing time, no more than `cap_myr`."""
    return np.minimum(eta * crossing_time_myr(positions, velocities, pot), cap_myr)


def block_levels(wanted_myr, span_myr, max_level=MAX_LEVEL):
    """
    `(levels, capped)`: for each star the `k` that cuts `span_myr` into `2^k` equal steps of at most `wanted_myr` (0 when the span is within the step), and
    whether the star needed more than `max_level` halvings (its level is then `max_level`).
    """
    wanted = np.asarray(wanted_myr, dtype=float)
    ratio = np.maximum(span_myr / np.maximum(wanted, TINY), 1.0)
    levels = np.ceil(np.log2(ratio)).astype(int)
    capped = levels > max_level
    return np.minimum(levels, max_level), capped


def kdk(positions, velocities, step_myr_, pot, extra=None, time_myr=0.0):
    """
    One kick-drift-kick step of `step_myr_` Myr for every star: half a kick, a drift, half a kick, with the pull of `pot` plus `extra(positions, velocities,
    time_myr)` (a callable, or `None`). Returns new `(positions, velocities)` arrays; the inputs are not changed.
    """
    positions = np.array(positions, dtype=float)
    velocities = np.array(velocities, dtype=float)
    dt = step_myr_ * galactic.KMS_TO_PC_PER_MYR

    def pull_at(x, v, t):
        a = _pull(pot, x)
        return a if extra is None else a + np.asarray(extra(x, v, t), dtype=float)

    velocities = velocities + 0.5 * dt * pull_at(positions, velocities, time_myr)
    positions = positions + dt * velocities
    velocities = velocities + 0.5 * dt * pull_at(positions, velocities, time_myr + step_myr_)
    return positions, velocities


def advance(positions, velocities, span_myr, pot, extra=None, eta=STEP_ETA, cap_myr=STEP_CAP_MYR, max_level=MAX_LEVEL, time_myr=0.0):
    """
    Moves every star `span_myr` Myr along its path: each star's steps are the block step its own timescale asks for (`step_myr`, `block_levels`), and the stars
    of one level go through their steps together. Returns `(positions, velocities, capped)`, `capped` a boolean array for the stars whose step was held at
    the cap `2^max_level` steps (their answer is the best the cap allows). A span of 0 returns the stars as they were.
    """
    positions = np.array(positions, dtype=float).reshape(-1, 3)
    velocities = np.array(velocities, dtype=float).reshape(-1, 3)
    capped = np.zeros(len(positions), dtype=bool)
    if span_myr < 0.0:
        raise ValueError("the span must not be negative")
    if span_myr == 0.0 or not len(positions):
        return positions, velocities, capped
    wanted = step_myr(positions, velocities, pot, eta, cap_myr)
    levels, capped = block_levels(wanted, span_myr, max_level)
    # `wanted` is no more than the cap, so no step is longer than the cap.
    for level in np.unique(levels):
        picked = np.flatnonzero(levels == level)
        count = 2 ** int(level)
        step = span_myr / count
        x, v = positions[picked], velocities[picked]
        for n in range(count):
            x, v = kdk(x, v, step, pot, extra, time_myr + n * step)
        positions[picked], velocities[picked] = x, v
    return positions, velocities, capped
