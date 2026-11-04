# planetgen/galaxy/rogue_density.py

"""
The smoothed star density a rogue planet is drawn from (GEN.1115; docs/design/systems/rogue-planet-drift.md sections 3.3, 5 and 7)
====================================================================================================================================

A rogue planet is nowhere near the star it left. Thrown out with a few km/s, it slides along the galactic orbit by about a kiloparsec every half Gyr (`physics/rogue_drift.py`), so an old rogue has lost the donor's
azimuth: its density depends on galactocentric radius and height, not on which spiral arm, cluster or bar end the donor was in. The number of rogues in a sector is therefore the yield per star times the sector's
volume times a *smoothed* star density, not the sector's own. This module is that smoothing, a pure function of any star-density function, shared by the sector draw (the Poisson mean of a sector) and the low-zoom
map (which can read the smooth density in place of drawing rogues one by one far out, so the zoom levels agree).

1. **Ring average.** The azimuthal mean of the density at the sector's radius `R` and height `z`: `RING_STEPS` points at equal angles, starting half a step off the axes. A pattern with `m` arms or bar ends is
   averaged out exactly for `m` below `RING_STEPS` (the mean of `cos(m (phi - phi0))` over the points is zero).
2. **Radial and vertical smoothing.** The donor was not on this ring and this height: the guiding radius shifted by `2 Omega v / kappa^2` and the epicycle and the vertical oscillation carry it a bounded way more. The density is
   averaged over a Gaussian of `RADIAL_SMOOTH_PC` in `R` (reflected at the axis) and `VERTICAL_SMOOTH_PC` in `z`, three nodes each (the 3-point Gauss-Hermite rule of a normal distribution, exact for a polynomial of degree 5).
   The two widths are the rms radial and vertical offsets `rogue_drift.drift` gives for the stand-in kick distribution (median 3 km/s, ln-sigma 0.8, isotropic) at 1 to 10 Gyr (193 and 30 pc; `tests/test_rogue_density.py`
   recomputes them), a Gaussian standing for a tail that reaches 0.7 kpc at the 99th percentile.
3. **The young share.** Rogues ejected within `YOUNG_AGE_MYR` have not lost the donor's structure: they follow the density at the sector itself. With ejections uniform over the population's history, that share is
   `YOUNG_AGE_MYR / history_myr` (a few percent); the rest follows the smoothed density.

The density function is the caller's: any `density_at(x, y, z)` that takes numpy arrays of galaxy-frame parsecs and returns the density (or the relative density) at each, for example
`lambda x, y, z: galaxy.density.total_density_array(x, y, z, shape)`. All `RING_STEPS` x 3 x 3 points are asked for in one call (so a vectorised function costs about one evaluation of 144 points).

Use:

    old = old_density(density_at, radius_pc, z_pc)                               # the smooth part, a function of (R, z) only: what the low-zoom map reads
    mean = rogue_density(density_at, (x, y, z), history_myr)                     # the density the sector's Poisson mean takes: young share local, rest smooth

Limits. A stationary disc is assumed (the ejection-time density equals today's; radial migration and the disc's ageing break that at the 10 to 30 percent level, the design note); the metallicity gradient across the
displacement is ignored; nothing here knows the rogue yield or the sector volume, which are the caller's (`lambda = yield x V x density`).
"""

import math

import numpy as np

from planetgen.util import detmath
from planetgen.util.checks import finite_domain

RING_STEPS = 16
"""int: Points round a ring of the average (the design note's "about 16"). [Numerical, stand-in RSM4, tier T9, tweak L0.]"""

RADIAL_SMOOTH_PC = 190.0
"""float: The Gaussian width in galactocentric radius of the displacement between a donor and an old rogue: the rms radial offset of `rogue_drift.drift` for the stand-in kick distribution at 1 to 10 Gyr
(the design note says "0.1 to 0.7 kpc" and guessed 0.5). [Stand-in RSM1, tier T6, tweak L1.]"""

VERTICAL_SMOOTH_PC = 30.0
"""float: The same in height (rms 30 pc; the design note says "about 20 to 50 pc"). [Stand-in RSM2, tier T6, tweak L1.]"""

YOUNG_AGE_MYR = 100.0
"""float: Rogues ejected within this long keep the donor's structure (the design note: "younger than about 100 Myr"), the share of the sector's own density. [Stand-in RSM3, tier T6, tweak L1.]"""

_NODES = (-math.sqrt(3.0), 0.0, math.sqrt(3.0))
_WEIGHTS = (1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0)
"""The 3-point Gauss-Hermite rule for a standard normal (nodes and weights)."""

_RING_ANGLES = (np.arange(RING_STEPS) + 0.5) * (2.0 * math.pi / RING_STEPS)
_COS, _SIN = np.cos(_RING_ANGLES), np.sin(_RING_ANGLES)


def _check_position(radius_pc, z_pc):
    if not (radius_pc >= 0.0 and math.isfinite(radius_pc) and math.isfinite(z_pc)):
        raise ValueError("the radius must be zero or more and the height finite")


def ring_average(density_at, radius_pc, z_pc):
    """The mean of `density_at` over the ring of radius `radius_pc` at height `z_pc` (`RING_STEPS` points)."""
    _check_position(radius_pc, z_pc)
    values = np.asarray(density_at(radius_pc * _COS, radius_pc * _SIN, np.full(RING_STEPS, float(z_pc))), dtype=float)
    return float(values.mean())


@finite_domain()
def old_density(density_at, radius_pc, z_pc, radial_smooth_pc=RADIAL_SMOOTH_PC, vertical_smooth_pc=VERTICAL_SMOOTH_PC):
    """
    The smoothed density an old rogue (one that has lost its donor's azimuth) is drawn from, at galactocentric radius `radius_pc` and height `z_pc`: the ring average, averaged over the radial and vertical
    Gaussians. It depends on the radius and the height only.

    Args:
        density_at (callable): `density_at(x, y, z)` for numpy arrays of galaxy-frame parsecs, giving the density at each point.
        radius_pc, z_pc (float): The sector's galactocentric radius and height.
        radial_smooth_pc, vertical_smooth_pc (float): The Gaussian widths (zero means no smoothing in that direction).
    """
    _check_position(radius_pc, z_pc)
    if radial_smooth_pc < 0.0 or vertical_smooth_pc < 0.0:
        raise ValueError("the smoothing widths must not be negative")
    radii = [abs(radius_pc + radial_smooth_pc * t) for t in _NODES] if radial_smooth_pc > 0.0 else [radius_pc]
    heights = [z_pc + vertical_smooth_pc * t for t in _NODES] if vertical_smooth_pc > 0.0 else [z_pc]
    radial_w = _WEIGHTS if radial_smooth_pc > 0.0 else (1.0,)
    vertical_w = _WEIGHTS if vertical_smooth_pc > 0.0 else (1.0,)
    x = np.concatenate([r * _COS for r in radii for _ in heights])
    y = np.concatenate([r * _SIN for r in radii for _ in heights])
    z = np.concatenate([np.full(RING_STEPS, h) for _ in radii for h in heights])
    values = np.asarray(density_at(x, y, z), dtype=float).reshape(len(radii) * len(heights), RING_STEPS).mean(axis=1)
    weights = np.array([wr * wz for wr in radial_w for wz in vertical_w])
    return float(detmath.matmul(weights, values[:, None])[0])


def young_share(history_myr):
    """The share of ejections that are younger than `YOUNG_AGE_MYR` when they are uniform over a population's history of `history_myr`."""
    if not history_myr > 0.0:
        raise ValueError("the history must be positive")
    return min(1.0, YOUNG_AGE_MYR / history_myr)


@finite_domain()
def rogue_density(density_at, position_pc, history_myr, radial_smooth_pc=RADIAL_SMOOTH_PC, vertical_smooth_pc=VERTICAL_SMOOTH_PC):
    """
    The density a sector's rogue count takes: `young_share` of the density at the sector itself plus the rest of the smoothed density of `old_density`.

    Args:
        density_at (callable): As for `old_density`.
        position_pc (tuple): The sector's galaxy-frame position `(x, y, z)`, pc.
        history_myr (float): How long the population has been ejecting planets (the donors' age span).
    """
    x, y, z = position_pc
    share = young_share(history_myr)
    local = float(np.asarray(density_at(np.array([x], dtype=float), np.array([y], dtype=float), np.array([z], dtype=float)), dtype=float)[0]) if share > 0.0 else 0.0
    return share * local + (1.0 - share) * old_density(density_at, detmath.hypot(x, y), z, radial_smooth_pc, vertical_smooth_pc)
