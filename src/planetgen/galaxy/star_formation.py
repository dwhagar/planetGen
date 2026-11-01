# planetgen/galaxy/star_formation.py

"""
Recent Star Formation by Galactic Radius
========================================

How much of the thin disk's young and intermediate stars (the ones formed in the last few billion years) a place
holds, against the disk's own average (GEN.134, `docs/design/galaxy/star-types-by-galactic-radius.md`). Derived from the
galaxy's own shape, with no value taken from the Milky Way's star-formation profile (those are the test targets of
`tests/test_star_formation_profile.py`, never inputs):

1. The star-formation law [R2321]: the rate per area follows the gas per area to a power `n` (1.4). The
   thin disk's stars are the sum of its past rates, so their surface density is `exp(-R / L)` (`L` the disk scale
   length). With a gas disk of scale length `g * L` the present rate falls as `exp(-n R / (g L))`, and the share
   of the stars that are recent as `exp(-(n / g - 1) R / L)`: the outer disk has fewer young stars per star.
   Where the gas thins the index rises (GEN.450): molecular gas forms stars at a constant efficiency and atomic gas,
   which dominates where the gas is under about 9 solar masses per pc2, at a falling one ([R0392], Leroy et
   al. 2008, [R3951]), so the rate is `A gas^n (gas / (gas + gas_s))^(m - n)` with `gas_s` the saturation density
   (`tuning.STAR_FORMATION_SATURATION_GAS_MSUN_PC2`) and `m` the steeper index
   (`tuning.STAR_FORMATION_LOW_DENSITY_EXPONENT`), and the gas at each radius is the disk's own, from the rings of its
   star-formation history (`galaxy.disk_history`). The extra factor is 1 where the gas is dense and
   `(gas / gas_s)^(m - n)` where it is thin.
2. The bar gap: a bar drives the gas inside its corotation inward and shocks and shears what is left, so the
   star-forming gas inside the bar's length is depleted (the star-formation desert of a barred galaxy; see the
   design note for the sources). The rate there is scaled by `S(R)`, a floor `f` at the centre rising to 1 outside
   the gap radius `R_gap` (a multiple of the bulge's scale radius, `tuning.BAR_GAP_RADIUS_RATIO`).
3. The disk makes as many stars over its life as before: the factor is normalised so the disk-wide average of the
   recent share is 1 (the populations still split the thin disk's density, so no cell's total changes).

The disk-wide shares the factor multiplies come from the disk's own gas history (GEN.469, `galaxy.disk_history`,
`disk_population_share`), not from the age windows' lengths.

The factor multiplies the weights of the `young` and `intermediate` populations in
`galaxy.density.population_densities`. It reads only the shape's scale lengths, so a galaxy of any size gets one.
"""

import functools

import numpy as np

from planetgen import tuning
from planetgen.galaxy import disk_history

RECENT_POPULATIONS = ("young", "intermediate")
"""tuple: The populations the factor applies to: those formed within the bar's and the disk's recent history."""

_GRID_POINTS = 6000
"""int: Points of the grid the normalisation integral is taken on (it is smooth; 6000 is converged to 1e-6)."""


def disk_population_share(name, shape):
    """
    The share of the thin disk's stars in the disk population `name` (`young`, `intermediate` or `old`) over the whole
    disk, from the disk's gas history (GEN.469, `galaxy.disk_history`): how much of its star formation fell in the
    population's age window (`tuning.STELLAR_POPULATION_AGE_RANGES_GY`). The recent factor then moves the young and
    intermediate shares from place to place; the three always sum to 1.
    """
    return disk_history.population_shares(shape)[name]


def gap_radius_pc(shape):
    """The radius inside which the bar has depleted the star-forming gas: `tuning.BAR_GAP_RADIUS_RATIO` bulge scale
    radii (the bar's own length scale)."""
    return tuning.BAR_GAP_RADIUS_RATIO * shape.bulge_scale_radius_pc


def _gap_floor(shape):
    """What is left of the star formation at the centre of the bar gap: `tuning.BAR_GAP_STAR_FORMATION_FLOOR`, or 1
    (no gap) for a galaxy without a bar (`shape.structure.bar_present`, GEN.533)."""
    return tuning.BAR_GAP_STAR_FORMATION_FLOOR if shape.structure.bar_present else 1.0


def _gap_suppression(radius_pc, gap_pc, floor):
    """`S(R)`: `floor` at the centre, half-way up at the gap radius, 1 well outside it (works on arrays)."""
    up = np.power(np.asarray(radius_pc, dtype=float) / gap_pc, tuning.BAR_GAP_EDGE_POWER)
    return floor + (1.0 - floor) * up / (1.0 + up)


def _log_gas(x, shape):
    """`(ln gas, d ln gas / dx)` at `x` disk scale lengths (arrays): the log of the disk's gas per area now (solar masses
    per pc2, the rings of `disk_history.history_of`) joined by straight lines, and past the outermost ring continued with
    the slope of the last segment (the gas of an exponential disk)."""
    found = disk_history.history_of(shape)
    rings, log_gas = found.radius_scale_lengths, np.log(found.gas_now)
    x = np.asarray(x, dtype=float)
    slopes = np.diff(log_gas) / np.diff(rings)
    segment = np.clip(np.searchsorted(rings, x) - 1, 0, len(slopes) - 1)
    slope = slopes[segment]
    return log_gas[segment] + slope * (x - rings[segment]), slope


def _thin_gas_factor(x, shape):
    """
    The factor the rising index puts on the recent share at `x` disk scale lengths, with its log slope per scale length:
    `(gas / (gas + gas_s))^(m - n)`. 1 where the gas is dense, falling as the gas thins.
    """
    log_gas, slope = _log_gas(x, shape)
    gas = np.exp(log_gas)
    saturation = tuning.STAR_FORMATION_SATURATION_GAS_MSUN_PC2
    power = tuning.STAR_FORMATION_LOW_DENSITY_EXPONENT - tuning.STAR_FORMATION_LAW_EXPONENT
    return (gas / (gas + saturation)) ** power, power * saturation / (gas + saturation) * slope


def _decline_per_scale_length():
    """`n / g - 1`: how many e-folds per disk scale length the recent share falls (the rate's steepening over the
    stars'), 0.4 for `n = 1.4` and a gas disk as long as the stellar one."""
    return tuning.STAR_FORMATION_LAW_EXPONENT / tuning.GAS_TO_STELLAR_SCALE_LENGTH_RATIO - 1.0


@functools.lru_cache(maxsize=64)
def _normalisation(decline, gap_over_length, floor, shape):
    """The factor `A` that makes the disk-wide mean of `exp(-decline x) * S * G` 1, `x` in disk scale lengths, taking
    each radius as many stars as the disk has there (`x exp(-x)`); `G` is the thin-gas factor of `shape` (whose
    `k_norm` the caller has set to 1, as the history does not read it)."""
    top = 80.0 / (1.0 + decline) if decline > -1.0 else 80.0
    x = np.linspace(0.0, max(top, 40.0), _GRID_POINTS)
    weight = x * np.exp(-x)
    shape_term = np.exp(-decline * x) * _gap_suppression(x, gap_over_length, floor) * _thin_gas_factor(x, shape)[0]
    return float(np.trapezoid(weight, x) / np.trapezoid(weight * shape_term, x))


def recent_factor(radius_pc, shape):
    """
    The recent-star-formation factor at cylindrical radius `radius_pc`: the share of the stars here that formed
    recently, over the disk's average. Above 1 where the disk forms stars faster than its past average, below 1 in
    the outer disk and the bar gap.

    Args:
        radius_pc (float or numpy.ndarray): The distance from the galaxy's axis.
        shape (GalaxyShape): The galaxy's shape (its disk scale length and bulge scale radius set the scales).

    Returns:
        float or numpy.ndarray: Matches `radius_pc`.
    """
    length = shape.disk_scale_length_pc
    gap = gap_radius_pc(shape)
    floor = _gap_floor(shape)
    decline = _decline_per_scale_length()
    flat = shape._replace(k_norm=1.0)
    norm = _normalisation(round(decline, 9), round(gap / length, 9), floor, flat)
    radius = np.asarray(radius_pc, dtype=float)
    value = (norm * np.exp(-decline * radius / length) * _gap_suppression(radius, gap, floor)
             * _thin_gas_factor(radius / length, flat)[0])
    return float(value) if value.ndim == 0 else value


def star_formation_rate_profile(radius_pc, shape):
    """The present star-formation rate per area at `radius_pc`, in units of the disk's lifetime-average rate at the
    centre of a pure exponential: the recent factor times the stars' own surface density (`exp(-R / L)`). The
    quantity the plausibility test compares with the observed profile."""
    radius = np.asarray(radius_pc, dtype=float)
    return recent_factor(radius, shape) * np.exp(-radius / shape.disk_scale_length_pc)


def _log_slope(radius_pc, shape):
    """`d ln(recent_factor) / dR` at `radius_pc` (closed form, arrays)."""
    length = shape.disk_scale_length_pc
    gap = gap_radius_pc(shape)
    floor = _gap_floor(shape)
    power = tuning.BAR_GAP_EDGE_POWER
    radius = np.asarray(radius_pc, dtype=float)
    up = np.power(radius / gap, power)
    rising = (1.0 - floor) * power * up / (1.0 + up) ** 2 / radius
    thin_slope = _thin_gas_factor(radius / length, shape._replace(k_norm=1.0))[1]
    return rising / _gap_suppression(radius, gap, floor) - _decline_per_scale_length() / length + thin_slope / length


@functools.lru_cache(maxsize=64)
def _turning_points(shape):
    """The radii at which the factor stops rising or falling (it falls a little from the centre, rises to the gap's
    edge and falls beyond: a local minimum and a local maximum), found on a log grid and bisected."""
    length = shape.disk_scale_length_pc
    grid = np.geomspace(1e-4 * length, 1e4 * length, 4000)
    slope = _log_slope(grid, shape)
    found = []
    for i in np.nonzero(np.sign(slope[:-1]) != np.sign(slope[1:]))[0]:
        lo, hi = float(grid[i]), float(grid[i + 1])
        rising = slope[i] > 0.0
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if (_log_slope(mid, shape) > 0.0) == rising:
                lo = mid
            else:
                hi = mid
        found.append(0.5 * (lo + hi))
    return tuple(found)


def recent_factor_range(near_pc, far_pc, shape):
    """
    The smallest and largest `recent_factor` at any radius from `near_pc` to `far_pc` (arrays, or numbers): the
    factor is smooth with at most two turning points, so the extremes are at an end of the interval or at a
    turning point inside it. The band bounds of the object-first scatter widen the young and intermediate weights
    by this range (`generation.object_first`).

    Returns:
        tuple: `(smallest, largest)`, matching the shape of `near_pc`.
    """
    near = np.asarray(near_pc, dtype=float)
    far = np.asarray(far_pc, dtype=float)
    smallest = np.minimum(recent_factor(near, shape), recent_factor(far, shape))
    largest = np.maximum(recent_factor(near, shape), recent_factor(far, shape))
    for point in _turning_points(shape):
        inside = (near <= point) & (point <= far)
        at = recent_factor(point, shape)
        smallest = np.where(inside, np.minimum(smallest, at), smallest)
        largest = np.where(inside, np.maximum(largest, at), largest)
    return smallest, largest
