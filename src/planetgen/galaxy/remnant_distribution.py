# planetgen/galaxy/remnant_distribution.py

"""
Where neutron stars and black holes are, against where stars are (GEN.132).

The generator places remnants in proportion to the stellar density, which is
right on average but not by region. The research in
`docs/design/systems/compact-remnant-regions.md` finds, against stars:

- Neutron stars sit in a layer about three times as thick (birth kicks), black
  holes about two and a half times (smaller kicks), so each is rarer in the
  plane and commoner above it than the stars.
- Black holes are a larger share of the remnants toward the core.
- Active radio pulsars follow the [R2750] radial profile, which
  peaks near 3 to 4 kpc and is far thinner than the stars in the inner galaxy.

Each factor multiplies the density a remnant is placed by, so it moves where
remnants are without changing the rate per star at the reference radius and
height. Every size here is an estimate from those trends, not a measurement:
the sources give shapes, not percentages.
"""

import math

import numpy as np

from planetgen import tuning
from planetgen.galaxy import scales
from planetgen.util import detmath


def _exponential_profile(z_pc, scale_pc):
    """An exponential layer's density at height `z_pc`, normalized to a unit column (the stars' own vertical
    profile, `galaxy.density._vertical`, GEN.544)."""
    x = abs(z_pc) / scale_pc
    if x > 700.0:
        return 0.0
    return math.exp(-x) / (2.0 * scale_pc)


def vertical_factor(kind, z_pc, thin_height_pc):
    """
    The kind's density against the stars' at height `z_pc`: the same column
    spread over `tuning.REMNANT_SCALE_HEIGHT_RATIO[kind]` times the thin disk's
    height, so it is below 1 in the plane and above it beyond about one
    scale height. Capped at `tuning.REMNANT_VERTICAL_FACTOR_MAX`; a kind with
    no ratio gets 1.
    """
    ratio = tuning.REMNANT_SCALE_HEIGHT_RATIO.get(kind)
    if ratio is None or thin_height_pc <= 0.0:
        return 1.0
    stars = _exponential_profile(z_pc, thin_height_pc)
    if stars <= 0.0:
        return tuning.REMNANT_VERTICAL_FACTOR_MAX
    return min(_exponential_profile(z_pc, thin_height_pc * ratio) / stars, tuning.REMNANT_VERTICAL_FACTOR_MAX)


def vertical_factor_array(kind, z_pc, thin_height_pc):
    """`vertical_factor` for an array of heights (one per band of a shell, the phenomena scatter's shell at once)."""
    z_pc = np.asarray(z_pc, dtype=float)
    ratio = tuning.REMNANT_SCALE_HEIGHT_RATIO.get(kind)
    if ratio is None or thin_height_pc <= 0.0:
        return np.ones_like(z_pc)

    def profile(scale_pc):
        x = np.abs(z_pc) / scale_pc
        return np.where(x > 700.0, 0.0, np.exp(-np.minimum(x, 700.0)) / (2.0 * scale_pc))

    stars = profile(thin_height_pc)
    thick = profile(thin_height_pc * ratio)
    cap = tuning.REMNANT_VERTICAL_FACTOR_MAX
    return np.where(stars <= 0.0, cap, np.minimum(thick / np.where(stars <= 0.0, 1.0, stars), cap))


def radial_factor(kind, radius_pc):
    """The kind's density against the stars' at galactocentric `radius_pc`:
    black holes gain `tuning.BLACK_HOLE_CORE_EXCESS` extra at the core, falling
    off over `tuning.BLACK_HOLE_CORE_SCALE_PC`; others 1."""
    if kind != "black-hole":
        return 1.0
    return 1.0 + tuning.BLACK_HOLE_CORE_EXCESS * math.exp(-radius_pc / tuning.BLACK_HOLE_CORE_SCALE_PC)


def placement_factor(kind, point_pc, thin_height_pc):
    """Both factors for a galaxy-frame point `(x, y, z)` in parsecs."""
    radius = detmath.hypot(point_pc[0], point_pc[1])
    return vertical_factor(kind, point_pc[2], thin_height_pc) * radial_factor(kind, radius)


def radial_factor_array(kind, radius_pc):
    """`radial_factor` for an array of radii (the phenomena scatter's whole layer at once, PERF.63)."""
    radius_pc = np.asarray(radius_pc, dtype=float)
    if kind != "black-hole":
        return np.ones_like(radius_pc)
    return 1.0 + tuning.BLACK_HOLE_CORE_EXCESS * np.exp(-radius_pc / tuning.BLACK_HOLE_CORE_SCALE_PC)


def pulsar_radial_factor(radius_pc):
    """
    Active radio pulsars against stars at galactocentric `radius_pc`, 1 at
    the reference radius: the [R2750] profile `(R/R0)^a exp(-b (R-R0)/R0)`
    (`tuning.PULSAR_RADIAL_A`, `_B`) over the stellar disk's own exponential, with R0 the galaxy's reference radius
    and the disk's scale length the galaxy's (`galaxy.scales`, GEN.404). For the Milky Way's scales: about 0.1 at
    1 kpc, 0.4 at 3 kpc, 0.8 near 15 kpc. Capped at 2.
    """
    reference = scales.reference_radius_pc()
    if radius_pc <= 0.0:
        return 0.0
    # The ratio is taken in logarithms: far out the disk's exponential underflows to zero, and the old quotient divided by it.
    log_ratio = (tuning.PULSAR_RADIAL_A * math.log(radius_pc / reference)
                 - tuning.PULSAR_RADIAL_B * (radius_pc - reference) / reference
                 + (radius_pc - reference) / scales.disk_scale_length_pc())
    return min(math.exp(min(log_ratio, math.log(2.0))), 2.0)
