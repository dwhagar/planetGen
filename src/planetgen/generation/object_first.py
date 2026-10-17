# planetgen/generation/object_first.py

"""
The object-first scatter (PERF.58): decide how many objects a shell holds,
then where, instead of walking every band of every shell.

For one shell a certified density majorant `M_b` bounds the density anywhere
in band `b` of that shell. The pass draws one Poisson count `N` with mean
`e * sum(cells_b * M_b)` (`cells_b` the band's wedges weighted by their volume in
edge cubes), gives each of the `N` candidates a band in proportion to that
(a cumulative table), a wedge and a point in the cell, and keeps it with
probability `true density / M_b`. That is exactly a Poisson process with the
true density, so the counts and the spatial distribution match a band-by-band
walk; an empty shell costs only its majorant. Study:
`docs/design/scatter-queue-feasibility.md`.

A sector can take several objects (`sector_capacity`): independent draws, a
count per sector, and an object past the sector's capacity dropped. The
capacity is the smallest of `tuning.OBJECT_FIRST_CAPACITIES` for which a
Poisson count with the sector's expected objects exceeds it with a chance
below `tuning.OBJECT_FIRST_TAIL`.
"""

import itertools
import math

import numpy as np

from planetgen import tuning
from planetgen.galaxy import density, spherical
from planetgen.galaxy import star_formation as starFormation


def _sech_squared(x):
    """`density._sech_squared` on an array (`1 / cosh(x) ** 2`, never overflowing)."""
    e = np.exp(-2.0 * np.minimum(np.abs(x), 350.0))
    return np.minimum(1.0, 4.0 * e / (1.0 + e) ** 2)


def _vertical(z, scale_height_pc):
    """`density._vertical` on an array."""
    return _sech_squared(z / (2.0 * scale_height_pc))


def _disk_share_bound(shape, z_near, z_far, fractions, rho_near, rho_far):
    """
    An upper bound on `sum(share_p * fraction_p)` over the three thin-disk
    populations in a band's cells, wherever in them (`z_near` and `z_far`
    are arrays: the nearest and farthest `|z|` of each band, `rho_near` and `rho_far` the nearest and farthest
    cylindrical radius). A population's weight is its
    share of disk star formation, times its vertical profile (between the
    band's near and far edge: the profile falls with `|z|`), its arm
    contrast (between `1 - A` and `1 + A`) and, for the young and intermediate, the recent star-formation factor
    (GEN.134), whose true smallest and largest values over the radius range are taken
    (`star_formation.recent_factor_range`: it is not monotone, so the ends are not enough). The sum is a linear
    fraction of the weights, so over the box of their ranges it peaks at a corner.
    """
    recent_low, recent_high = starFormation.recent_factor_range(rho_near, rho_far, shape)
    ages = tuning.STELLAR_POPULATION_AGE_RANGES_GY
    span = tuning.STAR_FORMATION_AGE_RANGE_GY[1] - tuning.STAR_FORMATION_AGE_RANGE_GY[0]
    ranges, f = [], []
    for name, ratio in tuning.STELLAR_POPULATION_SCALE_HEIGHT_RATIO.items():
        height = shape.disk_scale_height_pc * ratio
        share = (ages[name][1] - ages[name][0]) / span
        arm = abs(tuning.STELLAR_POPULATION_ARM_AMPLITUDE[name])
        low = share * _vertical(z_far, height) / height * max(1.0 - arm, 0.0)
        high = share * _vertical(z_near, height) / height * (1.0 + arm)
        if name in starFormation.RECENT_POPULATIONS:
            low, high = low * recent_low, high * recent_high
        ranges.append((low, high))
        f.append(fractions[name])
    # Far off the plane every profile underflows and the old disk takes it all (`population_densities`).
    best = np.full(np.shape(z_near), float(fractions["old"]))
    for corner in itertools.product(*ranges):
        total = sum(corner)
        mix = sum(weight * fraction for weight, fraction in zip(corner, f))
        with np.errstate(divide="ignore", invalid="ignore"):
            best = np.maximum(best, np.where(total > 0.0, mix / np.where(total > 0.0, total, 1.0), 0.0))
    return best


def band_geometry(shell, high_band, edge_pc):
    """
    For the bands `-high_band .. high_band` of a shell (index `band + high_band`), the extremes of
    each band's cells that bound a density falling with the distance from the axis and with `|z|`:
    `(rho_near, z_near, z_far)` arrays, the smallest cylindrical radius, the smallest `|z|` and the
    largest `|z|` of any point in the band.
    """
    r_in, r_out = spherical.shell_bounds_pc(shell, edge_pc)
    north, south = [], []
    for band in range(-high_band, high_band + 1):
        n, s = spherical.band_colatitude_bounds(shell, band)
        north.append(n)
        south.append(s)
    north, south = np.array(north), np.array(south)
    cos_n, cos_s = np.cos(north), np.cos(south)
    rho_near = r_in * np.minimum(np.sin(north), np.sin(south))
    spans = (cos_n > 0.0) & (cos_s < 0.0)
    z_near = np.where(spans, 0.0, np.where(cos_s >= 0.0, r_in * cos_s, r_in * np.abs(cos_n)))
    z_far = np.where(spans, r_out * np.maximum(cos_n, -cos_s),
                     np.where(cos_s >= 0.0, r_out * cos_n, r_out * np.abs(cos_s)))
    return rho_near, z_near, z_far


def band_majorants(shape, shell, high_band, edge_pc, fractions):
    """
    Per band `-high_band .. high_band` of a shell (index `band + high_band`), an upper bound on
    `sum_p density_p * fractions[p]` (the stars of the pass expected per unit of `e`) anywhere in the
    band's cells. Every component falls with the distance from the axis and with `|z|`, so the band's
    smallest cylindrical radius and its nearest `|z|` bound it; the arm term is bounded by `1 + |A|`,
    the populations' split by `_disk_share_bound`, the halo by its value at the band's nearest point (it falls
    with the spheroidal radius, which both the nearest radius and the nearest height lower).

    Returns:
        numpy.ndarray: One bound per band.
    """
    terms = density.model_terms(shape)
    inner, z_near, z_far = band_geometry(shell, high_band, edge_pc)
    thin = np.exp(-inner / shape.disk_scale_length_pc) * _vertical(z_near, shape.disk_scale_height_pc)
    thick = (terms["thick_disk_amplitude"] * np.exp(-inner / terms["thick_disk_scale_length_pc"])
             * _vertical(z_near, terms["thick_disk_scale_height_pc"]))
    along = inner / shape.bulge_scale_radius_pc
    up = z_near / terms["bulge_scale_z_pc"]
    bulge = shape.bulge_amplitude * np.exp(-0.5 * np.sqrt(along ** 4 + up ** 4))
    arm = shape.arm_level * (1.0 + abs(shape.arm_amplitude))
    outer = spherical.shell_bounds_pc(shell, edge_pc)[1]
    disk_share = _disk_share_bound(shape, z_near, z_far, fractions, inner, outer)
    k = shape.k_norm
    halo = density._halo(inner, z_near, density.halo_terms(shape))
    # Both nuclear terms fall with the radius or the height alone, so the band's nearest point bounds them
    # (`sqrt(inner^2 + z_near^2)` is no farther from the centre than any point of the band).
    nuclear_cluster = density._nuclear_cluster_array(np.sqrt(inner ** 2 + z_near ** 2), terms)
    nuclear_disc = density._nuclear_disc_array(inner, z_near, terms)
    bound = (np.maximum(k * bulge, 0.0) * fractions["bulge"]
             + np.maximum(k * nuclear_cluster, 0.0) * fractions["nuclear_cluster"]
             + np.maximum(k * nuclear_disc, 0.0) * fractions["nuclear_disc"]
             + np.maximum(k * thick, 0.0) * fractions["old"]
             + np.maximum(k * halo, 0.0) * fractions["halo"]
             + np.maximum(k * thin, 0.0) * arm * disk_share)
    return bound * (1.0 + 1e-12)


def sector_capacity(expected_objects):
    """
    How many objects a sector may take: the smallest of
    `tuning.OBJECT_FIRST_CAPACITIES` that a Poisson count with mean
    `expected_objects` exceeds with a chance below `tuning.OBJECT_FIRST_TAIL`
    (the largest, if none does).
    """
    capacities = tuning.OBJECT_FIRST_CAPACITIES
    if expected_objects <= 0.0:
        return capacities[0]
    term = math.exp(-expected_objects)
    cumulative = term
    seen = 0
    for capacity in capacities:
        while seen < capacity:
            seen += 1
            term *= expected_objects / seen
            cumulative += term
        if 1.0 - cumulative < tuning.OBJECT_FIRST_TAIL:
            return capacity
    return capacities[-1]

