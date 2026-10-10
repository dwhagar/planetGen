# planetgen/generation/object_first.py

"""
The object-first scatter (PERF.58): decide how many objects a layer holds,
then where, instead of walking every ring of every layer.

For one layer a certified density majorant `M_r` bounds the density anywhere
in ring `r` of that layer. The pass draws one Poisson count `N` with mean
`e * sum(slots_r * M_r)`, gives each of the `N` candidates a ring in
proportion to `slots_r * M_r` (a cumulative table), a slot and a point in the
cell as before, and keeps it with probability `true density / M_r`. That is
exactly a Poisson process with the true density, so the counts and the spatial
distribution match the ring-by-ring walk; an empty layer costs only its
majorant. Study: `docs/design/scatter-queue-feasibility.md`.

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
from planetgen.galaxy import density


def _disk_share_bound(shape, z_near, z_far, fractions):
    """
    An upper bound on `sum(share_p * fraction_p)` over the three thin-disk
    populations in a layer, wherever in it. A population's weight is its
    share of disk star formation, times its vertical profile (between the
    layer's near and far edge: the profile falls with `|z|`) and its arm
    contrast (between `1 - A` and `1 + A`). The sum is a linear fraction of
    the weights, so over the box of their ranges it peaks at a corner.
    """
    ages = tuning.STELLAR_POPULATION_AGE_RANGES_GY
    span = tuning.STAR_FORMATION_AGE_RANGE_GY[1] - tuning.STAR_FORMATION_AGE_RANGE_GY[0]
    names, ranges, f = [], [], []
    for name, ratio in tuning.STELLAR_POPULATION_SCALE_HEIGHT_RATIO.items():
        height = shape.disk_scale_height_pc * ratio
        share = (ages[name][1] - ages[name][0]) / span
        arm = abs(tuning.STELLAR_POPULATION_ARM_AMPLITUDE[name])
        low = share * density._vertical(z_far, height) / height * max(1.0 - arm, 0.0)
        high = share * density._vertical(z_near, height) / height * (1.0 + arm)
        names.append(name)
        ranges.append((low, high))
        f.append(fractions[name])
    # Far off the plane every profile underflows and the old disk takes it all (`population_densities`).
    best = fractions["old"]
    for corner in itertools.product(*ranges):
        total = sum(corner)
        if total > 0.0:
            best = max(best, sum(weight * fraction for weight, fraction in zip(corner, f)) / total)
    return best


def ring_majorants(shape, layer_index, outer_ring, edge_pc, fractions):
    """
    Per ring `0..outer_ring`, an upper bound on `sum_p density_p * fractions[p]`
    (the stars of the pass expected per unit of `e`) anywhere in the ring at
    this layer. Every component falls with radius and `|z|`, so the ring's
    inner edge and the layer's edge nearest the plane bound it; the arm term
    is bounded by `1 + |A|`, the populations' split by `_disk_share_bound`,
    the halo floor by `MIN_RELATIVE_DENSITY`.

    Returns:
        numpy.ndarray: One bound per ring.
    """
    terms = density.model_terms(shape)
    z_low, z_high = (layer_index - 0.5) * edge_pc, (layer_index + 0.5) * edge_pc
    z_near = 0.0 if z_low <= 0.0 <= z_high else min(abs(z_low), abs(z_high))
    z_far = max(abs(z_low), abs(z_high))
    inner = np.arange(outer_ring + 1, dtype=float) * edge_pc
    thin = np.exp(-inner / shape.disk_scale_length_pc) * density._vertical(z_near, shape.disk_scale_height_pc)
    thick = (terms["thick_disk_amplitude"] * np.exp(-inner / terms["thick_disk_scale_length_pc"])
             * density._vertical(z_near, terms["thick_disk_scale_height_pc"]))
    along = inner / shape.bulge_scale_radius_pc
    up = z_near / terms["bulge_scale_z_pc"]
    bulge = shape.bulge_amplitude * np.exp(-0.5 * np.sqrt(along ** 4 + up ** 4))
    if shape.core_amplitude:
        bulge = bulge + shape.core_amplitude * np.exp(-0.5 * (inner ** 2 + z_near ** 2) / terms["core_scale_radius_pc"] ** 2)
    arm = shape.arm_level * (1.0 + abs(shape.arm_amplitude))
    disk_share = _disk_share_bound(shape, z_near, z_far, fractions)
    k = shape.k_norm
    bound = (np.maximum(k * bulge, 0.0) * fractions["bulge"]
             + (np.maximum(k * thick, 0.0) + tuning.MIN_RELATIVE_DENSITY) * fractions["old"]
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

