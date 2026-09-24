# stellarObjects/galaxySkeleton.py

"""
Galaxy-Wide Density Skeleton: Per-Ring Layer Bands
=====================================================

Finds, for a given `galaxyDensity.GalaxyShape`, the one structural fact
worth precomputing about the whole galaxy: for each ring of the
cylindrical sector grid (`galaxyGeometry`), which layers could hold any
content at all. A sector's own position, density and outline are pure
functions of its `(ring_index, layer_index, ring_slot_index)` address,
so none of that is stored -- `generate.ensure_sector_generated`
recomputes it on demand.

**Why a ring's qualifying region is one band of layers.** A sector
center's `relative_density` is `bulge(r_3d) + disk(R) * f_z(z) *
arm_factor(R, theta)` (`galaxyDensity._raw_density`). Every center in
one ring sits at the same cylindrical radius `R`, so at a fixed layer
height `z` the only thing that varies around the ring is `arm_factor`,
whose maximum over `theta` is exactly `1 + arm_amplitude`. That makes
`bound_relative_density_at` an *exact* upper bound on what any slot in
that ring and layer could reach. It is symmetric in `z` and falls
strictly as `|z|` grows (both `bulge` and `f_z` do), so the layers that
clear the threshold always form one contiguous band centered on the
plane -- `find_ring_band` walks outward from layer 0 until the bound
drops below it.

**Why this is a safe superset.** A slot inside the band may still fall
short (an inter-arm trough, say); that exact per-slot check is one cheap
`relative_density` call, made only when the slot is actually visited.

Pure and side-effect-free, like `galaxyGeometry.py`/`galaxyDensity.py`.
"""

import math
from collections import namedtuple

from . import program_constants
from .galaxyDensity import _sech_squared
from .galaxyGeometry import ring_radius_pc
from .spaceSector import SpaceSector

RingBand = namedtuple("RingBand", ["layer_index_min", "layer_index_max"])
"""The inclusive range of layers in one ring that could hold content."""

MAX_LAYER_SCAN = 1 << 12
"""int: The farthest layer `find_ring_band` will walk from the plane --
matches the designation's layer range, and is 20x the ~340 layers a
Milky-Way-scale galaxy actually needs."""


def expected_system_count_at_density_1(edge_ly=program_constants.DEFAULT_SECTOR_EDGE_LY):
    """
    The `E` constant every qualification check compares against: how many
    systems a sector this size would hold at `relative_density = 1.0`
    (real local stellar density).

    Args:
        edge_ly (float): The sector edge length, in light-years.

    Returns:
        float: Expected system count at density 1.
    """
    return SpaceSector(name="galaxySkeleton-calibration", edge_ly=edge_ly).expected_system_count()


def _bound_raw_density_at(shape, r_cyl, z):
    """The exact maximum over `theta` of `galaxyDensity._raw_density` at
    cylindrical radius `r_cyl` and height `z` (unnormalized)."""
    r_3d = math.hypot(r_cyl, z)
    bulge = shape.bulge_amplitude * math.exp(-r_3d / shape.bulge_scale_radius_pc)
    disk_radial = math.exp(-r_cyl / shape.disk_scale_length_pc)
    f_z = _sech_squared(z / shape.disk_scale_height_pc)
    return bulge + disk_radial * f_z * (1.0 + shape.arm_amplitude)


def bound_relative_density_at(shape, r_cyl, z):
    """
    `_bound_raw_density_at` normalized by `shape.k_norm` -- the highest
    `relative_density` any point at `(r_cyl, z)` can have, whatever its
    angle. Directly comparable to a qualification threshold.
    """
    return shape.k_norm * _bound_raw_density_at(shape, r_cyl, z)


def find_ring_band(shape, edge_pc, ring_index, threshold_rho, max_layer=MAX_LAYER_SCAN):
    """
    The layers of ring `ring_index` whose sector centers could clear
    `threshold_rho`, or `None` if even the plane can't.

    Args:
        shape (galaxyDensity.GalaxyShape): The galaxy's shape parameters.
        edge_pc (float): The sector edge length (ring width and layer
            height), parsecs.
        ring_index (int): The ring to scan.
        threshold_rho (float): The `relative_density` a sector center must
            reach (`1 / expected_system_count_at_density_1(...)`).
        max_layer (int): Stop walking outward past this layer.

    Returns:
        RingBand or None: Symmetric about the plane
            (`layer_index_min == -layer_index_max`).
    """
    r_cyl = ring_radius_pc(ring_index, edge_pc)
    if bound_relative_density_at(shape, r_cyl, 0.0) < threshold_rho:
        return None
    top = 0
    while top < max_layer and bound_relative_density_at(shape, r_cyl, (top + 1) * edge_pc) >= threshold_rho:
        top += 1
    return RingBand(-top, top)


DEFAULT_EMPTY_STREAK_TO_STOP = 50
"""int: How many consecutive empty rings confirm the galaxy's edge has
been reached -- the qualifying region shrinks outward, and the margin
guards against one razor-thin ring a coarse scan could miss."""

DEFAULT_MAX_RING = 100000
"""int: Hard cap on rings scanned, so a pathological shape with no real
outward decay can't scan forever (real Milky-Way-scale parameters end
around ring 4,100)."""


def build_ring_bands(shape, edge_pc, threshold_rho,
                     empty_streak_to_stop=DEFAULT_EMPTY_STREAK_TO_STOP, max_ring=DEFAULT_MAX_RING):
    """
    Scans rings outward from the core until `empty_streak_to_stop`
    consecutive rings hold nothing (or `max_ring` is passed), returning
    every ring's band.

    Returns:
        tuple: `(bands, outer_ring_index, edge_confirmed)` -- `bands` is a
            list of `(ring_index, layer_index_min, layer_index_max)`;
            `outer_ring_index` is the last ring with a band (`-1` if none
            does); `edge_confirmed` is `False` when `max_ring` stopped the
            scan before a full empty streak did.
    """
    bands = []
    outer_ring_index = -1
    streak = 0
    for ring_index in range(max_ring + 1):
        band = find_ring_band(shape, edge_pc, ring_index, threshold_rho)
        if band is None:
            streak += 1
            if streak >= empty_streak_to_stop:
                return bands, outer_ring_index, True
            continue
        streak = 0
        outer_ring_index = ring_index
        bands.append((ring_index, band.layer_index_min, band.layer_index_max))
    return bands, outer_ring_index, False
