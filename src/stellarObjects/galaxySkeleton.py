# stellarObjects/galaxySkeleton.py

"""
Galaxy-Wide Density Skeleton: Per-Layer Extents
==================================================

Finds, for a given `galaxyDensity.GalaxyShape`, the one structural fact
worth precomputing about the whole galaxy: which layers of the
cylindrical sector grid (`galaxyGeometry`) hold any content, and how far
out each one reaches. A sector's own position, density and outline are
pure functions of its `(ring_index, layer_index, ring_slot_index)`
address, so none of that is stored -- `generate.ensure_sector_generated`
recomputes it on demand.

**The galaxy as a stack of slices.** The grid's layers run from the
highest layer whose density could clear the qualification threshold (a
sector expected to hold at least one star) to the lowest. Each layer is a
flat circular slice -- rings `0 .. outer_ring_index` around the axis, each
ring cut into the same wedges on every layer -- and it ends at the first
ring where that threshold can no longer be met. Nothing "zero density"
exists in the model (it is exponential), so the threshold is the edge.

**Why one outer ring per layer is exact.** A sector center's
`relative_density` is `bulge(r_3d) + disk(R) * f_z(z) * arm_factor(R,
theta)` (`galaxyDensity._raw_density`). At a fixed height `z` and
cylindrical radius `R` only `arm_factor` varies with angle, and its
maximum is exactly `1 + arm_amplitude`, so `bound_relative_density_at` is
an exact upper bound on any slot in that ring and layer. That bound falls
strictly as `R` grows (both `bulge` and `disk` do) and as `|z|` grows
(`bulge` and `f_z` do), and it is symmetric in `z`. So each layer's
qualifying rings are one run from ring 0 outward, the layers form one
band symmetric about the plane, and each layer reaches no farther out
than the layer below it (toward the plane) -- `build_layer_extents`
walks the whole outline in one pass.

**Why this is a safe superset.** A slot inside a layer's extent may still
fall short (an inter-arm trough, say); that exact per-slot check is one
cheap `relative_density` call, made only when the slot is actually
visited.

Pure and side-effect-free, like `galaxyGeometry.py`/`galaxyDensity.py`.
"""

import bisect
import math
from collections import namedtuple

from planetgen import tuning
from .galaxyDensity import _sech_squared
from .galaxyGeometry import layer_center_z_pc, ring_radius_pc, ring_sector_count
from .spaceSector import SpaceSector

MAX_LAYER_SCAN = 1 << 12
"""int: The farthest layer `build_layer_extents` will walk from the plane
-- matches the designation's layer range, and is 13x the ~320 layers a
Milky-Way-scale galaxy actually needs."""


def expected_system_count_at_density_1(edge_ly=tuning.DEFAULT_SECTOR_EDGE_LY):
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


DEFAULT_MAX_RING = 100000
"""int: Hard cap on how far out a layer is walked, so a pathological
shape with no real outward decay can't scan forever (real Milky-Way-scale
parameters end around ring 3,900)."""


def _qualifies(shape, edge_pc, ring_index, layer_index, threshold_rho):
    return bound_relative_density_at(
        shape, ring_radius_pc(ring_index, edge_pc), layer_center_z_pc(layer_index, edge_pc),
    ) >= threshold_rho


def build_layer_extents(shape, edge_pc, threshold_rho, max_ring=DEFAULT_MAX_RING, max_layer=MAX_LAYER_SCAN):
    """
    The galaxy's outline as a stack of layers: every layer from the
    highest to the lowest that holds any content, each with the last ring
    it reaches.

    Walks the plane's layer outward from ring 0, then climbs one layer at
    a time, pulling the outer ring inward from the layer below (see the
    module docstring for why that is exact), so the cost is about one
    check per ring plus one per layer.

    Args:
        shape (galaxyDensity.GalaxyShape): The galaxy's shape parameters.
        edge_pc (float): The sector edge length (ring width and layer
            height), parsecs.
        threshold_rho (float): The `relative_density` a sector center must
            reach (`1 / expected_system_count_at_density_1(...)`).
        max_ring (int): Stop walking a layer outward past this ring.
        max_layer (int): Stop climbing past this layer.

    Returns:
        tuple: `(extents, outer_ring_index, edge_confirmed)` -- `extents`
            is a list of `(layer_index, outer_ring_index)`, highest layer
            first, symmetric about layer 0 (empty if even the center of
            the plane can't qualify); `outer_ring_index` is the plane's
            (the galaxy's widest) outer ring, `-1` if there are no layers;
            `edge_confirmed` is `False` when `max_ring` cut the plane's
            walk short, so the galaxy's true edge wasn't reached.
    """
    if not _qualifies(shape, edge_pc, 0, 0, threshold_rho):
        return [], -1, True

    outer = 0
    while outer < max_ring and _qualifies(shape, edge_pc, outer + 1, 0, threshold_rho):
        outer += 1
    edge_confirmed = outer < max_ring
    upper = [(0, outer)]

    layer_index = 0
    ring = outer
    while layer_index < max_layer and _qualifies(shape, edge_pc, 0, layer_index + 1, threshold_rho):
        layer_index += 1
        while not _qualifies(shape, edge_pc, ring, layer_index, threshold_rho):
            ring -= 1
        upper.append((layer_index, ring))

    extents = [(j, r) for j, r in reversed(upper)] + [(-j, r) for j, r in upper[1:]]
    return extents, outer, edge_confirmed


def candidate_sector_count(extents):
    """
    Total candidate sectors in a set of layer extents: every slot of
    every ring each layer reaches.

    Args:
        extents (iterable): `(layer_index, outer_ring_index)` pairs.

    Returns:
        int: The count (an upper bound on qualifying sectors).
    """
    cumulative = [0]
    total = 0
    for layer_index, outer_ring_index in extents:
        while len(cumulative) <= outer_ring_index + 1:
            cumulative.append(cumulative[-1] + ring_sector_count(len(cumulative) - 1))
        total += cumulative[outer_ring_index + 1]
    return total


def column_extents(layer_extents):
    """
    The same outline seen from the side: for each ring, the highest and
    lowest layer its column of sectors reaches (its stack bound). A ring
    reaches layer `j` exactly when layer `j`'s outer ring is at least that
    ring, so this is `layer_extents` turned on its side.

    Args:
        layer_extents (iterable): `(layer_index, outer_ring_index)` pairs,
            as `build_layer_extents` returns them.

    Returns:
        list[tuple]: `(ring_index, layer_index_min, layer_index_max)`, ring
            0 outward, one per ring that any layer reaches.
    """
    extents = sorted(layer_extents)
    top = {}
    covered = -1
    for layer_index, outer_ring_index in reversed(extents):
        for ring_index in range(covered + 1, outer_ring_index + 1):
            top[ring_index] = layer_index
        covered = max(covered, outer_ring_index)
    bottom = {}
    covered = -1
    for layer_index, outer_ring_index in extents:
        for ring_index in range(covered + 1, outer_ring_index + 1):
            bottom[ring_index] = layer_index
        covered = max(covered, outer_ring_index)
    return [(ring_index, bottom[ring_index], top[ring_index]) for ring_index in sorted(top)]


class GalaxyBounds:
    """
    The galaxy's stored outline (`galaxy_layer`), for checking an address
    *before* anything is generated there: every `generate.py galaxy` path
    and visit-time generation go through `contains`.

    Attributes:
        edge_pc (float): The sector edge the outline was built at.
        outer_ring (dict): `{layer_index: outer_ring_index}`.
        top_layer_index (int): The highest layer (`-top_layer_index` is the
            lowest); `-1` when no layer holds content.
        outer_ring_index (int): The widest layer's (the plane's) outer
            ring; `-1` when no layer holds content.
    """

    def __init__(self, layer_extents, edge_pc):
        self.edge_pc = edge_pc
        self.outer_ring = dict(layer_extents)
        self.top_layer_index = max(self.outer_ring) if self.outer_ring else -1
        self.outer_ring_index = max(self.outer_ring.values()) if self.outer_ring else -1
        # Sectors in rings 0 .. r-1 of any one layer, for r = 0 .. outer + 1:
        # the per-ring lookup table (slot counts only depend on the ring).
        self._cumulative = [0]
        for ring_index in range(self.outer_ring_index + 1):
            self._cumulative.append(self._cumulative[-1] + ring_sector_count(ring_index))

    def __bool__(self):
        return bool(self.outer_ring)

    def contains(self, ring_index, layer_index):
        """Whether cell `(ring_index, layer_index, any slot)` lies inside the
        galaxy's outline."""
        outer = self.outer_ring.get(layer_index)
        return outer is not None and 0 <= ring_index <= outer

    def describe_miss(self, ring_index, layer_index):
        """A one-line reason `contains` said no, for error messages."""
        outer = self.outer_ring.get(layer_index)
        if outer is None:
            return (f"layer {layer_index} is outside the galaxy (layers run from {self.top_layer_index} "
                    f"down to {-self.top_layer_index})")
        return f"ring {ring_index} is outside the galaxy (layer {layer_index} ends at ring {outer})"

    def cell_count(self):
        """How many sectors (every slot of every ring of every layer) lie
        inside the outline."""
        return sum(self._cumulative[outer + 1] for outer in self.outer_ring.values())

    def random_address(self, rng, max_ring=None):
        """
        A random `(ring, layer, slot)` inside the outline, every sector
        equally likely (so uniform by volume), optionally only from rings
        `0 .. max_ring`.

        Returns:
            tuple or None: `None` when nothing lies inside.
        """
        layers = []
        total = 0
        for layer_index, outer in sorted(self.outer_ring.items()):
            top = outer if max_ring is None else min(outer, max_ring)
            if top >= 0:
                total += self._cumulative[top + 1]
                layers.append((total, layer_index))
        if total == 0:
            return None
        pick = rng.randrange(total)
        index = bisect.bisect_right([end for end, _layer in layers], pick)
        end, layer_index = layers[index]
        start = layers[index - 1][0] if index else 0
        pick -= start
        ring_index = bisect.bisect_right(self._cumulative, pick) - 1
        return ring_index, layer_index, pick - self._cumulative[ring_index]
