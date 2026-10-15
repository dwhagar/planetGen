# planetgen/galaxy/skeleton.py

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

**Why one outer ring per layer is enough.** A sector center's
`relative_density` is `bulge(x, y, z) + thin(R, z) * arm_factor(R, theta)
+ thick(R, z)` (`galaxyDensity._raw_density`). At a fixed height `z` and
cylindrical radius `R`, the bar bulge is largest along its long axis and
`arm_factor` is at most `arm_level * (1 + arm_amplitude)`, so
`bound_relative_density_at` (each term at its own maximum) is an upper
bound on any slot in that ring and layer -- exact wherever the bulge or
the arms are negligible, a hair high in the bar's arm-crossing rings. That
bound falls strictly as `R` grows and as `|z|` grows (every term does),
and it is symmetric in `z`. So each layer's qualifying rings are one run
from ring 0 outward, the layers form one band symmetric about the plane,
and each layer reaches no farther out than the layer below it (toward the
plane) -- `build_layer_extents` walks the whole outline in one pass.

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
from planetgen.galaxy.density import _vertical, bulge_bound, model_terms
from planetgen.galaxy import object_uid, spherical
from planetgen.galaxy.geometry import layer_center_z_pc, ring_radius_pc, ring_sector_count
from planetgen.galaxy.sector import SpaceSector

MAX_LAYER_SCAN = 1 << 12
"""int: The farthest layer `build_layer_extents` will walk from the plane
-- matches the designation's layer range, and is 4x the ~1,020 layers a
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
    """The maximum over `theta` of `galaxyDensity._raw_density` at
    cylindrical radius `r_cyl` and height `z` (unnormalized), or a hair
    above it: each term at its own maximum (the bar along its long axis,
    the thin disk on an arm crest)."""
    terms = model_terms(shape)
    thin = math.exp(-r_cyl / shape.disk_scale_length_pc) * _vertical(z, shape.disk_scale_height_pc)
    thick = (terms["thick_disk_amplitude"] * math.exp(-r_cyl / terms["thick_disk_scale_length_pc"])
             * _vertical(z, terms["thick_disk_scale_height_pc"]))
    return bulge_bound(r_cyl, z, shape) + thin * shape.arm_level * (1.0 + shape.arm_amplitude) + thick


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
    module docstring for why that is safe), so the cost is about one
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
    *before* anything is generated there: every `planetgen galaxy` path
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


# ---------------------------------------------------------------------------
# The outline on spherical shells (`galaxy/spherical.py`)
# ---------------------------------------------------------------------------

DEFAULT_MAX_SHELL = 100000
"""int: Hard cap on how far out the shells are walked, so a pathological
shape with no real outward decay can't scan forever (the default Milky Way
ends at shell 3,764)."""


def _band_qualifies(shape, edge_pc, shell, band, threshold_rho):
    plane, z = spherical.band_center_pc(shell, band, edge_pc)
    return bound_relative_density_at(shape, plane, z) >= threshold_rho


def _highest_qualifying_band(shape, edge_pc, shell, threshold_rho):
    """
    The highest `|band|` of shell `shell` with a cell that can qualify.
    Walking north from the equatorial band the bound only falls, until the
    last step toward the pole, where the thin-disk exponent flattens and the
    density rises again; so the pole cap is checked first (if it qualifies the
    outline reaches the axis and the whole shell is in), and otherwise the
    qualifying bands are one run from the equator, found by bisection.
    Assumes the equatorial band qualifies.
    """
    top = spherical.band_range(shell)[1]
    if _band_qualifies(shape, edge_pc, shell, top, threshold_rho):
        return top
    low, high = 0, top - 1      # `low` qualifies; `high` is the last band that might
    while low < high:
        middle = (low + high + 1) // 2
        if _band_qualifies(shape, edge_pc, shell, middle, threshold_rho):
            low = middle
        else:
            high = middle - 1
    return low


def build_shell_extents(shape, edge_pc, threshold_rho, max_shell=DEFAULT_MAX_SHELL):
    """
    The galaxy's outline on the spherical grid: every shell from the core
    out to the last one whose equatorial band can qualify, each with the
    highest band it reaches. A shell's qualifying bands are one run around
    the plane, so `-high_band .. high_band` (both caps in when `high_band` is
    the shell's top band) holds every cell that can qualify. Like the
    cylinder's outline it is a safe superset: the bound takes the best
    longitude, so a wedge inside may still fall short, and that exact check is
    made when the wedge is visited.

    The cost is a bisection per shell (about 14 bound evaluations), so the
    whole default galaxy is a few tens of thousands of calls.

    Args:
        shape (galaxyDensity.GalaxyShape): The galaxy's shape parameters.
        edge_pc (float): The sector edge length (shell depth), parsecs.
        threshold_rho (float): The `relative_density` a sector center must
            reach (`1 / expected_system_count_at_density_1(...)`).
        max_shell (int): Stop walking outward past this shell.

    Returns:
        tuple: `(extents, outer_shell_index, edge_confirmed)` -- `extents` is
            a list of `(shell_index, high_band)`, the core first (empty if
            even the core can't qualify); `outer_shell_index` the last shell
            (`-1` if none); `edge_confirmed` is `False` when `max_shell` cut
            the walk short.
    """
    if not _band_qualifies(shape, edge_pc, 0, 0, threshold_rho):
        return [], -1, True
    extents = [(0, 0)]
    shell = 0
    while shell < max_shell and _band_qualifies(shape, edge_pc, shell + 1, 0, threshold_rho):
        shell += 1
        extents.append((shell, _highest_qualifying_band(shape, edge_pc, shell, threshold_rho)))
    return extents, shell, shell < max_shell


def shell_candidate_counts(extents):
    """
    The cells each shell of an outline holds, for storing beside it: shell
    `k`'s bands `-high_band .. high_band`, every wedge.

    Args:
        extents (iterable): `(shell_index, high_band)` pairs, as
            `build_shell_extents` returns them.

    Returns:
        list[tuple]: `(shell_index, high_band, cell_count)`, same order.
    """
    return [(shell, high, spherical.cells_in_band_range(shell, high)) for shell, high in extents]


def shell_candidate_sector_count(extents):
    """Total candidate sectors in a set of shell extents: every wedge of
    every band each shell reaches (an upper bound on qualifying sectors)."""
    return sum(count for _shell, _high, count in shell_candidate_counts(extents))


def outline_object_id_layout(extents):
    """
    The object-ID layout (`object_uid.Layout`) an outline needs: the
    default 80 bits when its shells, bands and widest band's wedges fit the
    12/12/16-bit birth-sector field (the default galaxy does), else the field
    widened. The plan step calls this, so a galaxy whose addresses cannot be
    written down is refused before anything is stored.

    Args:
        extents (iterable): `(shell_index, high_band)` pairs, as
            `build_shell_extents` returns them (a count may follow).

    Returns:
        object_uid.Layout: The layout.

    Raises:
        ValueError: If even 128 bits do not hold the outline's addresses.
    """
    rows = list(extents)
    if not rows:
        return object_uid.DEFAULT_LAYOUT
    max_shell = max(row[0] for row in rows)
    max_band = max(row[1] for row in rows)
    max_wedges = max(spherical.wedge_count(row[0], 0) for row in rows)
    return object_uid.fit_layout(max_shell, -max_band, max_band, max_wedges)


class ShellBounds:
    """
    The galaxy's stored outline on spherical shells, for checking an address
    *before* anything is generated there: every `planetgen galaxy` path and
    visit-time generation go through `contains`.

    Attributes:
        edge_pc (float): The sector edge the outline was built at.
        high_band (dict): `{shell_index: highest band reached}`; the shell
            holds bands `-high .. high`.
        outer_shell_index (int): The outermost shell; `-1` when nothing holds
            content.
    """

    def __init__(self, shell_extents, edge_pc):
        """`shell_extents` is `(shell_index, high_band)` pairs, or
        `(shell_index, high_band, cell_count)` triples when the count is
        stored (it saves counting every band of every shell again)."""
        self.edge_pc = edge_pc
        self.high_band = {}
        self._cells = {}
        for row in shell_extents:
            self.high_band[row[0]] = row[1]
            self._cells[row[0]] = row[2] if len(row) > 2 else spherical.cells_in_band_range(row[0], row[1])
        self.outer_shell_index = max(self.high_band) if self.high_band else -1
        # Cells in shells 0 .. k-1, for k = 0 .. outer + 1.
        self._cumulative = [0]
        for shell in range(self.outer_shell_index + 1):
            self._cumulative.append(self._cumulative[-1] + self._cells.get(shell, 0))

    def __bool__(self):
        return bool(self.high_band)

    def contains(self, shell, band):
        """Whether cell `(shell, band, any wedge)` lies inside the galaxy's
        outline."""
        high = self.high_band.get(shell)
        return high is not None and abs(band) <= high

    def describe_miss(self, shell, band):
        """A one-line reason `contains` said no, for error messages."""
        high = self.high_band.get(shell)
        if high is None:
            return f"shell {shell} is outside the galaxy (shells run from 0 to {self.outer_shell_index})"
        return f"band {band} is outside the galaxy (shell {shell} reaches bands -{high} to {high})"

    def cell_count(self):
        """How many sectors (every wedge of every band of every shell) lie
        inside the outline."""
        return self._cumulative[-1]

    def random_address(self, rng, max_shell=None):
        """
        A random `(shell, band, wedge)` inside the outline, every sector
        equally likely (so uniform by volume), optionally only from shells
        `0 .. max_shell`.

        Returns:
            tuple or None: `None` when nothing lies inside.
        """
        top = self.outer_shell_index if max_shell is None else min(self.outer_shell_index, max_shell)
        if top < 0:
            return None
        total = self._cumulative[top + 1]
        pick = rng.randrange(total)
        shell = bisect.bisect_right(self._cumulative, pick) - 1
        return spherical.nth_cell_in_band_range(shell, self.high_band[shell], pick - self._cumulative[shell])
