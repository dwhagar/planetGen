# stellarObjects/galaxyGeometry.py

"""
Galaxy-Scale Cylindrical Sector Grid
=======================================

Every galaxy-placed sector is one cell of a cylindrical grid centered on
the galactic origin (see `docs/design/galaxy-coordinate-system.md`,
"Cylindrical sector grid"):

- **Ring** `i` (radial): cylindrical radius `R` in `[i*w, (i+1)*w)`.
- **Layer** `j` (height, signed): `z` in `[(j - 1/2)*h, (j + 1/2)*h)`, so
  layer 0 is centered on the galactic plane.
- **Slot** `k` (angle): ring `i` is cut into `ring_sector_count(i)` equal
  wedges, slot `k` spanning `theta` in `[2*pi*k/N, 2*pi*(k+1)/N)` from
  `+X`, counterclockwise.

`w` and `h` are both the sector edge length (`edge_pc`, a whole number of
parsecs: `program_constants.DEFAULT_SECTOR_EDGE_PC`, 4 pc, by default),
and `N` is chosen so the arc along a ring's centerline is also about one
edge, so every cell is close to an `edge_pc` cube (its arc within about
6% of an edge; `N` is a multiple of the ring's master wedge count so slot
boundaries line up from the center out -- see `ring_sector_count`).

Ring and slot counts are the same on every layer, so the grid is a stack
of identical circular slices: cell `(i, j, k)` sits directly above
`(i, j - 1, k)` and every slot boundary is one vertical plane through the
whole stack ("aligned prism" columns). A layer only holds the rings out to
where the galaxy's density gives out (`galaxySkeleton.build_layer_extents`),
so the slices shrink away from the plane, but any ring a layer has is cut
exactly like that ring on every other layer.

The cells tile space exactly -- no gaps, no overlaps -- and every lookup
is closed-form: `sector_address_at` maps a point straight to its cell,
`neighbor_addresses` lists a cell's face neighbors, and
`enumerate_sectors_within_radius` walks only the cells near a point.

Every distance here is in parsecs, matching `sectors.center_x/y/z_pc`.
Pure geometry: no database, no I/O.
"""

import math

DESIGNATION_SLOT_BITS = 20
"""int: Low bits of a packed designation holding the slot -- enough for
the ~411,000 slots of ring 65,535, far past any galaxy this project
builds (the default Milky Way reaches ring ~3,900, ~24,000 slots)."""

DESIGNATION_LAYER_BITS = 13
"""int: Bits above the slot holding the layer, biased by
`DESIGNATION_LAYER_BIAS` so negative layers pack as plain integers."""

DESIGNATION_LAYER_BIAS = 1 << (DESIGNATION_LAYER_BITS - 1)
"""int: Added to a layer index before packing (layers -4096..4095)."""


MASTER_WEDGE_BASE_COUNT = 3
"""int: Master wedges at the center (ring 0's own slot count)."""

MASTER_WEDGE_MIN_SLOTS = 8
"""int: A new master zone (twice the master wedges) starts at the first
ring where each doubled master wedge would still hold at least this many
slots."""


def _ring_centerline_edges(ring_index):
    """`2*pi*(i + 1/2)`, ring `i`'s centerline circumference in edges.

    Raises:
        ValueError: For a negative ring.
        OverflowError: For a ring so far out the float overflows (as
            `round` of it would), rather than looping forever.
    """
    if ring_index < 0:
        raise ValueError(f"ring_index must be >= 0, got {ring_index}")
    c = 2 * math.pi * (ring_index + 0.5)
    if not math.isfinite(c):
        raise OverflowError(f"ring_index {ring_index} is too far out")
    return c


def ring_master_count(ring_index):
    """
    How many master wedges ring `ring_index` sits in: 3 at the center,
    doubling (6, 12, ... 1,536 at the default galaxy's edge) at the first
    ring where each doubled wedge would hold at least
    `MASTER_WEDGE_MIN_SLOTS` slots. Every master line, once started, runs
    out to the edge, and `ring_sector_count` is always a multiple of this,
    so master lines are slot boundaries in every ring outward -- the
    Galaxy Map's mega-blocks cut their wedges on them.

    Args:
        ring_index (int): `i >= 0`.

    Returns:
        int: `3 * 2**z`.
    """
    c = _ring_centerline_edges(ring_index)
    master = MASTER_WEDGE_BASE_COUNT
    # 2*pi*(i + 1/2) is never a whole number, so the comparison is exact.
    while c >= 2 * master * MASTER_WEDGE_MIN_SLOTS:
        master *= 2
    return master


def ring_sector_count(ring_index):
    """
    How many slots ring `ring_index` holds: the multiple of its master
    wedge count (`ring_master_count`) nearest `2*pi*(i + 1/2)` (the ring's
    centerline circumference in edge lengths) -- 3, 9, 15, 21, 27, 36, 42,
    48, ... (the hybrid master-wedge rule Boss chose on 2026-09-30). Each
    slot's centerline arc stays within about 0.94-1.06 of one edge, the
    total count is within 0.1% of plain rounding, and slot boundaries line
    up on the master lines from the center out. Independent of `edge_pc`,
    which cancels out, and of the layer, so columns line up through every
    layer. `galaxyprisms.js`'s `ringSectorCount` mirrors it.

    Args:
        ring_index (int): `i >= 0`.

    Returns:
        int: `N_i >= 3`.
    """
    c = _ring_centerline_edges(ring_index)
    master = ring_master_count(ring_index)
    return max(master, master * round(c / master))


def ring_radius_pc(ring_index, edge_pc):
    """The cylindrical radius of ring `ring_index`'s centerline, where its
    sector centers sit: `(i + 1/2) * edge_pc`."""
    return (ring_index + 0.5) * edge_pc


def ring_bounds_pc(ring_index, edge_pc):
    """`(inner, outer)` cylindrical radius of ring `ring_index`, parsecs."""
    return ring_index * edge_pc, (ring_index + 1) * edge_pc


def layer_center_z_pc(layer_index, edge_pc):
    """The `z` of layer `layer_index`'s midplane: `j * edge_pc`."""
    return layer_index * edge_pc


def layer_bounds_pc(layer_index, edge_pc):
    """`(bottom, top)` `z` of layer `layer_index`, parsecs."""
    return (layer_index - 0.5) * edge_pc, (layer_index + 0.5) * edge_pc


def slot_angle_bounds(ring_index, slot_index):
    """`(theta_start, theta_end)` of one slot, radians from `+X`."""
    n = ring_sector_count(ring_index)
    step = 2 * math.pi / n
    return slot_index * step, (slot_index + 1) * step


def _check_slot(ring_index, slot_index):
    n = ring_sector_count(ring_index)
    if not (0 <= slot_index < n):
        raise ValueError(
            f"ring_slot_index {slot_index} out of range for ring {ring_index} "
            f"(holds {n} slots, 0..{n - 1})"
        )
    return n


def sector_position_pc(ring_index, layer_index, slot_index, edge_pc):
    """
    The `(x, y, z)` center of one sector, in parsecs: on its ring's
    centerline, at its slot's middle angle, on its layer's midplane.

    Raises:
        ValueError: If `slot_index` is out of range for the ring.
    """
    n = _check_slot(ring_index, slot_index)
    r = ring_radius_pc(ring_index, edge_pc)
    theta = (slot_index + 0.5) * 2 * math.pi / n
    return (r * math.cos(theta), r * math.sin(theta), layer_center_z_pc(layer_index, edge_pc))


def cylindrical_radius_pc(position):
    """`sqrt(x^2 + y^2)` -- distance from the galactic axis."""
    return math.hypot(position[0], position[1])


def galactic_radius_pc(position):
    """
    `sqrt(x^2 + y^2 + z^2)` -- straight-line distance from the galactic
    center, matching `sectors.galactic_radius_pc`.
    """
    x, y, z = position
    return math.sqrt(x * x + y * y + z * z)


def sector_address_at(position, edge_pc):
    """
    The `(ring_index, layer_index, ring_slot_index)` cell containing a
    galaxy-frame point.

    Args:
        position (tuple): `(x, y, z)`, parsecs.
        edge_pc (float): The sector edge length, parsecs.

    Returns:
        tuple: `(ring_index, layer_index, ring_slot_index)`.
    """
    x, y, z = position
    ring_index = int(math.floor(math.hypot(x, y) / edge_pc))
    layer_index = int(math.floor(z / edge_pc + 0.5))
    n = ring_sector_count(ring_index)
    theta = math.atan2(y, x) % (2 * math.pi)
    slot_index = min(n - 1, int(theta * n / (2 * math.pi)))
    return ring_index, layer_index, slot_index


def sector_orientation(center_pc):
    """
    A sector's local axes in the galaxy frame, the frame
    `star_systems.position_x/y/z_mpc` offsets are expressed in: local
    `+X` points radially outward from the galactic axis, local `+Y` points
    toward increasing `theta` (counterclockwise seen from galactic north),
    and local `+Z` is galactic north. The same convention for every
    sector, so a system's octant means the same thing everywhere.

    Args:
        center_pc (tuple): The sector's own `(x, y, z)` center, parsecs.
            Every grid cell's center is off the galactic axis by at least
            half an edge; a point exactly on the axis (only a hand-placed
            sector can be) gets the galaxy's own axes.

    Returns:
        tuple: `(local_x, local_y, local_z)`, each a unit `(x, y, z)`.
    """
    x, y, _z = center_pc
    r = math.hypot(x, y)
    if r < 1e-12:
        return (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)
    cos_t, sin_t = x / r, y / r
    return (cos_t, sin_t, 0.0), (-sin_t, cos_t, 0.0), (0.0, 0.0, 1.0)


def local_to_galaxy_pc(center_pc, offset_pc):
    """Converts a sector-local offset (parsecs, `sector_orientation`'s
    axes) into an absolute galaxy-frame point."""
    ax, ay, az = sector_orientation(center_pc)
    ox, oy, oz = offset_pc
    return tuple(center_pc[i] + ox * ax[i] + oy * ay[i] + oz * az[i] for i in range(3))


def galaxy_to_local_pc(center_pc, point_pc):
    """`local_to_galaxy_pc`'s inverse: an absolute galaxy-frame point as
    an offset (parsecs) from the sector's center along its local axes."""
    axes = sector_orientation(center_pc)
    delta = tuple(point_pc[i] - center_pc[i] for i in range(3))
    return tuple(sum(delta[i] * axis[i] for i in range(3)) for axis in axes)


def sector_cell_vertices_pc(ring_index, layer_index, slot_index, edge_pc):
    """
    The 8 corners of one cell, galaxy frame, parsecs. List index is
    `4*r_bit + 2*z_bit + theta_bit` (each bit 0 for the low bound, 1 for
    the high one), so corners `i` and `i ^ 1`, `i ^ 2`, `i ^ 4` share an
    edge. The two curved faces (inner and outer ring surfaces) bow
    slightly between their corners; at 4 pc the bow is under 0.2 ly
    from ring 10 outward. Ring 0's cells are pie wedges, so their four
    inner corners all sit on the galactic axis.
    """
    _check_slot(ring_index, slot_index)
    r_bounds = ring_bounds_pc(ring_index, edge_pc)
    z_bounds = layer_bounds_pc(layer_index, edge_pc)
    t_bounds = slot_angle_bounds(ring_index, slot_index)
    vertices = []
    for r in r_bounds:
        for z in z_bounds:
            for theta in t_bounds:
                vertices.append((r * math.cos(theta), r * math.sin(theta), z))
    return vertices


class SectorCell:
    """
    One grid cell in its own sector-local frame (`sector_orientation`'s
    axes, origin at the sector center), in whatever length unit it was
    built with -- `SpaceSector` uses light-years. Lets sector generation
    place systems inside the real cylindrical cell instead of a cube.

    Attributes:
        r_inner, r_outer (float): The ring's inner/outer cylindrical radius.
        r_center (float): The ring centerline radius (the local origin sits
            this far from the galactic axis).
        half_angle (float): Half the slot's angular width, radians.
        half_height (float): Half the layer's height.
    """

    def __init__(self, r_inner, r_outer, half_angle, half_height):
        self.r_inner = r_inner
        self.r_outer = r_outer
        self.r_center = (r_inner + r_outer) / 2
        self.half_angle = half_angle
        self.half_height = half_height

    @classmethod
    def for_ring(cls, ring_index, edge):
        """The cell shape of any sector in ring `ring_index`, for sector
        edge length `edge` (any unit)."""
        n = ring_sector_count(ring_index)
        return cls(ring_index * edge, (ring_index + 1) * edge, math.pi / n, edge / 2)

    @property
    def volume(self):
        """The cell's volume, in the unit cubed."""
        return (self.r_outer ** 2 - self.r_inner ** 2) * self.half_angle * 2 * self.half_height

    def contains(self, point):
        """Whether a sector-local point lies inside the cell (boundary
        inclusive, with a rounding allowance of about a trillionth of the
        cell's size, so every point `sample` draws counts as inside even
        when it lands exactly on a face)."""
        lx, ly, lz = point
        tolerance = 1e-12 * max(1.0, self.r_outer)
        if abs(lz) > self.half_height + tolerance:
            return False
        px = self.r_center + lx
        r = math.hypot(px, ly)
        if r < self.r_inner - tolerance or r > self.r_outer + tolerance:
            return False
        if r <= tolerance:
            return True
        return abs(math.atan2(ly, px)) <= self.half_angle + 1e-12

    def sample(self, rng):
        """A uniformly random sector-local point inside the cell."""
        r = math.sqrt(rng.uniform(self.r_inner ** 2, self.r_outer ** 2))
        theta = rng.uniform(-self.half_angle, self.half_angle)
        z = rng.uniform(-self.half_height, self.half_height)
        return (r * math.cos(theta) - self.r_center, r * math.sin(theta), z)


def _overlapping_slots(from_ring, slot_index, to_ring):
    """Slots of `to_ring` whose angular span overlaps slot `slot_index`
    of `from_ring` (touching at a single boundary angle doesn't count).
    Integer arithmetic, so coinciding boundaries -- every master line
    (`ring_master_count`) is one -- are exact."""
    n_from = ring_sector_count(from_ring)
    n_to = ring_sector_count(to_ring)
    first = (slot_index * n_to) // n_from
    last = -((-(slot_index + 1) * n_to) // n_from) - 1  # ceil(...) - 1
    return list(range(first, last + 1))


def neighbor_addresses(ring_index, layer_index, slot_index):
    """
    Every cell sharing a face (of nonzero area) with this one: the two
    slots either side in the same ring and layer, the cells directly
    above and below, and the one or two overlapping slots in each
    adjacent ring (ring slot counts differ, so cells in neighboring rings
    are offset like brick courses). Ring 0 has no inward neighbor.

    Returns:
        list[tuple]: `(ring_index, layer_index, ring_slot_index)` tuples,
            no duplicates, not including the cell itself.
    """
    n = _check_slot(ring_index, slot_index)
    result = [
        (ring_index, layer_index, (slot_index - 1) % n),
        (ring_index, layer_index, (slot_index + 1) % n),
        (ring_index, layer_index - 1, slot_index),
        (ring_index, layer_index + 1, slot_index),
    ]
    for other in (ring_index - 1, ring_index + 1):
        if other < 0:
            continue
        result.extend((other, layer_index, s) for s in _overlapping_slots(ring_index, slot_index, other))
    seen = set()
    unique = []
    for address in result:
        if address not in seen and address != (ring_index, layer_index, slot_index):
            seen.add(address)
            unique.append(address)
    return unique


def sector_zone(ring_index, edge_ly, zone_target_ly=100.0):
    """
    The Zone index -- a fixed-width band of consecutive rings about
    `zone_target_ly` light-years wide that the Galaxy pages group sectors
    by (`html/lib/galaxymap.py`'s identically-named concept). Takes
    light-years, unlike the rest of this module.

    Returns:
        int: `ring_index // round(zone_target_ly / edge_ly)`.
    """
    zone_ring_width = max(1, round(zone_target_ly / edge_ly))
    return ring_index // zone_ring_width


def sector_quadrant(x_pc, y_pc):
    """
    Classifies a galaxy-frame `(x, y)` position into one of 4 azimuthal
    Quadrants, numbered 1-4 counterclockwise from `+X` -- the same split
    as `html/lib/galaxymap.py`'s `sector_quadrant` (labelled "I"-"IV").

    Returns:
        int: 1, 2, 3, or 4.
    """
    theta = math.atan2(y_pc, x_pc) % (2 * math.pi)
    return min(3, int(theta // (math.pi / 2))) + 1


def provisional_sector_designation(ring_index, layer_index, slot_index):
    """
    A short, reversible code for a sector address, for referring to one
    before (or without) generating it: ring, biased layer and slot
    bit-packed into one integer and printed as uppercase hex, e.g.
    `"FE81000A2B"`. `parse_sector_designation` undoes it.

    Raises:
        ValueError: If the ring is negative, the slot is out of range for
            the ring (`parse_sector_designation` would reject the result),
            or the layer or slot doesn't fit its bit field.
    """
    if ring_index < 0:
        raise ValueError(f"ring_index must be >= 0, got {ring_index}")
    _check_slot(ring_index, slot_index)
    biased_layer = layer_index + DESIGNATION_LAYER_BIAS
    if not (0 <= biased_layer < (1 << DESIGNATION_LAYER_BITS)):
        raise ValueError(f"layer_index {layer_index} out of designation range")
    if not (0 <= slot_index < (1 << DESIGNATION_SLOT_BITS)):
        raise ValueError(f"ring_slot_index {slot_index} out of designation range")
    packed = (
        (ring_index << (DESIGNATION_LAYER_BITS + DESIGNATION_SLOT_BITS))
        | (biased_layer << DESIGNATION_SLOT_BITS)
        | slot_index
    )
    return f"{packed:X}"


def parse_sector_designation(designation):
    """
    The `(ring_index, layer_index, ring_slot_index)` a designation from
    `provisional_sector_designation` encodes.

    Raises:
        ValueError: If it isn't valid hex or names an out-of-range slot.
    """
    packed = int(designation.strip(), 16)
    slot_index = packed & ((1 << DESIGNATION_SLOT_BITS) - 1)
    biased_layer = (packed >> DESIGNATION_SLOT_BITS) & ((1 << DESIGNATION_LAYER_BITS) - 1)
    ring_index = packed >> (DESIGNATION_LAYER_BITS + DESIGNATION_SLOT_BITS)
    _check_slot(ring_index, slot_index)
    return ring_index, biased_layer - DESIGNATION_LAYER_BIAS, slot_index


def _slots_near_angle(n, theta_center, half_width):
    """Slot indices (of an `n`-slot ring) whose center angle lies within
    `half_width` of `theta_center`, wrapping around; every slot when the
    window covers the whole ring."""
    step = 2 * math.pi / n
    if half_width * 2 >= 2 * math.pi - step:
        return range(n)
    # Slot k's center is (k + 0.5) * step.
    first = math.ceil((theta_center - half_width) / step - 0.5 - 1e-9)
    last = math.floor((theta_center + half_width) / step - 0.5 + 1e-9)
    return sorted({k % n for k in range(first, last + 1)})


def enumerate_sectors_within_radius(center, radius_pc, edge_pc):
    """
    Every sector address whose center lies within `radius_pc` of `center`
    -- the primitive behind `generate.py galaxy`'s neighborhood modes and
    the Galaxy Map's planned-sector tier. Visits only the rings and
    layers that can reach the sphere and, within each, only the slots in
    the matching angular window, so the cost scales with the answer, not
    with a ring's slot count.

    Args:
        center (tuple): `(x, y, z)`, parsecs -- any point.
        radius_pc (float): `>= 0`.
        edge_pc (float): `> 0`.

    Yields:
        tuple: `(ring_index, layer_index, ring_slot_index, x, y, z,
            distance_pc)` for every qualifying cell, in no particular
            order; `(x, y, z)` is that cell's own center.

    Raises:
        ValueError: If `radius_pc < 0` or `edge_pc <= 0`.
    """
    if radius_pc < 0:
        raise ValueError(f"radius_pc must be >= 0, got {radius_pc}")
    if edge_pc <= 0:
        raise ValueError(f"edge_pc must be > 0, got {edge_pc}")

    cx, cy, cz = center
    r_c = math.hypot(cx, cy)
    theta_c = math.atan2(cy, cx)
    eps = 1e-9

    ring_min = max(0, math.ceil((r_c - radius_pc) / edge_pc - 0.5 - eps))
    ring_max = math.floor((r_c + radius_pc) / edge_pc - 0.5 + eps)
    layer_min = math.ceil((cz - radius_pc) / edge_pc - eps)
    layer_max = math.floor((cz + radius_pc) / edge_pc + eps)

    for ring_index in range(ring_min, ring_max + 1):
        r_m = ring_radius_pc(ring_index, edge_pc)
        n = ring_sector_count(ring_index)
        for layer_index in range(layer_min, layer_max + 1):
            z_m = layer_center_z_pc(layer_index, edge_pc)
            dz = z_m - cz
            planar_sq = radius_pc * radius_pc - dz * dz
            if planar_sq < 0:
                continue
            planar = math.sqrt(planar_sq)
            if abs(r_m - r_c) > planar + eps:
                continue
            if r_c < 1e-9 or r_m + r_c <= planar:
                slots = range(n)
            else:
                cos_max = (r_m * r_m + r_c * r_c - planar_sq) / (2 * r_m * r_c)
                half_width = math.acos(max(-1.0, min(1.0, cos_max)))
                slots = _slots_near_angle(n, theta_c, half_width)
            for slot_index in slots:
                theta = (slot_index + 0.5) * 2 * math.pi / n
                x, y = r_m * math.cos(theta), r_m * math.sin(theta)
                dist = math.dist((x, y, z_m), center)
                if dist <= radius_pc + eps:
                    yield (ring_index, layer_index, slot_index, x, y, z_m, dist)


def describe_sector_cell(ring_index, layer_index, slot_index, edge_pc):
    """
    Everything the grid knows about one cell, for display and the API
    (`GET /api/galaxy/cell`): its address and designation, its center in
    Cartesian `(x, y, z)`, cylindrical `(R, theta, z)` and spherical
    `(r, theta, polar)` coordinates (angles in radians; polar measured
    from galactic north), its bounds, mean arc length, volume and 8
    corners (`sector_cell_vertices_pc`'s order). Parsecs throughout.

    Raises:
        ValueError: If `slot_index` is out of range for the ring.
    """
    _check_slot(ring_index, slot_index)
    center = sector_position_pc(ring_index, layer_index, slot_index, edge_pc)
    r_inner, r_outer = ring_bounds_pc(ring_index, edge_pc)
    z_bottom, z_top = layer_bounds_pc(layer_index, edge_pc)
    theta_start, theta_end = slot_angle_bounds(ring_index, slot_index)
    r_center = ring_radius_pc(ring_index, edge_pc)
    theta = (theta_start + theta_end) / 2
    r_3d = galactic_radius_pc(center)
    cell = SectorCell.for_ring(ring_index, edge_pc)
    return {
        "ring_index": ring_index, "layer_index": layer_index, "ring_slot_index": slot_index,
        "designation": provisional_sector_designation(ring_index, layer_index, slot_index),
        "cartesian_pc": list(center),
        "cylindrical": {"r_pc": r_center, "theta_rad": theta, "z_pc": center[2]},
        "spherical": {
            "r_pc": r_3d, "theta_rad": theta,
            "polar_rad": math.acos(center[2] / r_3d) if r_3d > 0 else 0.0,
        },
        "bounds": {
            "r_inner_pc": r_inner, "r_outer_pc": r_outer,
            "z_bottom_pc": z_bottom, "z_top_pc": z_top,
            "theta_start_rad": theta_start, "theta_end_rad": theta_end,
        },
        "mean_arc_length_pc": r_center * (theta_end - theta_start),
        "volume_pc3": cell.volume,
        "vertices_pc": [list(v) for v in sector_cell_vertices_pc(ring_index, layer_index, slot_index, edge_pc)],
    }
