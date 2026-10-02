# stellarObjects/galaxyViewport.py

"""
Interactive 3D Galaxy Map viewport queries.

The flat, server-rendered-once SVG Galaxy Map (`html/lib/galaxymap.py`)
plots every galaxy-placed sector in one shot, since the whole point of
that map is a single bird's-eye overview. The interactive 3D map
(`html/lib/galaxymap3d.py`/`static/galaxymap3d.js`) is the opposite: a
real camera that moves freely through the galaxy, so what it needs to
plot changes every time it moves -- this module is the data layer behind
that, called fresh (via `queryDb.galaxy_tiles`) each time the camera's
viewport changes.

Two tiers of content, matching the sprite kinds
`static/galaxymap3d.js` draws:

- **Placed**: real, already-generated sectors -- `queryDb.
  galaxy_sectors_in_view` handles this tier directly (it needs the
  database; nothing here does).
- **Planned**: real, not-yet-generated sector *addresses* -- exact
  `(ring_index, layer_index, ring_slot_index)` grid cells this galaxy's own shape model
  predicts would qualify (see `galaxyDensity.predicted_star_count`, the
  same >= 1 threshold `generate.py`'s `ensure_sector_generated`/
  `_BatchDensity` already gate real generation on), enumerated exactly
  via `galaxyGeometry.enumerate_sectors_within_radius` -- `planned_slots_in_view`.
  Deliberately capped to a small view radius (`PLANNED_RADIUS_CAP_PC`):
  see that function's own docstring for why an unbounded radius here
  would be catastrophic.
Predicted density itself is no longer served from here: the map's prisms
evaluate `galaxyDensity.relative_density` in the browser
(`static/galaxyprisms.js`).

Pure and side-effect-free, like `galaxyGeometry.py`/`galaxyDensity.py` --
no database or I/O; `queryDb.py` owns combining this with the database's
own placed-sector rows.
"""

import heapq
import math

from .galaxyDensity import predicted_star_count, relative_density
from .galaxyGeometry import enumerate_sectors_within_radius, provisional_sector_designation

PLANNED_RADIUS_CAP_PC = 40.0
"""float: The largest view radius `planned_slots_in_view` will actually
enumerate individual addresses for, however large a `radius_pc` its
caller asks for. `enumerate_sectors_within_radius` costs about one step
per address it returns, and a view spanning the whole galaxy holds
billions. Beyond this radius the map shows predicted density only (its
prisms, computed in the browser) -- exact addresses only become worth enumerating once the
view has actually narrowed down to something close to a "handful to a
few dozen real sectors" neighborhood, the scale this primitive was always
designed for (see `galaxyGeometry`'s own module docstring).

Was 200 pc, which with 11.5 ly sectors meant enumerating ~770 thousand
slots (and, before only the closest `PLANNED_MAX_RESULTS` were kept
lazily, building a dict for every one of them) just to return the
closest 4,000 -- ~400 MB and up to 140 s per request, which is what
OOM-killed Apache when several camera moves' requests overlapped. 40 pc
holds ~6,000 slots, already more than `PLANNED_MAX_RESULTS`."""

PLANNED_MAX_RESULTS = 4000
"""int: Hard cap on how many planned-slot entries `planned_slots_in_view`
returns, closest-first -- even within `PLANNED_RADIUS_CAP_PC`, a view
centered deep in the bulge (where nearly every slot qualifies) could
still enumerate more addresses than are worth serializing/rendering in
one response."""


def qualifying_threshold_star_count():
    """
    The predicted-star-count threshold a slot must clear to be treated as
    "planned" (worth showing as a real, addressable, not-yet-generated
    sector) rather than empty space -- exactly `generate.py`'s own
    `ensure_sector_generated`/`_BatchDensity` gate (`predicted_star_count
    >= 1.0`, i.e. "this position is predicted to hold at least one real
    star system"), restated here so this module's own callers don't need
    to import `generate.py` (a script, not a library module) just for one
    constant.

    Returns:
        float: `1.0`.
    """
    return 1.0


def planned_slots_in_view(center_pc, radius_pc, edge_pc, shape,
                           expected_system_count_at_density_1, exclude_addresses,
                           cap=PLANNED_MAX_RESULTS):
    """
    Every real, not-yet-generated sector address within `radius_pc` of
    `center_pc` that this galaxy's own shape model predicts would qualify
    (`predicted_star_count >= 1`) -- the "planned" tier `static/
    galaxymap3d.js` renders as small, dim, clickable dots (see the module
    docstring's own tier breakdown), each carrying its real
    `(ring_index, layer_index, ring_slot_index)` address and a copyable
    provisional designation, ready to feed straight into `generate.py
    galaxy --ring I --layer J --slot K`.

    Args:
        center_pc (tuple): `(x, y, z)`, galaxy-frame parsecs -- the
                           camera's current view center.
        radius_pc (float): The requested view radius, parsecs -- internally
                           clamped to `PLANNED_RADIUS_CAP_PC` (see that
                           constant's own docstring); a caller with a wider
                           view gets no addresses beyond that cap.
        edge_pc (float): The sector edge length, parsecs.
        shape (galaxyDensity.GalaxyShape or None): The galaxy's stored
            shape parameters. `None` (no `generate.py plan` has been run
            against this database yet) means there's no density model to
            qualify slots against, so every enumerated address is
            returned unfiltered (predicted star count/density fields are
            `None` on each entry) -- addressable, even if this galaxy's
            real content hasn't been planned out yet.
        expected_system_count_at_density_1 (float or None): See
            `galaxySkeleton.expected_system_count_at_density_1` -- required
            (and used) only when `shape` is given.
        exclude_addresses (set): `{(ring_index, layer_index, ring_slot_index), ...}`
            -- addresses to skip because a real sector already exists
            there (the caller already queried those separately as the
            "placed" tier; this avoids listing the same address twice
            under two different tiers).
        cap (int): See `PLANNED_MAX_RESULTS`.

    Returns:
        list[dict]: Closest-first, each with `ring_index`, `layer_index`,
            `ring_slot_index`, `x`/`y`/`z` (parsecs), `distance_pc`
            (from `center_pc`), `designation`
            (`provisional_sector_designation`), and `predicted_star_count`/
            `relative_density` (both `None` when `shape` is `None`).
    """
    r = min(radius_pc, PLANNED_RADIUS_CAP_PC)
    # heapq.nsmallest over a generator keeps only `cap` candidates alive at
    # once, and the per-entry dict (with its designation string) is built
    # only for the ones actually returned.
    closest = heapq.nsmallest(
        cap,
        _qualifying_slots(
            enumerate_sectors_within_radius(center_pc, r, edge_pc),
            shape, expected_system_count_at_density_1, exclude_addresses,
        ),
        key=lambda slot: (slot[6], slot[0], slot[1], slot[2]),
    )
    return [_planned_entry(slot, distance_pc=slot[6]) for slot in closest]


def _qualifying_slots(slots, shape, expected_system_count_at_density_1, exclude_addresses):
    """Filters `enumerate_sectors_within_radius`-shaped tuples down to the
    ones worth showing as "planned" (see `planned_slots_in_view`), yielding
    `(ring_index, layer_index, ring_slot_index, x, y, z, distance_pc,
    star_count, density)` -- the last two `None` when `shape` is `None`."""
    threshold = qualifying_threshold_star_count()
    for ring_index, layer_index, slot_index, x, y, z, distance_pc in slots:
        if (ring_index, layer_index, slot_index) in exclude_addresses:
            continue
        if shape is not None:
            star_count = predicted_star_count((x, y, z), shape, expected_system_count_at_density_1)
            if star_count < threshold:
                continue
            density = relative_density((x, y, z), shape)
        else:
            star_count = None
            density = None
        yield (ring_index, layer_index, slot_index, x, y, z, distance_pc, star_count, density)


def _planned_entry(slot, distance_pc=None):
    """One planned-tier dict from a `_qualifying_slots` tuple."""
    ring_index, layer_index, slot_index, x, y, z, _distance, star_count, density = slot
    entry = {
        "ring_index": ring_index, "layer_index": layer_index, "ring_slot_index": slot_index,
        "x": x, "y": y, "z": z,
        "designation": provisional_sector_designation(ring_index, layer_index, slot_index),
        "predicted_star_count": star_count, "relative_density": density,
    }
    if distance_pc is not None:
        entry["distance_pc"] = distance_pc
    return entry


# ---------------------------------------------------------------------
# Cube tiles -- the map-tile model behind GET /api/galaxy/tiles.
#
# Instead of asking "everything within R of this point" on every camera
# move (unbounded work, nothing reusable between moves), the 3D map asks
# for fixed cubes of space, the way a web map asks for fixed image tiles.
# Space is an octree: level 0 is one cube `TILE_ROOT_EDGE_PC` on a side,
# centered on the galactic origin, and each level halves the edge.
# A tile's key is `"level/ix/iy/iz"`, where `ix` counts cubes along x from
# the root cube's -x face (same for y/z). A tile's contents depend only on
# its key and the database, so the browser and the web layer can both
# cache it (see `html/lib/tilecache.py` and `static/galaxymap3d.js`).
#
# Each tile's work is bounded by construction: placed sectors are capped
# per tile (`queryDb.GALAXY_TILE_MAX_PLACED`), planned slots are only
# enumerated for small tiles (`PLANNED_TILE_MAX_EDGE_PC`), and density
# clouds are a fixed point count.
# ---------------------------------------------------------------------

TILE_ROOT_EDGE_PC = 65536.0
"""float: Edge of the level-0 cube, parsecs -- comfortably larger than
any galaxy this project generates (the default is 15,000 pc in radius),
and a power of two so every level's edge is an exact float."""

TILE_MAX_LEVEL = 12
"""int: Finest tile level -- `TILE_ROOT_EDGE_PC / 2**12` = 16 pc."""

PLANNED_TILE_MAX_EDGE_PC = 16.0
"""float: Planned slots are only listed in tiles at most this big. A
16 pc cube holds about 100 slots of the default 11.5 ly sector size; a
bigger tile would mean far more dots than a view at that zoom can show."""

PLANNED_MAX_SLOTS_PER_TILE = 250
"""int: Safety net for a galaxy with a much smaller sector edge than the
default: a tile that could hold more slots than this lists none, rather
than enumerating an unbounded number."""


def tile_key(level, ix, iy, iz):
    """The canonical `"level/ix/iy/iz"` string for a tile."""
    return f"{level}/{ix}/{iy}/{iz}"


def parse_tile_key(key):
    """
    Parses and validates a `"level/ix/iy/iz"` tile key.

    Args:
        key (str): The tile key.

    Returns:
        tuple: `(level, ix, iy, iz)`, all ints.

    Raises:
        ValueError: If `key` is malformed or out of range.
    """
    parts = str(key).split("/")
    if len(parts) != 4:
        raise ValueError(f"tile key {key!r} must look like level/ix/iy/iz")
    try:
        level, ix, iy, iz = (int(part) for part in parts)
    except ValueError:
        raise ValueError(f"tile key {key!r} must be four integers")
    if not 0 <= level <= TILE_MAX_LEVEL:
        raise ValueError(f"tile level must be 0..{TILE_MAX_LEVEL}, got {level}")
    span = 2 ** level
    if not all(0 <= index < span for index in (ix, iy, iz)):
        raise ValueError(f"tile indices at level {level} must be 0..{span - 1}")
    return level, ix, iy, iz


def tile_edge_pc(level):
    """Edge length of a level-`level` tile, parsecs."""
    return TILE_ROOT_EDGE_PC / (2 ** level)


def tile_bounds_pc(level, ix, iy, iz):
    """
    A tile's half-open box, `[lo, hi)` on each axis, parsecs.

    Returns:
        tuple: `((x_lo, y_lo, z_lo), (x_hi, y_hi, z_hi))`.
    """
    edge = tile_edge_pc(level)
    origin = -TILE_ROOT_EDGE_PC / 2.0
    lo = (origin + ix * edge, origin + iy * edge, origin + iz * edge)
    hi = (lo[0] + edge, lo[1] + edge, lo[2] + edge)
    return lo, hi


def tile_keys_containing(point_pc):
    """
    The key of the one tile per level whose half-open box holds
    `point_pc` -- every tile a sector centered there appears in, so the
    tiles to refetch when that sector changes (`queryDb.galaxy_changes`).

    Args:
        point_pc (tuple): `(x, y, z)`, parsecs.

    Returns:
        list[str]: One key per level, coarsest first, or `[]` for a point
            outside the root cube (which no tile holds).
    """
    origin = -TILE_ROOT_EDGE_PC / 2.0
    keys = []
    for level in range(TILE_MAX_LEVEL + 1):
        edge = tile_edge_pc(level)
        span = 2 ** level
        index = []
        for value in point_pc:
            i = math.floor((value - origin) / edge)
            # Settle a rounding tie against the same `lo <= v < hi` test
            # `queryDb.galaxy_sectors_in_box` uses.
            if origin + i * edge > value:
                i -= 1
            elif origin + (i + 1) * edge <= value:
                i += 1
            if not 0 <= i < span:
                return []
            index.append(i)
        keys.append(tile_key(level, *index))
    return keys


def tile_level_for_view_radius(radius_pc):
    """
    The level whose tiles are the smallest still at least `radius_pc` on
    a side, so a view sphere of that radius touches at most 3 tiles along
    each axis (27 in all). `static/galaxymap3d.js` computes the same thing
    client-side; this is the Python twin the `/galaxy` page uses for the
    first frame.
    """
    if not radius_pc > 0:  # also NaN
        return TILE_MAX_LEVEL
    if math.isinf(radius_pc):
        return 0
    ratio = TILE_ROOT_EDGE_PC / radius_pc
    if math.isinf(ratio):  # a subnormal radius (MAP.90): the finest level
        return TILE_MAX_LEVEL
    level = math.floor(math.log2(ratio))
    return max(0, min(TILE_MAX_LEVEL, level))


def tiles_intersecting_sphere(level, center_pc, radius_pc):
    """
    Every tile key at `level` whose box comes within `radius_pc` of
    `center_pc` -- what a view of that radius needs.

    Returns:
        list[str]: Tile keys, nearest first.
    """
    edge = tile_edge_pc(level)
    origin = -TILE_ROOT_EDGE_PC / 2.0
    span = 2 ** level
    ranges = []
    for axis in range(3):
        first = math.floor((center_pc[axis] - radius_pc - origin) / edge)
        last = math.floor((center_pc[axis] + radius_pc - origin) / edge)
        ranges.append(range(max(0, first), min(span - 1, last) + 1))

    found = []
    radius_sq = radius_pc * radius_pc
    for ix in ranges[0]:
        for iy in ranges[1]:
            for iz in ranges[2]:
                lo, hi = tile_bounds_pc(level, ix, iy, iz)
                distance_sq = 0.0
                for axis in range(3):
                    nearest = min(max(center_pc[axis], lo[axis]), hi[axis])
                    distance_sq += (center_pc[axis] - nearest) ** 2
                if distance_sq <= radius_sq:
                    found.append((distance_sq, tile_key(level, ix, iy, iz)))
    found.sort()
    return [key for _distance, key in found]


def _in_box(point, lo, hi):
    return all(lo[axis] <= point[axis] < hi[axis] for axis in range(3))


def planned_slots_in_tile(level, ix, iy, iz, edge_pc, shape,
                          expected_system_count_at_density_1, exclude_addresses):
    """
    Every planned slot (see `planned_slots_in_view` for what qualifies)
    whose center lies in this tile's half-open box -- so each slot belongs
    to exactly one tile per level. Empty for tiles bigger than
    `PLANNED_TILE_MAX_EDGE_PC`, or when the sector edge is so small the
    tile could hold more than `PLANNED_MAX_SLOTS_PER_TILE` slots.

    Returns:
        list[dict]: Same entries as `planned_slots_in_view`, minus
            `distance_pc` (there's no view center), ordered by address.
    """
    tile_edge = tile_edge_pc(level)
    if tile_edge > PLANNED_TILE_MAX_EDGE_PC:
        return []
    if (tile_edge / edge_pc + 1) ** 3 > PLANNED_MAX_SLOTS_PER_TILE:
        return []

    lo, hi = tile_bounds_pc(level, ix, iy, iz)
    center = tuple((lo[axis] + hi[axis]) / 2.0 for axis in range(3))
    circumradius = tile_edge * math.sqrt(3.0) / 2.0
    in_tile = (
        slot for slot in enumerate_sectors_within_radius(center, circumradius, edge_pc)
        if _in_box(slot[3:6], lo, hi)
    )
    slots = sorted(
        _qualifying_slots(in_tile, shape, expected_system_count_at_density_1, exclude_addresses),
        key=lambda slot: slot[:3],
    )
    return [_planned_entry(slot) for slot in slots]
