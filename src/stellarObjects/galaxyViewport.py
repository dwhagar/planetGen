# stellarObjects/galaxyViewport.py

"""
Interactive 3D Galaxy Map viewport queries.

The flat, server-rendered-once SVG Galaxy Map (`html/lib/galaxymap.py`)
plots every galaxy-placed sector in one shot, since the whole point of
that map is a single bird's-eye overview. The interactive 3D map
(`html/lib/galaxymap3d.py`/`static/galaxymap3d.js`) is the opposite: a
real camera that moves freely through the galaxy, so what it needs to
plot changes every time it moves -- this module is the data layer behind
that, called fresh (via `queryDb.galaxy_view`) each time the camera's
viewport changes.

Three tiers of content, matching the three sprite kinds
`static/galaxymap3d.js` draws:

- **Placed**: real, already-generated sectors -- `queryDb.
  galaxy_sectors_in_view` handles this tier directly (it needs the
  database; nothing here does).
- **Planned**: real, not-yet-generated sector *addresses* -- exact
  `(shell_index, shell_slot_index)` slots this galaxy's own shape model
  predicts would qualify (see `galaxyDensity.predicted_star_count`, the
  same >= 1 threshold `generate.py`'s `ensure_sector_generated`/
  `_BatchDensity` already gate real generation on), enumerated exactly
  via `galaxyGeometry.enumerate_sectors_within_radius` -- `planned_slots_in_view`.
  Deliberately capped to a small view radius (`PLANNED_RADIUS_CAP_PC`):
  see that function's own docstring for why an unbounded radius here
  would be catastrophic.
- **Density**: for a view too wide to enumerate individual planned
  slots, a coarse illustrative point cloud sampled directly from
  `galaxyDensity.relative_density` -- `density_sample_points`. Not real
  addresses (nothing here is clickable/generatable), purely a "here's
  roughly how the predicted spiral looks from this far out" visual,
  resampled fresh (deterministically, for a stable view while the camera
  holds still) on every viewport query so it stays as high-resolution as
  the current view warrants, unlike the flat map's own one-shot density
  grid (`html/lib/galaxymap.py`'s `_expected_density_elements`).

Pure and side-effect-free, like `galaxyGeometry.py`/`galaxyDensity.py` --
no database or I/O; `queryDb.py` owns combining this with the database's
own placed-sector rows.
"""

import heapq
import math
import random

from .galaxyDensity import predicted_star_count, relative_density
from .galaxyGeometry import enumerate_sectors_within_radius, provisional_sector_designation

PLANNED_RADIUS_CAP_PC = 40.0
"""float: The largest view radius `planned_slots_in_view` will actually
enumerate individual slot addresses for, however large a `radius_pc` its
caller asks for. `enumerate_sectors_within_radius` is only cheap for a
"local neighborhood" -- a view spanning the whole galaxy (tens of
thousands of parsecs) would enumerate a shell's entire slot count, which
runs into the hundreds of millions for an outer shell (see
`docs/design/galaxy-coordinate-system.md` section 3's own worked table).
Beyond this radius, callers get `density_sample_points`'s illustrative
cloud instead -- exact addresses only become worth enumerating once the
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

DENSITY_SAMPLE_COUNT = 1200
"""int: How many points `density_sample_points` draws per call -- enough
to read as a recognizable cloud (comparable to the flat map's own
`_DENSITY_GRID_CELLS ** 2` tile budget, minus the ~21% corner-crop) without
being an unreasonable per-request payload for something purely
illustrative."""


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


def planned_slots_in_view(center_pc, radius_pc, edge_pc, edge_ly, shape,
                           expected_system_count_at_density_1, exclude_addresses,
                           cap=PLANNED_MAX_RESULTS):
    """
    Every real, not-yet-generated sector address within `radius_pc` of
    `center_pc` that this galaxy's own shape model predicts would qualify
    (`predicted_star_count >= 1`) -- the "planned" tier `static/
    galaxymap3d.js` renders as small, dim, clickable dots (see the module
    docstring's own tier breakdown), each carrying its real
    `(shell_index, shell_slot_index)` address and a copyable provisional
    designation, ready to feed straight into `generate.py galaxy --shell
    K --slot N`.

    Args:
        center_pc (tuple): `(x, y, z)`, galaxy-frame parsecs -- the
                           camera's current view center.
        radius_pc (float): The requested view radius, parsecs -- internally
                           clamped to `PLANNED_RADIUS_CAP_PC` (see that
                           constant's own docstring); a caller with a wider
                           view should use `density_sample_points` instead
                           for anything beyond that cap.
        edge_pc (float): The sector edge length, parsecs.
        edge_ly (float): The same edge length, light-years (for
                         `provisional_sector_designation`'s own Ring
                         lookup -- see that function's docstring for why
                         it takes both units rather than converting one
                         from the other itself).
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
        exclude_addresses (set): `{(shell_index, shell_slot_index), ...}`
            -- addresses to skip because a real sector already exists
            there (the caller already queried those separately as the
            "placed" tier; this avoids listing the same address twice
            under two different tiers).
        cap (int): See `PLANNED_MAX_RESULTS`.

    Returns:
        list[dict]: Closest-first, each with `shell_index`,
            `shell_slot_index`, `x`/`y`/`z` (parsecs), `distance_pc`
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
        key=lambda slot: (slot[5], slot[0], slot[1]),
    )
    return [
        _planned_entry(slot, edge_pc, edge_ly, distance_pc=slot[5])
        for slot in closest
    ]


def _qualifying_slots(slots, shape, expected_system_count_at_density_1, exclude_addresses):
    """Filters `enumerate_sectors_within_radius`-shaped tuples down to the
    ones worth showing as "planned" (see `planned_slots_in_view`), yielding
    `(shell_index, shell_slot_index, x, y, z, distance_pc, star_count,
    density)` -- the last two `None` when `shape` is `None`."""
    threshold = qualifying_threshold_star_count()
    for shell_index, shell_slot_index, x, y, z, distance_pc in slots:
        if (shell_index, shell_slot_index) in exclude_addresses:
            continue
        if shape is not None:
            star_count = predicted_star_count((x, y, z), shape, expected_system_count_at_density_1)
            if star_count < threshold:
                continue
            density = relative_density((x, y, z), shape)
        else:
            star_count = None
            density = None
        yield (shell_index, shell_slot_index, x, y, z, distance_pc, star_count, density)


def _planned_entry(slot, edge_pc, edge_ly, distance_pc=None):
    """One planned-tier dict from a `_qualifying_slots` tuple."""
    shell_index, shell_slot_index, x, y, z, _distance, star_count, density = slot
    entry = {
        "shell_index": shell_index, "shell_slot_index": shell_slot_index,
        "x": x, "y": y, "z": z,
        "designation": provisional_sector_designation(shell_index, shell_slot_index, edge_pc, edge_ly),
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

DENSITY_TILE_SAMPLE_COUNT = 1600
"""int: Points in one density cloud (`density_points_for_tile`). A density
cloud covers twice its anchor tile's edge in radius (see that function),
so it spreads over more space than the old per-view cloud and gets a few
more points than `DENSITY_SAMPLE_COUNT` to stay about as dense on screen."""


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
    client-side; this is the Python twin `html/galaxy.py` uses for the
    first frame.
    """
    if radius_pc <= 0:
        return TILE_MAX_LEVEL
    level = math.floor(math.log2(TILE_ROOT_EDGE_PC / radius_pc))
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


def planned_slots_in_tile(level, ix, iy, iz, edge_pc, edge_ly, shape,
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
        if _in_box(slot[2:5], lo, hi)
    )
    slots = sorted(
        _qualifying_slots(in_tile, shape, expected_system_count_at_density_1, exclude_addresses),
        key=lambda slot: (slot[0], slot[1]),
    )
    return [_planned_entry(slot, edge_pc, edge_ly) for slot in slots]


def density_points_for_tile(level, ix, iy, iz, shape, count=DENSITY_TILE_SAMPLE_COUNT):
    """
    The illustrative density cloud (see `density_sample_points`) anchored
    on one tile: sampled within twice the tile's edge of its center, which
    covers any view of radius up to one tile edge centered anywhere inside
    the tile. The map uses the tile holding its camera target at the
    current zoom level (`tile_level_for_view_radius`), so the cloud only
    changes when the target crosses into another tile or the zoom level
    changes -- and each cloud is cacheable by its tile key.

    Returns:
        list[dict]: `density_sample_points`' entries.
    """
    lo, hi = tile_bounds_pc(level, ix, iy, iz)
    center = tuple((lo[axis] + hi[axis]) / 2.0 for axis in range(3))
    return density_sample_points(center, 2.0 * tile_edge_pc(level), shape, count=count)


def _sample_point_in_sphere(rng, center_pc, radius_pc):
    """One point drawn uniformly by volume from the sphere of `radius_pc`
    around `center_pc` -- the standard cube-root-of-uniform radial draw
    (matching `spaceSector._random_unit_direction`'s own "uniform area,
    not uniform angle" concern, extended to volume here) paired with a
    uniform-on-the-sphere direction."""
    cx, cy, cz = center_pc
    u = rng.random()
    r = radius_pc * (u ** (1.0 / 3.0))
    theta = rng.uniform(0.0, 2.0 * math.pi)
    cos_phi = rng.uniform(-1.0, 1.0)
    sin_phi = math.sqrt(max(0.0, 1.0 - cos_phi * cos_phi))
    return (
        cx + r * sin_phi * math.cos(theta),
        cy + r * sin_phi * math.sin(theta),
        cz + r * cos_phi,
    )


BULGE_SAMPLE_FRACTION = 0.18
"""float: Fraction of `density_sample_points`' draws taken from the bulge
proposal (`_sample_bulge_point_pc`) rather than the disk proposal
(`_sample_disk_point_pc`) -- a fixed mixture weight, not derived from the
shape's own `bulge_amplitude` (which would need integrating both terms'
real mass over volume, more precision than a purely illustrative cloud
needs). Chosen so a real Milky-Way-scale shape's bright core still reads
as a visible, denser clump without swamping the disk's own point budget."""

MAX_DENSITY_SAMPLE_ATTEMPTS_FACTOR = 40
"""int: `density_sample_points` draws at most `count *
MAX_DENSITY_SAMPLE_ATTEMPTS_FACTOR` candidates from either sampler before
falling back to plain uniform points for whatever's still missing -- see
that function's own docstring."""

DENSITY_PILOT_ATTEMPTS = 2000
"""int: How many mixture-proposal candidates `density_sample_points`
draws before judging whether that proposal can fill the view at all.
A view small next to the galaxy (a few hundred parsecs out in the disk)
keeps only a handful of galaxy-wide candidates, so it switches to
`_sample_local_density_points` instead."""


def _sample_disk_point_pc(rng, shape):
    """
    One point drawn from an exponential-radial / Laplace-vertical
    proposal shaped like this galaxy's own disk envelope
    (`disk_scale_length_pc`/`disk_scale_height_pc`) -- deliberately blind
    to spiral-arm structure in *placement* (arm contrast still comes
    through each accepted point's own real `relative_density`-driven
    color in `density_sample_points`, exactly the way the flat map's own
    `_DENSITY_ARM_CONTRAST` layered arm shading on top of a smooth radial
    glow rather than trying to bias tile placement by arm). Standard
    inverse-CDF sampling for both the exponential radial and Laplace
    vertical distributions.

    Args:
        rng (random.Random): This call's own seeded generator.
        shape (galaxyDensity.GalaxyShape): The galaxy's shape parameters.

    Returns:
        tuple: `(x, y, z)`, galaxy-frame parsecs, centered on the
              *galactic origin* (not any view center -- the real mass
              this proposal approximates is fixed there regardless of
              where a viewer's camera happens to be looking).
    """
    r_cyl = -shape.disk_scale_length_pc * math.log(max(1e-12, 1.0 - rng.random()))
    theta = rng.uniform(0.0, 2.0 * math.pi)
    v = rng.random() - 0.5
    z = -shape.disk_scale_height_pc * math.copysign(1.0, v) * math.log(max(1e-12, 1.0 - 2.0 * abs(v)))
    return (r_cyl * math.cos(theta), r_cyl * math.sin(theta), z)


def _sample_bulge_point_pc(rng, shape):
    """
    One point drawn from an isotropic exponential-radius proposal shaped
    like this galaxy's own bulge envelope (`bulge_scale_radius_pc`) --
    an approximation (a true 3D exponential-density-profile draw needs an
    `r^2` Jacobian correction this skips), acceptable for a purely
    illustrative point cloud rather than a physically exact sampler.

    Args:
        rng (random.Random): This call's own seeded generator.
        shape (galaxyDensity.GalaxyShape): The galaxy's shape parameters.

    Returns:
        tuple: `(x, y, z)`, galaxy-frame parsecs, centered on the
              galactic origin.
    """
    r = -shape.bulge_scale_radius_pc * math.log(max(1e-12, 1.0 - rng.random()))
    theta = rng.uniform(0.0, 2.0 * math.pi)
    cos_phi = rng.uniform(-1.0, 1.0)
    sin_phi = math.sqrt(max(0.0, 1.0 - cos_phi * cos_phi))
    return (r * sin_phi * math.cos(theta), r * sin_phi * math.sin(theta), r * cos_phi)


def _sample_local_density_points(rng, center_pc, radius_pc, shape, count, max_attempts):
    """
    Up to `count` points inside the view sphere, drawn by rejection
    sampling against the real `relative_density`: uniform candidates,
    each kept with probability `density / ceiling`. The ceiling is the
    highest density among a first batch of `count` candidates, so a view
    where density barely varies (deep in the disk) keeps nearly every
    candidate, and one straddling the disk's edge keeps mostly the dense
    side. Stops after `max_attempts` candidates.

    Returns:
        list[dict]: `density_sample_points`' entries.
    """
    probe = []
    for _ in range(count):
        point = _sample_point_in_sphere(rng, center_pc, radius_pc)
        probe.append((point, relative_density(point, shape)))
    ceiling = max(density for _point, density in probe)
    if ceiling <= 0:
        return []

    points = []

    def offer(point, density):
        if rng.random() * ceiling < density:
            points.append({"x": point[0], "y": point[1], "z": point[2], "relative_density": density})

    for point, density in probe:
        offer(point, density)
    attempts = len(probe)
    while len(points) < count and attempts < max_attempts:
        attempts += 1
        point = _sample_point_in_sphere(rng, center_pc, radius_pc)
        offer(point, relative_density(point, shape))
    return points[:count]


def density_sample_points(center_pc, radius_pc, shape, count=DENSITY_SAMPLE_COUNT):
    """
    A coarse, illustrative point cloud of this galaxy's real predicted
    density within `radius_pc` of `center_pc` -- the "density" tier for a
    view too wide to enumerate individual planned-slot addresses
    (`planned_slots_in_view`'s own `PLANNED_RADIUS_CAP_PC`), see the
    module docstring's tier breakdown.

    Points are drawn by **importance sampling** from a bulge+disk mixture
    proposal shaped like the galaxy's own real mass distribution
    (`_sample_bulge_point_pc`/`_sample_disk_point_pc`, mixed by
    `BULGE_SAMPLE_FRACTION`), each candidate kept only if it lands within
    `radius_pc` of `center_pc` -- not a uniform draw over the query
    volume. A view spanning thousands of parsecs is overwhelmingly empty
    halo by volume, so a uniform draw would waste nearly its whole budget
    out there and read as a sparse, shapeless scatter rather than a
    galaxy; sampling from where the real mass actually concentrates (a
    thin disk plus a bright core) is what makes the resulting cloud read
    as a recognizable bulge+disk shape at a glance, with each point's own
    real `relative_density` (which *does* include the spiral-arm term)
    still driving its brightness/color, so arm structure still shows
    through as contrast within that shape.

    That proposal only works for a view big enough to catch a fair share
    of the galaxy's mass. A zoomed-in view (a few hundred parsecs out in
    the disk) keeps almost none of its candidates, so when a short pilot
    run (`DENSITY_PILOT_ATTEMPTS`) projects it can't fill `count` within
    budget, the view is sampled locally instead
    (`_sample_local_density_points`): uniform candidates inside the view,
    each kept in proportion to its real `relative_density`. That keeps a
    zoomed-in cloud shaped like the disk around it (a slab thinning away
    from the plane) rather than the uniform ball the old top-up drew.

    Deterministically seeded from `(center_pc, radius_pc)` (rounded to
    avoid reseeding on floating-point noise) rather than a fresh random
    seed per call, so repeated queries at the same view (e.g. a debounced
    re-fetch that lands while the camera is holding still) render the same
    cloud instead of visibly flickering between independent random draws.

    Args:
        center_pc (tuple): `(x, y, z)`, galaxy-frame parsecs.
        radius_pc (float): The view radius to sample within, parsecs.
        shape (galaxyDensity.GalaxyShape or None): The galaxy's stored
            shape parameters -- `None` means there's nothing to sample a
            density from, so this returns an empty list (the caller falls
            back to whatever illustrative default it already has, same as
            the flat map's own `galaxy_shape is None` fallback).
        count (int): See `DENSITY_SAMPLE_COUNT`.

    Returns:
        list[dict]: Each with `x`/`y`/`z` (parsecs) and `relative_density`
            (`galaxyDensity.relative_density` at that point) -- empty if
            `shape` is `None` or `radius_pc <= 0`. Always exactly `count`
            long otherwise (the bounded uniform-sphere fallback below
            tops up whatever the mixture proposal couldn't fill within
            its own attempt budget).
    """
    if shape is None or radius_pc <= 0:
        return []

    seed = (round(center_pc[0], 1), round(center_pc[1], 1), round(center_pc[2], 1), round(radius_pc, 1))
    rng = random.Random(str(seed))
    radius_sq = radius_pc * radius_pc

    points = []
    attempts = 0
    max_attempts = count * MAX_DENSITY_SAMPLE_ATTEMPTS_FACTOR
    while len(points) < count and attempts < max_attempts:
        attempts += 1
        if attempts == DENSITY_PILOT_ATTEMPTS and len(points) * max_attempts < count * attempts:
            # This view is too small for the galaxy-wide proposal; sample
            # it locally instead (see the docstring).
            points = _sample_local_density_points(rng, center_pc, radius_pc, shape, count, max_attempts)
            break
        if rng.random() < BULGE_SAMPLE_FRACTION:
            x, y, z = _sample_bulge_point_pc(rng, shape)
        else:
            x, y, z = _sample_disk_point_pc(rng, shape)
        dx, dy, dz = x - center_pc[0], y - center_pc[1], z - center_pc[2]
        if dx * dx + dy * dy + dz * dz > radius_sq:
            continue
        points.append({"x": x, "y": y, "z": z, "relative_density": relative_density((x, y, z), shape)})

    # Bounded fallback -- if neither sampler filled `count` within its
    # budget (a view almost entirely in near-empty halo), top up with
    # uniform points in the view so the cloud never silently thins out.
    while len(points) < count:
        x, y, z = _sample_point_in_sphere(rng, center_pc, radius_pc)
        points.append({"x": x, "y": y, "z": z, "relative_density": relative_density((x, y, z), shape)})

    return points
