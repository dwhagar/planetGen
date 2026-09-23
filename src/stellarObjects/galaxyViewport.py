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

import math
import random

from .galaxyDensity import predicted_star_count, relative_density
from .galaxyGeometry import enumerate_sectors_within_radius, provisional_sector_designation

PLANNED_RADIUS_CAP_PC = 200.0
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
designed for (see `galaxyGeometry`'s own module docstring)."""

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
    results = []
    for shell_index, shell_slot_index, x, y, z, distance_pc in enumerate_sectors_within_radius(
        center_pc, r, edge_pc,
    ):
        if (shell_index, shell_slot_index) in exclude_addresses:
            continue

        if shape is not None:
            star_count = predicted_star_count((x, y, z), shape, expected_system_count_at_density_1)
            if star_count < qualifying_threshold_star_count():
                continue
            density = relative_density((x, y, z), shape)
        else:
            star_count = None
            density = None

        designation = provisional_sector_designation(shell_index, shell_slot_index, edge_pc, edge_ly)
        results.append({
            "shell_index": shell_index, "shell_slot_index": shell_slot_index,
            "x": x, "y": y, "z": z, "distance_pc": distance_pc,
            "designation": designation,
            "predicted_star_count": star_count, "relative_density": density,
        })

    results.sort(key=lambda entry: entry["distance_pc"])
    return results[:cap]


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
MAX_DENSITY_SAMPLE_ATTEMPTS_FACTOR` mixture-proposal candidates before
falling back to `_sample_point_in_sphere` for whatever's still missing --
see that function's own docstring for why a view far from where the
mixture proposal actually places mass (rare, e.g. deep in the halo) needs
a bounded fallback rather than an unbounded/empty result."""


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
        if rng.random() < BULGE_SAMPLE_FRACTION:
            x, y, z = _sample_bulge_point_pc(rng, shape)
        else:
            x, y, z = _sample_disk_point_pc(rng, shape)
        dx, dy, dz = x - center_pc[0], y - center_pc[1], z - center_pc[2]
        if dx * dx + dy * dy + dz * dz > radius_sq:
            continue
        points.append({"x": x, "y": y, "z": z, "relative_density": relative_density((x, y, z), shape)})

    # Bounded fallback -- a view far from where the mixture proposal above
    # actually places mass (rare: e.g. a camera panned deep into the
    # halo) can exhaust max_attempts short of `count` accepted points;
    # top up with the old uniform-in-view-sphere draw so the cloud never
    # silently thins out, even though those extra points are less likely
    # to land somewhere bright.
    while len(points) < count:
        x, y, z = _sample_point_in_sphere(rng, center_pc, radius_pc)
        points.append({"x": x, "y": y, "z": z, "relative_density": relative_density((x, y, z), shape)})

    return points
