# planetgen/galaxy/nebula_volume.py

"""
A nebula group's shape, sampled from the galaxy's gas and star-forming density (MAP.365)
========================================================================================

Boss (2026-11-04): the nebulae on the Galaxy Map must be irregular, bulbous, wispy shapes taken from the galaxy model's real
gas and star-forming density, never a circle and never a random blob per object ID. This module turns the model's density
(`nebula_density.density_at`) into the points the map draws.

* **An absolute lattice.** The samples sit at whole multiples of one spacing (`spacing_pc(level)`, a power of two times
  `LATTICE_BASE_PC`), measured from the galaxy's origin, never from the nebula or the block. Two groups that ask for the same
  place at the same level get the same point and the same value, so a nebula that crosses a block or sector border has no seam,
  and a finer level refines the coarser one on the same grid.
* **One owner for every point.** A point belongs to the nebula it is nearest to in its own radii (`d / r`), and a block draws
  only the lattice points that lie in it (a nebula that crosses a border is drawn half by each block). So two blocks never draw
  the same place twice (adding their colors would double it), no point is lost where nebulae touch, and the shape needs no home
  block.
* **Soft, uneven edges.** The value of a point is the density near the nebula's centre pushed through a soft threshold at the
  model's `FORM_LEVEL` and a window that fades to nothing at the nebula's radius, so a core is dense, the rim is wispy and the
  edge follows the density, not a circle. A shell family (a planetary nebula, a supernova remnant) uses a shell profile.
* **Detail follows the lattice.** The finer the lattice, the more noise octaves the density is asked for
  (`detail_for_spacing`), so a closer view has more filaments and the coarse blobs stay where they were.

The drawing choices here (the window, the softness of the threshold, the number of samples across a radius) are display
settings, not physics; the physics is the density they read. No database and no three.js here: `db/nebula_volumes.py` reads the
nebulae and `static/nebulavolume.js` draws the points.
"""

import math
from collections import namedtuple

import numpy as np

LATTICE_BASE_PC = 0.25
"""float: The finest lattice spacing, parsecs (level 0); level n is this times 2**n."""

MAX_LEVEL = 10
"""int: The coarsest level (256 pc)."""

SAMPLES_ACROSS_RADIUS = 6.0
"""float: About how many lattice steps a group's radius is, which sets its level: a group is drawn from a few thousand points
whatever its size (the finer groups of a closer node hold the finer detail)."""

MAX_POINTS = 60000
"""int: The most lattice points looked at for one group; a group needing more is raised to a coarser level."""

MIN_STRENGTH = 2
"""int: A point weaker than this (of 255) is left out."""

DETAIL_WAVELENGTH_PC = 100.0
"""float: The wavelength of noise octave 0 in `nebula_density` (`W0`); octave n resolves features of about W0 / 2**(n - 1)."""

SOFT_BELOW = 0.12
SOFT_ABOVE = 0.30
"""floats: The soft threshold on the density: nothing below `FORM_LEVEL - SOFT_BELOW`, full strength from `FORM_LEVEL + SOFT_ABOVE`."""

WINDOW_FROM = 0.6
"""float: A solid nebula's window is 1 inside this share of its radius and fades to 0 at the radius."""

SHELL_PROFILE = (0.35, 0.7, 0.8, 1.0)
"""tuple: A shell nebula's profile in shares of its radius: rises over the first two numbers, falls over the last two."""

SHELL_FAMILIES = ("planetary", "supernova-remnant")
"""tuple: The families drawn as a shell, not a filled ball."""

FIELD_OF_FAMILY = {
    "diffuse": "young", "emission": "young", "reflection": "young", "planetary": "young",
    "dark": "gas", "supernova-remnant": "young",
}
"""dict: Which density a family is drawn from: the dark family (molecular gas) from the gas, the rest from the young stars."""

Volume = namedtuple("Volume", ["level", "spacing_pc", "ix", "iy", "iz", "strength", "family"])
"""The sampled points: `ix`, `iy`, `iz` (int32 lattice steps from the origin, so a point is at `i * spacing_pc` parsecs),
`strength` (uint8, 0 to 255) and `family` (uint8: the index in `FAMILIES` of the nebula that owns the point)."""


def spacing_pc(level):
    """The lattice spacing of `level`, parsecs."""
    return LATTICE_BASE_PC * 2.0 ** int(level)


def level_for_radius(radius_pc):
    """The lattice level for nebulae of about `radius_pc`: `SAMPLES_ACROSS_RADIUS` steps to the radius."""
    if not radius_pc > 0:
        return MAX_LEVEL
    level = round(math.log2(radius_pc / SAMPLES_ACROSS_RADIUS / LATTICE_BASE_PC))
    return int(min(MAX_LEVEL, max(0, level)))


def typical_radius(nebulae):
    """The radius a group's lattice is chosen for: the volume-weighted mean of its members' radii, so the big nebulae that fill
    the picture get their steps and the small ones around them, a few pixels across, get what is left."""
    radii = np.array([n["radius_pc"] for n in nebulae], dtype=np.float64)
    return float((radii ** 4).sum() / (radii ** 3).sum())


def detail_for_spacing(spacing, max_detail):
    """How many noise octaves to ask for at lattice spacing `spacing` (pc): enough to resolve features of two steps."""
    n = 1.0 + math.log2(DETAIL_WAVELENGTH_PC / (2.0 * spacing))
    return int(min(max_detail, max(0, math.ceil(n))))


def _smoothstep(edge0, edge1, x):
    t = np.clip((x - edge0) / (edge1 - edge0), 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def profile(u, shell):
    """The window of a nebula at `u` (the distance from its centre in its radii): 1 in a solid nebula's core, down to 0 at the edge;
    a shell's rises to a ring and falls again."""
    u = np.asarray(u, dtype=np.float64)
    if shell:
        rise0, rise1, fall0, fall1 = SHELL_PROFILE
        return _smoothstep(rise0, rise1, u) * (1.0 - _smoothstep(fall0, fall1, u))
    return 1.0 - _smoothstep(WINDOW_FROM, 1.0, u)


def _lattice_in_spheres(spheres, spacing):
    """The lattice points (an (n, 3) int64 array, no repeats, sorted) inside any of `spheres` (rows `x, y, z, r`)."""
    found = []
    for x, y, z, r in spheres:
        lo = np.floor((np.array([x, y, z]) - r) / spacing).astype(np.int64)
        hi = np.ceil((np.array([x, y, z]) + r) / spacing).astype(np.int64)
        axes = [np.arange(lo[i], hi[i] + 1, dtype=np.int64) for i in range(3)]
        grid = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)
        inside = np.sum((grid * spacing - np.array([x, y, z])) ** 2, axis=1) <= r * r
        found.append(grid[inside])
    if not found:
        return np.zeros((0, 3), dtype=np.int64)
    return np.unique(np.concatenate(found), axis=0)


def _points_needed(spheres, spacing):
    return sum(((2.0 * r / spacing) + 2.0) ** 3 for _x, _y, _z, r in spheres)


def _empty_volume(level):
    return Volume(level, spacing_pc(level), *(np.zeros(0, dtype=t) for t in (np.int32, np.int32, np.int32, np.uint8, np.uint8)))


def sample_group(nebulae, inside, density, form_level, max_detail, families, level=None):
    """
    The points that draw a block's nebulae at its lattice level.

    Args:
        nebulae (list[dict]): Every nebula whose sphere reaches the block: `x`, `y`, `z`, `radius_pc` (parsecs) and `family` (a name in
            `families`). A point is drawn only where one of them holds it (inside its radius); where several do, the one it is
            nearest to in radii owns it.
        inside (callable or None): `inside(points_pc) -> bool array`, which of the (n, 3) lattice points lie in the block; the block's
            own, so two blocks never draw the same place. `None` keeps every point (a test).
        density (callable): `density(points_pc, field, detail) -> array in 0..1`, one value per row of the (n, 3) `points_pc`;
            `field` is `"gas"` or `"young"`. `nebula_density.density_at` with the galaxy's seed and shape filled in.
        form_level (float): The value above which a point is inside a nebula (`nebula_density.FORM_LEVEL`).
        max_detail (int): The most noise octaves the density gives (`nebula_density.MAX_DETAIL`).
        families (tuple[str]): The family names, in the order of the `family` byte of the result.
        level (int, optional): The lattice level; default `level_for_radius` of the nebulae's `typical_radius`, raised until the
            points looked at fit `MAX_POINTS`.

    Returns:
        Volume: The points, in lattice order. Empty arrays for a block no nebula reaches.
    """
    if not nebulae:
        return _empty_volume(0)
    spheres = [(n["x"], n["y"], n["z"], n["radius_pc"]) for n in nebulae]
    if level is None:
        level = level_for_radius(typical_radius(nebulae))
    while level < MAX_LEVEL and _points_needed(spheres, spacing_pc(level)) > MAX_POINTS:
        level += 1
    spacing = spacing_pc(level)
    lattice = _lattice_in_spheres(spheres, spacing)
    if len(lattice) and inside is not None:
        lattice = lattice[np.asarray(inside(lattice * spacing), dtype=bool)]
    if not len(lattice):
        return _empty_volume(level)
    points = lattice * spacing
    centres = np.array([[n["x"], n["y"], n["z"]] for n in nebulae])
    radii = np.array([n["radius_pc"] for n in nebulae], dtype=np.float64)
    nearest = np.zeros(len(points), dtype=np.int64)
    nearest_u = np.full(len(points), np.inf)
    for start in range(0, len(nebulae), 64):    # in slices of nebulae, so the (points x nebulae) array stays small
        block = slice(start, start + 64)
        distance = np.sqrt(((points[:, None, :] - centres[None, block, :]) ** 2).sum(axis=2)) / radii[None, block]
        best = np.argmin(distance, axis=1)
        best_u = distance[np.arange(len(points)), best]
        better = best_u < nearest_u
        nearest[better] = start + best[better]
        nearest_u[better] = best_u[better]
    owned = nearest_u < 1.0
    lattice, points, nearest, nearest_u = lattice[owned], points[owned], nearest[owned], nearest_u[owned]
    if not len(lattice):
        return _empty_volume(level)
    detail = detail_for_spacing(spacing, max_detail)
    owner_family = np.array([families.index(nebulae[i]["family"]) for i in nearest], dtype=np.uint8)
    strength = np.zeros(len(points), dtype=np.float64)
    for field in sorted(set(FIELD_OF_FAMILY[families[f]] for f in set(owner_family.tolist()))):
        pick = np.array([FIELD_OF_FAMILY[families[f]] == field for f in owner_family])
        value = np.asarray(density(points[pick], field, detail), dtype=np.float64)
        strength[pick] = _smoothstep(form_level - SOFT_BELOW, form_level + SOFT_ABOVE, value)
    for family_index in sorted(set(owner_family.tolist())):
        pick = owner_family == family_index
        strength[pick] *= profile(nearest_u[pick], families[family_index] in SHELL_FAMILIES)
    quantised = np.rint(strength * 255.0).astype(np.int64)
    keep = quantised >= MIN_STRENGTH
    return Volume(level, spacing, lattice[keep, 0].astype(np.int32), lattice[keep, 1].astype(np.int32),
                  lattice[keep, 2].astype(np.int32), quantised[keep].astype(np.uint8), owner_family[keep])


MAGIC = b"NVOL"
"""bytes: The start of a packed volume."""

PACK_VERSION = 1


def pack(volume):
    """A volume as bytes: `NVOL`, the version (u8), the level (u8), the count (u32) and the spacing (f32), then the columns
    `ix`, `iy`, `iz` (int32 each), `strength` and `family` (u8 each), little-endian, so it compresses well."""
    n = len(volume.ix)
    head = MAGIC + bytes([PACK_VERSION, volume.level]) + np.array([n], dtype="<u4").tobytes() \
        + np.array([volume.spacing_pc], dtype="<f4").tobytes()
    return head + b"".join(
        np.ascontiguousarray(column).astype(dtype).tobytes()
        for column, dtype in ((volume.ix, "<i4"), (volume.iy, "<i4"), (volume.iz, "<i4"),
                              (volume.strength, "u1"), (volume.family, "u1")))
