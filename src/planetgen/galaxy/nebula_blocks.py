"""
planetgen/galaxy/nebula_blocks.py

Nebula aggregates for the Galaxy Map's block ladder (MAP.172, docs/design/interface/nebula-map-visibility.md sections 3.1
and 5.1). Every block of the ladder (`galaxy/drill.py`, levels 243, 81, 27, 9 and 3) keeps what the nebulae in it
add up to, so the map can shade a block by how much of it is gas and group the clouds of a block into one region,
for blocks no sector has been generated in too: the numbers come from the seed's nebulae, not from the sectors.

A nebula is a sphere (centre, radius, class). Per block:

* `nebula_count`, `volume_pc3`, the volume-weighted `centroid`, the `bound_radius_pc` that holds every member
  around that centroid, the volume of each family and the class holding the most volume (`dominant_class`): of
  the nebulae whose *centre* lies in the block, so a nebula is in one block and the counts add up through the
  ladder.
* `cover`: how much of the block's volume lies inside a nebula. This takes the part of each sphere that lies *in*
  the block, not the whole sphere, so a cloud wider than a block shades every block it reaches. Overlaps are not
  looked for one pair at a time: the summed volume `cover_volume_pc3` (which adds up exactly, block by block and
  nebula by nebula) turns into a share the way randomly placed clouds cover a volume, `1 - exp(-v / V)`. That is
  within a few percent of the true union where the field is 11% full (the design's own figure) and cannot pass 1.

The pure rules live here (numpy over arrays of spheres); `db/nebula_blocks.py` stores them.
"""

import math
from collections import namedtuple
from collections.abc import Sequence

import numpy as np

from planetgen import tuning
from planetgen.galaxy import drill

LEVELS = drill.STAGE_LEVELS[:-1]
"""tuple: The levels kept: every stage block size above a sector (`drill.STAGE_LEVELS`; the map's nodes above them hold no nebula cover yet)."""

FAMILIES = ("diffuse", "emission", "reflection", "planetary", "dark", "supernova-remnant")
"""tuple: The nebula families, in the order of the table's family columns."""

CLASS_LETTERS = tuple(sorted(tuning.NEBULA_CLASSES))
"""tuple: Every nebula class letter, sorted: the index of a class in the working arrays."""

COVER_SAMPLES = (64, 256, 1024)
"""tuple: How many points a sphere's volume is split by among the level-3 blocks it reaches: the smallest number that
`_samples_for` finds enough for its size. A sphere that lies wholly in one block puts every point in it, so the points
only matter for a cloud that crosses a block's edge; a bigger cloud crosses more of them."""

_CHUNK = 20000
"""int: Spheres of the smallest sample size worked on at a time, so the arrays of sample points stay small."""

_BAND_OFFSET = 1 << 21
_WEDGE_BITS = 22
_SHELL_SHIFT = 44

Block = drill.DrillBlock

LevelCover = namedtuple("LevelCover", ["keys", "totals"])
"""One level's cover (PERF.193): the packed keys of the blocks some sphere reaches (`pack`, an int64 array, sorted and
without repeats) and the cubic parsecs of sphere volume in each (a float64 array, the same positions)."""

LevelRows = namedtuple("LevelRows", ["blocks", "classes"])
"""One level's rows: `blocks` (a `BlockRows`, dicts of the block table's columns read like a list) and `classes`
(`(shell, band, wedge, class, volume)` tuples of the class table)."""


def _family_of(letter):
    return tuning.NEBULA_CLASSES[letter]["family"]


_FAMILY_INDEX_OF_CLASS = np.array([FAMILIES.index(_family_of(letter)) for letter in CLASS_LETTERS], dtype=np.int64)


def class_codes(letters):
    """The working-array index of each class letter in `letters`."""
    lookup = {letter: index for index, letter in enumerate(CLASS_LETTERS)}
    return np.array([lookup[letter] for letter in letters], dtype=np.int64)


def _unit_ball_points(count):
    """`count` points spread evenly through the unit ball, the same every call (a three-dimensional additive
    recurrence sequence, so the first `n` points are well spread for any `n`)."""
    g = 1.2207440846057596    # the root of x**4 = x + 1, which gives the 3-D R-sequence
    alpha = np.array([1 / g, 1 / g ** 2, 1 / g ** 3])
    u = (0.5 + alpha[None, :] * np.arange(1, count + 1)[:, None]) % 1.0
    cos_theta = 2 * u[:, 0] - 1
    sin_theta = np.sqrt(1 - cos_theta ** 2)
    phi = 2 * math.pi * u[:, 1]
    radius = np.cbrt(u[:, 2])
    return np.column_stack((radius * sin_theta * np.cos(phi), radius * sin_theta * np.sin(phi), radius * cos_theta))


_BALL = _unit_ball_points(max(COVER_SAMPLES))


def pack(shell, band, wedge):
    """One integer for each `(shell, band, wedge)` of arrays, so blocks can be grouped."""
    return (shell << _SHELL_SHIFT) | ((band + _BAND_OFFSET) << _WEDGE_BITS) | wedge


def unpack(key):
    """`pack`'s inverse for one integer: `(shell, band, wedge)`."""
    key = int(key)
    return (key >> _SHELL_SHIFT, ((key >> _WEDGE_BITS) & ((1 << (_SHELL_SHIFT - _WEDGE_BITS)) - 1)) - _BAND_OFFSET,
            key & ((1 << _WEDGE_BITS) - 1))


def unpack_array(keys):
    """`unpack` for an array of packed keys: `(shell, band, wedge)` arrays."""
    keys = np.asarray(keys, dtype=np.int64)
    return (keys >> _SHELL_SHIFT, ((keys >> _WEDGE_BITS) & ((1 << (_SHELL_SHIFT - _WEDGE_BITS)) - 1)) - _BAND_OFFSET,
            keys & ((1 << _WEDGE_BITS) - 1))


def locate(x, y, z, m, edge_pc):
    """
    The block of size `m` holding each point (arrays of parsecs): `(shell, band, wedge)` arrays. The same cut as
    `drill.drill_block_bounds` draws: block shell `s` holds radii from `(s*m - 1/2)` to `((s+1)*m - 1/2)` edges,
    the row `floor(colatitude * rows / pi)` of the shell's odd row count, and the wedge `floor(longitude * wedges
    / 2pi)` of that row's wedge count.
    """
    radius = np.sqrt(x * x + y * y + z * z)
    shell = np.floor(radius / edge_pc + 0.5).astype(np.int64) // m
    shells, shell_inverse = np.unique(shell, return_inverse=True)
    rows = np.array([drill.drill_band_count(m, int(s)) for s in shells], dtype=np.int64)[shell_inverse]
    colatitude = np.arctan2(np.hypot(x, y), z)
    row = np.clip(np.floor(colatitude * rows / math.pi).astype(np.int64), 0, rows - 1)
    band = (rows - 1) // 2 - row
    pairs, pair_inverse = np.unique((shell << 24) | row, return_inverse=True)
    wedges = np.empty(len(pairs), dtype=np.int64)
    for i, pair in enumerate(pairs):
        pair_shell = int(pair >> 24)
        pair_rows = drill.drill_band_count(m, pair_shell)
        wedges[i] = drill.drill_wedge_count(m, pair_shell, (pair_rows - 1) // 2 - int(pair & ((1 << 24) - 1)))
    wedges = wedges[pair_inverse]
    longitude = np.arctan2(y, x) % (2 * math.pi)
    wedge = np.minimum(wedges - 1, np.floor(longitude * wedges / (2 * math.pi)).astype(np.int64))
    return shell, band, wedge


def block_volume_pc3(m, shell, band, wedge, edge_pc):
    """The volume of block `(m, shell, band, wedge)`, cubic parsecs, from its analytic bounds."""
    (inner, outer), (north, south), (start, end) = drill.drill_block_bounds(Block(m, shell, band, wedge), edge_pc)
    return (outer ** 3 - inner ** 3) / 3 * (math.cos(north) - math.cos(south)) * (end - start)


def cover_of(cover_volume_pc3, volume_pc3):
    """The share of a block of `volume_pc3` that summed nebula volume `cover_volume_pc3` covers (see the module)."""
    if volume_pc3 <= 0:
        return 0.0
    return -math.expm1(-cover_volume_pc3 / volume_pc3)


def sphere_volumes(radius):
    """The volumes of spheres of `radius` (array, parsecs)."""
    return 4 / 3 * math.pi * radius ** 3


def _centre_keys(x, y, z, m, edge_pc):
    keys = np.empty(len(x), dtype=np.int64)
    for start in range(0, len(x), 200000):
        end = start + 200000
        keys[start:end] = pack(*locate(x[start:end], y[start:end], z[start:end], m, edge_pc))
    return keys


def _samples_for(radius, edge_pc):
    """The number of points (`COVER_SAMPLES`) each radius (array, parsecs) is split with: 64 up to one sector edge, 256
    up to two, 1024 beyond."""
    return np.where(radius <= edge_pc, COVER_SAMPLES[0], np.where(radius <= 2 * edge_pc, COVER_SAMPLES[1],
                                                                  COVER_SAMPLES[2]))


def fine_cover(x, y, z, radius, edge_pc):
    """
    How much sphere volume lies in each level-3 block the spheres reach: a `LevelCover` of the packed blocks and
    their cubic parsecs, each sphere's volume split by where evenly spread points inside it land. The bigger blocks'
    cover is these sums added up (`roll_up`), so it adds up through the ladder exactly.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    radius = np.asarray(radius, dtype=float)
    volume = sphere_volumes(radius)
    samples_of = _samples_for(radius, edge_pc)
    keys_out = []
    weights_out = []
    for samples in COVER_SAMPLES:
        which = np.flatnonzero(samples_of == samples)
        ball = _BALL[:samples]
        chunk = max(1, _CHUNK * COVER_SAMPLES[0] // samples)
        for start in range(0, len(which), chunk):
            part = which[start:start + chunk]
            px = (x[part, None] + radius[part, None] * ball[None, :, 0]).ravel()
            py = (y[part, None] + radius[part, None] * ball[None, :, 1]).ravel()
            pz = (z[part, None] + radius[part, None] * ball[None, :, 2]).ravel()
            keys = pack(*locate(px, py, pz, drill.DRILL_STEP, edge_pc))
            unique, inverse = np.unique(keys, return_inverse=True)
            keys_out.append(unique)
            weights_out.append(np.bincount(inverse, weights=np.repeat(volume[part] / samples, samples),
                                           minlength=len(unique)))
    if not keys_out:
        return LevelCover(np.zeros(0, dtype=np.int64), np.zeros(0))
    return _summed(np.concatenate(keys_out), np.concatenate(weights_out))


def _summed(keys, weights):
    """A `LevelCover` of `weights` added up by `keys` (arrays of one length, keys repeating)."""
    unique, inverse = np.unique(keys, return_inverse=True)
    return LevelCover(unique, np.bincount(inverse, weights=weights, minlength=len(unique)))


def _parents_of(keys, m):
    """The packed key of the level-`3m` block holding each level-`m` block of `keys` (an int64 array): the parent's shell,
    row and wedge from `drill.drill_parent`, worked out once for each distinct `(shell, band)` and applied to every
    wedge of it."""
    shell, band, wedge = unpack_array(keys)
    pairs, inverse = np.unique(keys >> _WEDGE_BITS, return_inverse=True)
    parent_shell = np.empty(len(pairs), dtype=np.int64)
    parent_band = np.empty(len(pairs), dtype=np.int64)
    divide = np.empty(len(pairs), dtype=np.int64)
    for i, pair in enumerate(pairs):
        pair_shell, pair_band, _ = unpack(int(pair) << _WEDGE_BITS)
        parent = drill.drill_parent(drill.DrillBlock(m, pair_shell, pair_band, 0))
        parent_shell[i], parent_band[i] = parent.shell, parent.band
        divide[i] = (drill.drill_wedge_count(m, pair_shell, pair_band)
                     // drill.drill_wedge_count(parent.m, parent.shell, parent.band))
    return pack(parent_shell[inverse], parent_band[inverse], wedge // divide[inverse])


def roll_up(fine):
    """`fine_cover`'s result as `{m: LevelCover}` for every level, each block's cover the sum of the level-3 blocks
    inside it."""
    result = {drill.DRILL_STEP: fine}
    level = fine
    m = drill.DRILL_STEP
    while m != drill.STAGE_TOP:
        level = _summed(_parents_of(level.keys, m), level.totals) if len(level.keys) else level
        m *= drill.DRILL_STEP
        result[m] = level
    return result


def level_rows(x, y, z, radius, codes, m, edge_pc, covers=None):
    """
    The aggregate rows of level `m` for the spheres at `x`, `y`, `z` (parsecs) of `radius` (parsecs) and class
    `codes` (`class_codes`): a `LevelRows`. Blocks that no sphere's centre is in but some sphere reaches get a row
    with a count of zero and the cover only. `covers` is the level's `LevelCover` of `roll_up`
    when the caller has it, else it is worked out here.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    radius = np.asarray(radius, dtype=float)
    codes = np.asarray(codes, dtype=np.int64)
    volume = sphere_volumes(radius)
    if covers is None:
        covers = roll_up(fine_cover(x, y, z, radius, edge_pc))[m]
    classes = []
    unique = np.zeros(0, dtype=np.int64)
    if len(x):
        keys = _centre_keys(x, y, z, m, edge_pc)
        unique, inverse = np.unique(keys, return_inverse=True)
        groups = len(unique)
        count = np.bincount(inverse, minlength=groups)
        total = np.bincount(inverse, weights=volume, minlength=groups)
        centroid = np.column_stack([np.bincount(inverse, weights=volume * axis, minlength=groups) / total
                                    for axis in (x, y, z)])
        distance = np.sqrt((x - centroid[inverse, 0]) ** 2 + (y - centroid[inverse, 1]) ** 2
                           + (z - centroid[inverse, 2]) ** 2) + radius
        bound = np.zeros(groups)
        np.maximum.at(bound, inverse, distance)
        families = np.zeros((groups, len(FAMILIES)))
        np.add.at(families, (inverse, _FAMILY_INDEX_OF_CLASS[codes]), volume)
        classes_n = len(CLASS_LETTERS)
        class_keys, class_inverse = np.unique(inverse * classes_n + codes, return_inverse=True)
        class_volume = np.bincount(class_inverse, weights=volume, minlength=len(class_keys))
        class_group = class_keys // classes_n
        class_code = class_keys % classes_n
        order = np.lexsort((class_volume, class_group))
        last = np.ones(len(order), dtype=bool)
        last[:-1] = class_group[order][1:] != class_group[order][:-1]
        dominant = np.empty(groups, dtype=np.int64)
        dominant[class_group[order][last]] = class_code[order][last]
        for group, code, vol in zip(class_group, class_code, class_volume):
            shell, band, wedge = unpack(unique[group])
            classes.append((shell, band, wedge, CLASS_LETTERS[int(code)], float(vol)))
    every = np.union1d(unique, covers.keys)
    cover_at = _lookup(covers.keys, every)
    cover_volumes = np.where(cover_at >= 0, covers.totals[np.maximum(cover_at, 0)], 0.0)
    group_at = _lookup(unique, every)
    if len(x):
        groups = {"nebula_count": count, "volume_pc3": total, "centroid": centroid, "bound": bound,
                  "families": families, "dominant": dominant}
    else:
        groups = None
    return LevelRows(BlockRows(m, edge_pc, every, cover_volumes, group_at, groups), classes)


class BlockRows(Sequence):
    """
    The block dicts of one level (`LevelRows.blocks`), made as they are read, from the arrays that hold them: a level
    has millions of blocks at the finest sizes, and a dict each is most of a run's memory (PERF.193). Reads like a
    list: `len`, an index, a slice (a list of dicts) and iteration.
    """

    def __init__(self, m, edge_pc, keys, cover_volumes, group_at, groups):
        self._m = m
        self._edge_pc = edge_pc
        self._keys = keys
        self._cover_volumes = cover_volumes
        self._group_at = group_at
        self._groups = groups

    def __len__(self):
        return len(self._keys)

    def __getitem__(self, index):
        if isinstance(index, slice):
            return [self._row(i) for i in range(*index.indices(len(self)))]
        if index < 0:
            index += len(self)
        if not 0 <= index < len(self):
            raise IndexError(index)
        return self._row(index)

    def _row(self, i):
        shell, band, wedge = unpack(int(self._keys[i]))
        g = int(self._group_at[i])
        if g >= 0:
            groups = self._groups
            row = {
                "nebula_count": int(groups["nebula_count"][g]), "volume_pc3": float(groups["volume_pc3"][g]),
                "centroid_x_pc": float(groups["centroid"][g, 0]), "centroid_y_pc": float(groups["centroid"][g, 1]),
                "centroid_z_pc": float(groups["centroid"][g, 2]), "bound_radius_pc": float(groups["bound"][g]),
                "families": [float(v) for v in groups["families"][g]],
                "dominant_class": CLASS_LETTERS[int(groups["dominant"][g])]}
        else:
            row = {
                "nebula_count": 0, "volume_pc3": 0.0, "centroid_x_pc": None, "centroid_y_pc": None,
                "centroid_z_pc": None, "bound_radius_pc": None, "families": [0.0] * len(FAMILIES),
                "dominant_class": None}
        cover_volume = float(self._cover_volumes[i])
        row.update(m=self._m, shell_index=shell, band_index=band, wedge_index=wedge, cover_volume_pc3=cover_volume,
                   cover=cover_of(cover_volume, block_volume_pc3(self._m, shell, band, wedge, self._edge_pc)))
        return row


def _lookup(sorted_keys, keys):
    """For each of `keys` its position in `sorted_keys` (sorted, no repeats), or -1 where it is not there."""
    if not len(sorted_keys):
        return np.full(len(keys), -1, dtype=np.int64)
    at = np.minimum(np.searchsorted(sorted_keys, keys), len(sorted_keys) - 1)
    return np.where(sorted_keys[at] == keys, at, -1)


def iter_rows(x, y, z, radius, codes, edge_pc, levels=LEVELS, on_level=None):
    """
    `level_rows` for every level, one `(m, LevelRows)` at a time, so a caller that writes each level as it comes keeps
    one level's rows in memory instead of every level's (PERF.187); a level's cover is dropped once its rows are made.
    `on_level(done, total)` is called after each.
    """
    covers = roll_up(fine_cover(x, y, z, radius, edge_pc))
    for done, m in enumerate(levels, start=1):
        level = level_rows(x, y, z, radius, codes, m, edge_pc, covers.pop(m))
        yield m, level
        del level
        if on_level is not None:
            on_level(done, len(levels))


def all_rows(x, y, z, radius, codes, edge_pc, levels=LEVELS, on_level=None):
    """`level_rows` for every level: `{m: LevelRows}`. `on_level(done, total)` is called after each."""
    return dict(iter_rows(x, y, z, radius, codes, edge_pc, levels, on_level))
