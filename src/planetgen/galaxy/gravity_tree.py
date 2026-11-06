# planetgen/galaxy/gravity_tree.py

"""
The gravity map's cell hierarchy and interaction lists (MAP.385; docs/design/orbits/gravity-map.md, `reports/perf72-gravity-grid-cost-2026-11-05.md`)
==========================================================================================================================================================

The gravity map sums, for every zone of a sector, the point masses around it. Summing every star of the galaxy is out of the question (15 to 17 s per
sector grid against 1.6 million stars), so the far sources are summed as one point mass per cell, the way a tree code does. This module holds the two
pure steps in front of `gravity_kernel`:

1. **`build(positions, masses)`** puts the sources into a hierarchy of cubic cells whose smallest ("leaf") cell is `LEAF_PC` (0.5 pc, PERF.72) wide and whose
   every level doubles the edge. A cell holds its total mass and its center of mass. The cells are the Morton (bit-interleaved) cells of the integer leaf
   coordinates, so a cell's parent is a shift of its code and the whole structure is a few sorted arrays, built in one pass.
2. **`interaction_list(hierarchy, center, radius, theta)`** takes a target (a ball round the zones to be filled: a sector's bounding sphere) and returns the
   point masses to hand the kernel: the stars of every leaf cell the target is too near to take as a point, and the center of mass of every cell that is
   far enough. A cell is taken as a point when `0.866 * edge / (|center of mass - target center| - radius) < theta`, i.e. its half-diagonal over its
   distance from the target ball's surface. `theta` is at most 1 (a cell then cannot reach into the ball); the PERF.72 errors of the pull against the exact
   sum are 0.15 percent at 0.5, 0.45 percent at 0.7 and 4.6 percent at 1.0, so the choice is capped at `MAX_THETA` (0.7).
3. **`choose_theta`** takes the smallest `theta` from `MIN_THETA` (0.3) to `MAX_THETA` (0.7) whose list fits a pair budget (`PAIR_BUDGET`, 50 million
   point-source pairs, about 0.1 s of the kernel on four cores), so a sparse place is summed more exactly than a dense one at the same cost; over budget even
   at the cap it takes the cap and says so.

`fields_at` runs the three and the kernel for one target and returns the ten numbers per point with what it did. The tree is not the sector grid's
storage and not a cache: it is derived from the sources of one run, built in memory, and thrown away.

Units: positions in parsecs, masses in solar masses. The module imports numba (the traversal is compiled), so the web app must not import it (OPS.42).
"""

import dataclasses
import math

import numpy as np
from numba import njit

from planetgen.galaxy import gravity_kernel
from planetgen.galaxy.gravity_settings import LEAF_PC, MAX_THETA, MIN_THETA, PAIR_BUDGET, THETA_STEP

HALF_DIAGONAL = math.sqrt(3.0) / 2.0
"""float: A cube's half-diagonal over its edge."""
_BITS = 21
"""int: Bits of a leaf coordinate in a Morton code (a galaxy of 2 million leaves a side, 1 Mpc at 0.5 pc)."""


@dataclasses.dataclass(frozen=True)
class Hierarchy:
    """
    The cells of one set of sources.

    Attributes:
        origin (numpy.ndarray): The corner the leaf coordinates count from, parsecs.
        leaf_pc (float): The leaf edge, parsecs.
        positions (numpy.ndarray): The sources sorted by leaf cell, `(n, 3)`.
        masses (numpy.ndarray): Their masses, solar masses.
        level_start (numpy.ndarray): Where level `l`'s cells begin in the cell arrays; level 0 holds the leaves, the last level the top; `len == levels + 1`.
        codes (numpy.ndarray): Each cell's Morton code at its own level.
        cell_mass (numpy.ndarray): Each cell's mass.
        cell_center (numpy.ndarray): Each cell's center of mass, `(cells, 3)`.
        child_lo (numpy.ndarray): For a cell above level 0, the first of its children in the cell arrays; for a leaf, the first of its sources in `positions`.
        child_hi (numpy.ndarray): One past the last of them.
    """

    origin: np.ndarray
    leaf_pc: float
    positions: np.ndarray
    masses: np.ndarray
    level_start: np.ndarray
    codes: np.ndarray
    cell_mass: np.ndarray
    cell_center: np.ndarray
    child_lo: np.ndarray
    child_hi: np.ndarray

    @property
    def levels(self):
        """int: The number of levels."""
        return len(self.level_start) - 1

    @property
    def source_count(self):
        """int: The number of sources."""
        return len(self.masses)


def _spread(values):
    """The bits of the integers `values` (below 2^21) spread to every third bit (the Morton interleave of one axis)."""
    v = values.astype(np.uint64) & np.uint64(0x1FFFFF)
    v = (v | (v << np.uint64(32))) & np.uint64(0x1F00000000FFFF)
    v = (v | (v << np.uint64(16))) & np.uint64(0x1F0000FF0000FF)
    v = (v | (v << np.uint64(8))) & np.uint64(0x100F00F00F00F00F)
    v = (v | (v << np.uint64(4))) & np.uint64(0x10C30C30C30C30C3)
    v = (v | (v << np.uint64(2))) & np.uint64(0x1249249249249249)
    return v


def build(positions, masses_solar, leaf_pc=LEAF_PC):
    """
    The `Hierarchy` of the point masses `positions` (parsecs, `(n, 3)`) of `masses_solar`.

    Raises:
        ValueError: If the shapes are wrong, a position or mass is not finite, a mass is negative, the leaf is not positive, or the sources span more than 2^21 leaves
            along an axis.
    """
    positions = np.ascontiguousarray(positions, dtype=np.float64)
    masses = np.ascontiguousarray(masses_solar, dtype=np.float64)
    if positions.ndim != 2 or positions.shape[1] != 3 or masses.shape != (positions.shape[0],):
        raise ValueError("positions is an (n, 3) array and masses_solar has one mass per position")
    if not (leaf_pc > 0.0 and math.isfinite(leaf_pc)):
        raise ValueError("the leaf edge must be positive and finite")
    if not (np.isfinite(positions).all() and np.isfinite(masses).all() and (masses >= 0.0).all()):
        raise ValueError("every position must be finite and every mass zero or more and finite")
    count = len(masses)
    if count == 0:
        empty = np.zeros(0, dtype=np.int64)
        return Hierarchy(np.zeros(3), leaf_pc, positions, masses, np.zeros(1, dtype=np.int64), np.zeros(0, dtype=np.uint64), np.zeros(0), np.zeros((0, 3)), empty, empty)
    origin = positions.min(axis=0)
    leaf = np.floor((positions - origin) / leaf_pc).astype(np.int64)
    if leaf.max() >= (1 << _BITS):
        raise ValueError("the sources span more than 2^21 leaves along an axis")
    code = _spread(leaf[:, 0]) | (_spread(leaf[:, 1]) << np.uint64(1)) | (_spread(leaf[:, 2]) << np.uint64(2))
    order = np.argsort(code, kind="stable")
    code, positions, masses = code[order], positions[order], masses[order]

    codes, mass, center, lo, hi = [], [], [], [], []
    unique, start = np.unique(code, return_index=True)
    mass0 = np.add.reduceat(masses, start)
    moment = np.stack([np.add.reduceat(masses * positions[:, k], start) for k in range(3)], 1)
    plain = np.stack([np.add.reduceat(positions[:, k], start) for k in range(3)], 1) / np.diff(np.append(start, count))[:, None]
    codes.append(unique)
    mass.append(mass0)
    center.append(np.where(mass0[:, None] > 0.0, moment / np.where(mass0 > 0.0, mass0, 1.0)[:, None], plain))
    lo.append(start.astype(np.int64))
    hi.append(np.append(start[1:], count).astype(np.int64))
    sizes = [len(unique)]
    while sizes[-1] > 1:
        below = codes[-1] >> np.uint64(3)
        unique, first = np.unique(below, return_index=True)
        total = np.add.reduceat(mass[-1], first)
        moment = np.stack([np.add.reduceat(mass[-1] * center[-1][:, k], first) for k in range(3)], 1)
        plain = np.stack([np.add.reduceat(center[-1][:, k], first) for k in range(3)], 1) / np.diff(np.append(first, len(below)))[:, None]
        codes.append(unique)
        mass.append(total)
        center.append(np.where(total[:, None] > 0.0, moment / np.where(total > 0.0, total, 1.0)[:, None], plain))
        base = int(np.sum(sizes[:-1]))
        lo.append(first.astype(np.int64) + base)
        hi.append(np.append(first[1:], len(below)).astype(np.int64) + base)
        sizes.append(len(unique))
    return Hierarchy(origin, leaf_pc, positions, masses, np.concatenate([[0], np.cumsum(sizes)]).astype(np.int64), np.concatenate(codes), np.concatenate(mass),
                     np.concatenate(center), np.concatenate(lo), np.concatenate(hi))


@njit(cache=True)
def _walk(level_start, cell_mass, cell_center, child_lo, child_hi, leaf_pc, positions, masses, center, radius, theta, exact_total, fill, out_positions, out_masses):
    """
    Walks the cells from the top level down for one target ball. Returns `(exact, points)`: the number of sources in the leaves too near to take as points, and the
    number of cells taken as points. When `fill`, writes the exact sources to the outputs from index 0 and the points from index `exact_total`.
    """
    levels = len(level_start) - 1
    top = levels - 1
    capacity = (level_start[top + 1] - level_start[top]) + 8 * levels + 8
    stack = np.empty(capacity, dtype=np.int64)
    stack_level = np.empty(capacity, dtype=np.int64)
    depth = 0
    for cell in range(level_start[top], level_start[top + 1]):
        stack[depth] = cell
        stack_level[depth] = top
        depth += 1
    exact = 0
    points = 0
    while depth > 0:
        depth -= 1
        cell = stack[depth]
        level = stack_level[depth]
        if cell_mass[cell] <= 0.0:
            continue
        dx = cell_center[cell, 0] - center[0]
        dy = cell_center[cell, 1] - center[1]
        dz = cell_center[cell, 2] - center[2]
        distance = math.sqrt(dx * dx + dy * dy + dz * dz) - radius
        if distance > 0.0 and 0.8660254037844386 * leaf_pc * (2.0 ** level) < theta * distance:
            if fill:
                out_positions[exact_total + points, 0] = cell_center[cell, 0]
                out_positions[exact_total + points, 1] = cell_center[cell, 1]
                out_positions[exact_total + points, 2] = cell_center[cell, 2]
                out_masses[exact_total + points] = cell_mass[cell]
            points += 1
        elif level == 0:
            for source in range(child_lo[cell], child_hi[cell]):
                if fill:
                    out_positions[exact, 0] = positions[source, 0]
                    out_positions[exact, 1] = positions[source, 1]
                    out_positions[exact, 2] = positions[source, 2]
                    out_masses[exact] = masses[source]
                exact += 1
        else:
            for child in range(child_lo[cell], child_hi[cell]):
                stack[depth] = child
                stack_level[depth] = level - 1
                depth += 1
    return exact, points


@dataclasses.dataclass(frozen=True)
class InteractionList:
    """
    What the kernel is handed for one target.

    Attributes:
        positions (numpy.ndarray): The point masses, parsecs, `(exact + points, 3)`: the exact sources first, then the cells' centers of mass.
        masses (numpy.ndarray): Their masses, solar masses.
        exact (int): How many of them are single sources.
        points (int): How many are cells.
        theta (float): The opening angle that made the list.
    """

    positions: np.ndarray
    masses: np.ndarray
    exact: int
    points: int
    theta: float

    @property
    def size(self):
        """int: The number of point masses."""
        return self.exact + self.points


def _check_target(center, radius, theta):
    center = np.ascontiguousarray(center, dtype=np.float64)
    if center.shape != (3,) or not np.isfinite(center).all():
        raise ValueError("the target center is three finite numbers")
    if not (math.isfinite(radius) and radius >= 0.0):
        raise ValueError("the target radius must be zero or more and finite")
    if theta is not None and not (0.0 < theta <= 1.0):
        raise ValueError("the opening angle must be above 0 and at most 1")
    return center


def count(hierarchy, center, radius, theta):
    """`(exact, points)`: the size of the list `interaction_list` would make, without making it."""
    center = _check_target(center, radius, theta)
    if hierarchy.source_count == 0:
        return 0, 0
    nothing = np.zeros((0, 3))
    return tuple(int(n) for n in _walk(hierarchy.level_start, hierarchy.cell_mass, hierarchy.cell_center, hierarchy.child_lo, hierarchy.child_hi, hierarchy.leaf_pc,
                                       hierarchy.positions, hierarchy.masses, center, float(radius), float(theta), 0, False, nothing, nothing[:, 0]))


def interaction_list(hierarchy, center, radius, theta):
    """
    The `InteractionList` for the target ball of `radius` parsecs round `center`, at opening angle `theta` (above 0, at most 1). The same hierarchy, target and
    angle give the same list in the same order.
    """
    exact, points = count(hierarchy, center, radius, theta)
    positions = np.empty((exact + points, 3))
    masses = np.empty(exact + points)
    if exact + points:
        _walk(hierarchy.level_start, hierarchy.cell_mass, hierarchy.cell_center, hierarchy.child_lo, hierarchy.child_hi, hierarchy.leaf_pc, hierarchy.positions,
              hierarchy.masses, np.ascontiguousarray(center, dtype=np.float64), float(radius), float(theta), exact, True, positions, masses)
    return InteractionList(positions, masses, exact, points, float(theta))


def choose_theta(hierarchy, center, radius, point_count, budget=PAIR_BUDGET, minimum=MIN_THETA, maximum=MAX_THETA):
    """
    The smallest opening angle from `minimum` to `maximum` (steps of `THETA_STEP`) whose list costs no more than `budget` pairs with `point_count` points; the
    `maximum` when none does. Returns `(theta, fits)`.
    """
    steps = int(round((maximum - minimum) / THETA_STEP))
    for step in range(steps + 1):
        theta = round(minimum + step * THETA_STEP, 10)
        exact, points = count(hierarchy, center, radius, theta)
        if (exact + points) * point_count <= budget:
            return theta, True
    return round(maximum, 10), False


def bounding_sphere(points):
    """`(center, radius)` of a ball round `points`, parsecs: the middle of their box and the distance to the farthest."""
    points = np.asarray(points, dtype=np.float64)
    center = (points.min(axis=0) + points.max(axis=0)) / 2.0
    return center, float(np.sqrt(((points - center) ** 2).sum(axis=1)).max())


@dataclasses.dataclass(frozen=True)
class Result:
    """What `fields_at` returns: the ten numbers per point (`gravity_kernel`'s columns) and what made them (`list`, with `fits` the budget answer)."""

    values: np.ndarray
    list: InteractionList
    fits: bool


def fields_at(hierarchy, points, theta=None, budget=PAIR_BUDGET, softening_pc=gravity_kernel.DEFAULT_SOFTENING_PC, threads=None):
    """
    The field at `points` (one target: the ball round them) from the hierarchy's sources. `theta` `None` takes `choose_theta` against `budget`.
    """
    points = np.ascontiguousarray(points, dtype=np.float64)
    center, radius = bounding_sphere(points)
    fits = True
    if theta is None:
        theta, fits = choose_theta(hierarchy, center, radius, len(points), budget)
    chosen = interaction_list(hierarchy, center, radius, theta)
    return Result(gravity_kernel.fields(points, chosen.positions, chosen.masses, softening_pc, threads), chosen, fits)
