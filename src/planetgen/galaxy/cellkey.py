# planetgen/galaxy/cellkey.py

"""
The cell key: where an object is *now*, as one integer (GEN.223, region queries).

The galaxy is cut into cubes of `CELL_PC` (0.25 pc) over a box of plus or minus `REACH_PC`
(262,144 pc, the grid's 262 kpc reach with room to spare), and the three cube numbers (21 bits each) are woven into one
63-bit number, x in the highest bit of each triple (a Morton, or Z-order, key). Cubes that sit inside one bigger aligned
cube have keys that form one unbroken run, so a sphere becomes a short list of key runs
(docs/design/spherical-region-queries.md).

It is not the object ID. The ID is the permanent birth address and never changes; the cell key is derived from the
position and changes when the object moves into another cube. MySQL computes it itself, as a stored generated column
(`sql_cell_key`) on every table that holds a position, so it cannot disagree with the position columns, and a covering
index `(cell_key, x, y, z)` answers a sphere without opening a row. This module holds the same function in Python
(`cell_key`, `cell_keys`) so the two can be tested against each other.

A position outside the box has a NULL key and is found by the query's `cell_key IS NULL` branch. A position just inside
the box edge whose rounded cube number would be `2 ** BITS` (the last double below `REACH_PC`) is put in the last cube,
`LAST_CELL`, so the key never wraps to the opposite edge (DB.65).
"""

import math

import numpy as np

REACH_PC = 262144.0
"""float: Half the width of the box the key covers, parsecs (2 ** 18)."""
CELL_PC = 0.25
"""float: The edge of one cube, parsecs."""
BITS = 21
"""int: Bits per axis; `2 * REACH_PC / CELL_PC` is `2 ** BITS`."""
KEY_BITS = 3 * BITS
LAST_CELL = (1 << BITS) - 1
"""int: The highest cube number; a coordinate whose floored number rounds up past it is held here."""

SPREAD_STEPS = (
    (32, 0x1F00000000FFFF),
    (16, 0x1F0000FF0000FF),
    (8, 0x100F00F00F00F00F),
    (4, 0x10C30C30C30C30C3),
    (2, 0x1249249249249249),
)
"""tuple: `(shift, mask)` of the five steps that move bit i of a 21-bit number to bit 3i."""


def spread(value):
    """Moves bit i of `value` (a 21-bit integer) to bit 3i, the way the key interleaves an axis."""
    for shift, mask in SPREAD_STEPS:
        value = (value | (value << shift)) & mask
    return value


def cell_index(coordinate_pc):
    """The cube number of one coordinate (parsecs), 0 to `2 ** BITS - 1`, or None outside the box."""
    if not abs(coordinate_pc) < REACH_PC:
        return None
    return min(int(math.floor((coordinate_pc + REACH_PC) / CELL_PC)), LAST_CELL)


def cell_key(x_pc, y_pc, z_pc):
    """The key of a point (parsecs, galaxy frame), or None when it lies outside the box."""
    cells = [cell_index(c) for c in (x_pc, y_pc, z_pc)]
    if None in cells:
        return None
    return (spread(cells[0]) << 2) | (spread(cells[1]) << 1) | spread(cells[2])


def cell_keys(points):
    """
    The keys of many points at once.

    Args:
        points (array-like): Shape `(n, 3)`, parsecs, galaxy frame.

    Returns:
        tuple: `(keys, valid)`: a `uint64` array and a boolean array; `keys` is 0 where `valid` is False (outside the box).
    """
    points = np.asarray(points, dtype=np.float64).reshape(-1, 3)
    valid = np.all(np.abs(points) < REACH_PC, axis=1)
    cells = np.floor((np.where(valid[:, None], points, 0.0) + REACH_PC) / CELL_PC).astype(np.uint64)
    cells = np.minimum(cells, np.uint64(LAST_CELL))
    keys = np.zeros(len(points), dtype=np.uint64)
    for axis, lift in enumerate((2, 1, 0)):
        value = cells[:, axis]
        for shift, mask in SPREAD_STEPS:
            value = (value | (value << np.uint64(shift))) & np.uint64(mask)
        keys |= value << np.uint64(lift)
    return np.where(valid, keys, np.uint64(0)), valid


def _chunk_steps():
    # How a 7-bit chunk is spread: bit i goes to 3i in three steps that move the bits whose index has the 4, then the 2,
    # then the 1 set (the classic Morton chain), each as a multiply by `2 ** shift + 1` and a mask of the positions the
    # chunk's bits hold after it.
    steps = []
    for weight in (4, 2, 1):
        mask = 0
        for i in range(7):
            mask |= 1 << (i + 2 * (i & ~(weight - 1)))
        steps.append((2 * weight, mask))
    return tuple(steps)


CHUNK_BITS = 7
CHUNK_STEPS = _chunk_steps()
"""tuple: `(shift, mask)` of the three steps that spread a 7-bit chunk of a cube number (see `sql_cell_key`)."""


def _spread_sql(cell):
    # MySQL has no variables in a generated column, and `(v | v << s)` names `v` twice, which would double the
    # expression at each of five steps. A multiply by `2 ** s + 1` names it once, but overflows a 64-bit product when the
    # top bits of a 21-bit number move, and MySQL raises an error on overflow instead of dropping the bits. So the cube
    # number is cut into three 7-bit chunks (21 bits in all), each spread with three multiplies that cannot overflow, and
    # the chunks are joined.
    parts = []
    for index in range(3):
        chunk = f"((({cell}) >> {index * CHUNK_BITS}) & {(1 << CHUNK_BITS) - 1})"
        for shift, mask in CHUNK_STEPS:
            chunk = f"(({chunk} * {(1 << shift) + 1}) & {mask})"
        parts.append(f"({chunk} << {index * CHUNK_BITS * 3})")
    return "(" + " | ".join(parts) + ")"


def sql_cell_key(x, y, z):
    """
    The MySQL expression for the key of a point whose coordinates (parsecs, galaxy frame) are the SQL expressions `x`,
    `y` and `z`, NULL outside the box. It is what a `GENERATED ALWAYS AS (...) STORED` column holds.
    """
    reach = f"{int(REACH_PC)}"
    axes = [_spread_sql(f"LEAST(CAST(FLOOR((({c}) + {reach}) / {CELL_PC}) AS UNSIGNED), {LAST_CELL})") for c in (x, y, z)]
    inside = " AND ".join(f"ABS({c}) < {reach}" for c in (x, y, z))
    return f"IF({inside}, ({axes[0]} << 2) | ({axes[1]} << 1) | {axes[2]}, NULL)"


# --- The sphere as runs of keys (GEN.276) ------------------------------------------------------------------

BUDGET = 256
"""int: The most runs `sphere_runs` returns by default (the flat part of the curve, design note 6.3)."""

_INFLATE_CELLS = 1e-6
"""float: The sphere is laid over the grid a millionth of a cell fatter than it is, so a point on its surface that
float rounding put on the wrong side of a cube's face is still inside a run."""


def _spread_array(values):
    values = values.astype(np.uint64)
    for shift, mask in SPREAD_STEPS:
        values = (values | (values << np.uint64(shift))) & np.uint64(mask)
    return values


def _weave(ix, iy, iz):
    return (_spread_array(ix) << np.uint64(2)) | (_spread_array(iy) << np.uint64(1)) | _spread_array(iz)


def _merged(lows, highs):
    """Sorts the runs `[lows[i], highs[i]]` and joins the ones that overlap or touch."""
    if len(lows) == 0:
        return lows, highs
    order = np.argsort(lows, kind="stable")
    lows, highs = lows[order], highs[order]
    reach = np.maximum.accumulate(highs)
    starts = np.ones(len(lows), dtype=bool)
    starts[1:] = lows[1:] > reach[:-1] + np.uint64(1)
    first = np.flatnonzero(starts)
    last = np.append(first[1:] - 1, len(lows) - 1)
    return lows[first], reach[last]


def sphere_runs(center, r, budget=BUDGET, reach_pc=REACH_PC):
    """
    The cell keys a sphere can hold, as at most `budget` ascending, merged runs `[(low, high), ...]` (inclusive):
    every point within `r` of `center` whose cell has a key lies in one of them. Used as the `VALUES` list of the
    region query (design note 6.3, 6.4); the runs hold some keys of cells outside the sphere, and the exact distance
    test is the query's own. Integers and numpy only: the web app never imports numba.

    The sphere is laid over the nested grid of cubes the key weaves: the cover starts at the coarsest level whose
    cubes are at least four radii wide, then each level splits the cubes the sphere's surface crosses (the cubes
    wholly inside it become runs and stay) for as long as the merged run count stays within `budget`, so the
    finest cubes are about a quarter to a sixth of the radius. A sphere that reaches past the box gets runs for
    the part inside it only (a point outside has no key, `cell_key`); one wholly outside gets none.

    Args:
        center (sequence): The sphere's center `(x, y, z)`, parsecs, galaxy frame.
        r (float): Its radius, parsecs, zero or more.
        budget (int): The most runs to return, 1 or more.
        reach_pc (float): Half the width of the box, when the plan stores another reach than `REACH_PC`; the cell
            is then `2 x reach / 2 ** BITS`, as for the key.

    Returns:
        list[tuple[int, int]]: The runs, ascending and not touching one another; `[]` for a sphere wholly outside
        the box.

    Raises:
        ValueError: A radius, reach or budget that is not a finite positive number (radius: not negative), or a
            center that is not three finite numbers.
    """
    center = [float(c) for c in center]
    if len(center) != 3 or not all(math.isfinite(c) for c in center):
        raise ValueError("the center is three finite numbers")
    r, reach_pc = float(r), float(reach_pc)
    if not (math.isfinite(r) and r >= 0.0) or not (math.isfinite(reach_pc) and reach_pc > 0.0):
        raise ValueError("the radius is zero or more and the reach more than zero, both finite")
    if int(budget) < 1:
        raise ValueError("the budget is at least one run")
    cell = 2.0 * reach_pc / (1 << BITS)
    where = [(c + reach_pc) / cell for c in center]
    radius = r / cell + _INFLATE_CELLS
    size = 1 << BITS
    # The cubes of the coarsest level whose edge is at least four radii, clipped to the box.
    level = BITS
    while level > 0 and (1 << (BITS - level)) < 4.0 * radius:
        level -= 1
    edge = 1 << (BITS - level)
    axes = []
    for c in where:
        first = max(int(math.floor((c - radius) / edge)), 0)
        last = min(int(math.floor((c + radius) / edge)), (1 << level) - 1)
        if last < first:
            return []
        axes.append(np.arange(first, last + 1, dtype=np.int64))
    nodes = [a.ravel() for a in np.meshgrid(*axes, indexing="ij")]
    kept_lows, kept_highs = [], []      # the cubes wholly inside, from the levels already split
    best = None                         # the last cover that fit the budget: (lows, highs)
    while True:
        edge = 1 << (BITS - level)
        near = np.zeros(len(nodes[0]))
        far = np.zeros(len(nodes[0]))
        for axis in range(3):
            low = nodes[axis] * float(edge)
            high = low + edge
            gap = np.maximum(np.maximum(low - where[axis], 0.0), where[axis] - high)
            near += gap * gap
            span = np.maximum(np.abs(where[axis] - low), np.abs(where[axis] - high))
            far += span * span
        touching = near <= radius * radius
        inside = touching & (far <= radius * radius)
        crossing = touching & ~inside
        span_mask = np.uint64((1 << (3 * (BITS - level))) - 1)

        def runs_of(chosen):
            if not chosen.any():
                return np.empty(0, dtype=np.uint64), np.empty(0, dtype=np.uint64)
            lows = _weave(*(nodes[axis][chosen] * edge for axis in range(3)))
            return lows, lows | span_mask

        inside_lows, inside_highs = runs_of(inside)
        crossing_lows, crossing_highs = runs_of(crossing)
        lows, highs = _merged(np.concatenate(kept_lows + [inside_lows, crossing_lows]),
                              np.concatenate(kept_highs + [inside_highs, crossing_highs]))
        if len(lows) > budget:
            if best is None:                # even the first cover is too fine: take a coarser one
                if level == 0:
                    break
                level -= 1
                nodes = [axis_nodes // 2 for axis_nodes in nodes]
                nodes = [np.unique(np.stack(nodes), axis=1)[axis] for axis in range(3)]
                kept_lows, kept_highs = [], []
                continue
            break
        best = (lows, highs)
        if not crossing.any() or level == BITS:
            break
        kept_lows.append(inside_lows)
        kept_highs.append(inside_highs)
        parent = [axis_nodes[crossing] for axis_nodes in nodes]
        offsets = [(a, b, c) for a in (0, 1) for b in (0, 1) for c in (0, 1)]
        nodes = [np.concatenate([parent[axis] * 2 + offset[axis] for offset in offsets]) for axis in range(3)]
        level += 1
    if best is None:                        # a budget of one run, and the first cover had more: the whole box
        return [(0, (1 << KEY_BITS) - 1)]
    return [(int(low), int(high)) for low, high in zip(*best)]
