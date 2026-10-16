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

A position outside the box has a NULL key and is found by the query's `cell_key IS NULL` branch.
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
    return int(math.floor((coordinate_pc + REACH_PC) / CELL_PC))


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
    axes = [_spread_sql(f"CAST(FLOOR((({c}) + {reach}) / {CELL_PC}) AS UNSIGNED)") for c in (x, y, z)]
    inside = " AND ".join(f"ABS({c}) < {reach}" for c in (x, y, z))
    return f"IF({inside}, ({axes[0]} << 2) | ({axes[1]} << 1) | {axes[2]}, NULL)"
