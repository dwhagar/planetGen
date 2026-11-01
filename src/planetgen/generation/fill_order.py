# planetgen/generation/fill_order.py

"""
The order a neighbourhood or radial fill generates its sectors in (GEN.101)
===========================================================================

Boss (2026-10-07): fill the sectors nearest the original sector first, circling outward, along a geometric
space-filling curve through a pruned Hilbert octree. A pruned ball is not always walkable with face-adjacent steps
(docs/design/galaxy/fill-order-curves-and-core.md: 16 of 1,335 balls pass the parity test), so the order is deterministic and
mostly adjacent, and every step that is not adjacent is counted and reported. The decided variant (the default of the
open question, 2026-10-10) is the 4-voxel-block one, so the first sector is the centre:

* The candidates are laid on a lattice of one sector edge around the centre (`lattice_offset`: the sector centre minus
  the fill's centre, in edges, rounded to the nearest voxel; the grid is spherical, so cells near a voxel share it
  and tie-break on distance and address).
* The lattice is cut into blocks of `BLOCK_EDGES` voxels a side, the centre voxel in the middle of the block at the
  origin. Blocks come nearest first (the distance of the block's middle from the centre), so the ball is filled
  from the centre outward in shells of blocks; blocks the same distance away come in the Hilbert order of their
  coordinates, so the walk goes from one to a neighbour.
* Inside a block the voxels come in the Hilbert order of their place in the block (a curve of two levels whose steps
  are all face-adjacent), then by distance from the centre and by address, so the whole order is a pure function of
  the set of cells and the centre: shuffling the candidates changes nothing. The sector nearest the centre comes
  first of all, whatever its place on the curve.

`fill_order` returns the candidates sorted; `step_statistics` counts the steps longer than `JUMP_EDGES` sector edges,
which `run_galaxy` logs once per run.
"""

import math
from planetgen.util import detmath

BLOCK_EDGES = 4
"""int: Voxels to a block side."""

JUMP_EDGES = 2.5
"""float: A step between consecutive sectors longer than this many sector edges is a jump (a face-adjacent step is one
edge, and a step over one voxel the grid has no cell in is two; the cells are not exactly cubes, so up to two and a half
edges is still a neighbour step)."""


def hilbert_index(point, bits):
    """
    The index along the three-dimensional Hilbert curve of `bits` levels of the lattice point `point` (three integers
    in `0 .. 2**bits - 1`), by Skilling's transpose algorithm ([R4113], "Programming the Hilbert curve"). Points
    one apart on the curve are face-adjacent.
    """
    x = [int(c) for c in point]
    if not all(0 <= c < 1 << bits for c in x):
        raise ValueError(f"{point} is outside a lattice of {bits} levels")
    top = 1 << (bits - 1)
    q = top
    while q > 1:
        p = q - 1
        for i in range(3):
            if x[i] & q:
                x[0] ^= p
            else:
                swap = (x[0] ^ x[i]) & p
                x[0] ^= swap
                x[i] ^= swap
        q >>= 1
    for i in range(1, 3):
        x[i] ^= x[i - 1]
    flip = 0
    q = top
    while q > 1:
        if x[2] & q:
            flip ^= q - 1
        q >>= 1
    for i in range(3):
        x[i] ^= flip
    index = 0
    for bit in range(bits - 1, -1, -1):
        for i in range(3):
            index = (index << 1) | ((x[i] >> bit) & 1)
    return index


def lattice_offset(position, center, edge_pc):
    """The lattice voxel of `position`: its offset from `center`, in sector edges, to the nearest integers."""
    return tuple(math.floor((p - c) / edge_pc + 0.5) for p, c in zip(position, center))


def _block_of(voxel):
    """`(block, local)` of a voxel: the block coordinates (the centre voxel is in the middle of block 0) and the
    place inside it, each in `0 .. BLOCK_EDGES - 1`."""
    half = BLOCK_EDGES // 2
    block = tuple((c + half) // BLOCK_EDGES for c in voxel)
    local = tuple(c + half - b * BLOCK_EDGES for c, b in zip(voxel, block))
    return block, local


def fill_order(candidates, center, edge_pc):
    """
    `candidates` in fill order.

    Args:
        candidates (iterable): `(shell, band, wedge, x, y, z, distance_pc)` as
            `spherical.enumerate_sectors_within_radius` yields them.
        center (tuple): The fill's centre, parsecs.
        edge_pc (float): The sector edge.

    Returns:
        list: The same tuples, the centre's sector first when it is among them.
    """
    cells = list(candidates)
    if not cells:
        return cells
    first = min(cells, key=lambda cell: (round(cell[6], 9), cell[:3]))[:3]
    placed = []
    for cell in cells:
        block, local = _block_of(lattice_offset(cell[3:6], center, edge_pc))
        placed.append((cell, block, local))
    reach = max(abs(c) for _cell, block, _local in placed for c in block)
    bits = max(1, (2 * reach).bit_length())
    local_bits = (BLOCK_EDGES - 1).bit_length()
    half = BLOCK_EDGES / 2.0

    def key(entry):
        cell, block, local = entry
        middle = sum((BLOCK_EDGES * c - half + 0.5 * (BLOCK_EDGES - 1)) ** 2 for c in block)
        return (cell[:3] != first, round(middle, 9), hilbert_index([c + reach for c in block], bits),
                hilbert_index(local, local_bits), round(cell[6], 9), cell[:3])

    return [cell for cell, _block, _local in sorted(placed, key=key)]


def step_statistics(ordered, edge_pc):
    """
    How well the order walks: `(steps, jumps, mean_jump_edges, longest_jump_edges)` over the steps between consecutive
    sectors, a jump being a step longer than `JUMP_EDGES` edges. `mean_jump_edges` and `longest_jump_edges` are 0 with
    no jump.
    """
    steps = max(len(ordered) - 1, 0)
    lengths = []
    for before, after in zip(ordered, ordered[1:]):
        length = detmath.dist(before[3:6], after[3:6]) / edge_pc
        if length > JUMP_EDGES:
            lengths.append(length)
    return steps, len(lengths), (sum(lengths) / len(lengths) if lengths else 0.0), max(lengths, default=0.0)
