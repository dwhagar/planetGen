"""
planetgen/galaxy/drill.py

The Galaxy Map's block ladder on the spherical sector grid
(docs/design/galaxy/spherical-block-ladder.md; the grid is `galaxy/spherical.py`).
Blocks are `m` shells deep, `m` = 177,147 down to 3 (every size three times its
child), and a block of one is a sector; each block sits wholly inside one block
of the next size up. The top, 177,147 shells (708 kpc, beyond the grid's 262 kpc
reach) is a root of three blocks: its one block shell is one row, cut into three
wedges, and the galaxy lies wholly inside it (MAP.258: the nodes of the
single-space map are these blocks, and the cut walks them from the root).
`static/galaxyprisms.js`'s `drill*` functions are the same rules for the page;
`tests/test_galaxydrill.py` checks that the two agree.

A block is `(m, shell, band, wedge)`:

- `shell`: block shell `s`, holding the sector shells `s*m` to `s*m + m - 1`
  (block shell 0 also holds the core ball). Its centreline radius is
  `R = s*m + (m - 1)/2` edges.
- `band`: a row of blocks between two colatitudes, signed like a sector's band
  (0 the equatorial row, positive north). A block shell has `drill_band_count`
  rows, an odd number of equal colatitude steps, so the galactic plane runs
  through the middle of the equatorial row. A row's blocks aim at `m` edges of
  arc.
- `wedge`: one of the row's `drill_wedge_count` equal longitude steps, from
  +X counterclockwise.

At `m = 1` the block is a sector, `(1, shell, band, wedge)` as
`spherical.py` names it.

Blocks nest by counts. A child block shell's row count is a whole odd
multiple of its parent's, so every parent row edge is a child row edge, and a
child row's wedge count is a whole multiple of the wedge count of the parent
row that holds it, so every parent wedge line is a child wedge line. A sector
joins the level-3 block holding its centre (Boss, 2026-10-15 06:58Z): the
block shell of its shell, the row holding its band's centre colatitude, and
the wedge holding its centre longitude. The core and the polar caps lie on the
axis: they join wedge 0 of the polar rows (the core, the equatorial row).

Everything is integer math apart from the row and wedge counts, which are
rounded from floats; `galaxyprisms.js` copies each expression step for step
(including JavaScript's round-half-up) so both sides pick the same counts.
"""

import math
from bisect import bisect_left
from collections import namedtuple
from functools import lru_cache

from planetgen.galaxy import spherical

DRILL_LEVELS = (177147, 59049, 19683, 6561, 2187, 729, 243, 81, 27, 9, 3, 1)
"""tuple: Block sizes from the root down to a sector. Each step divides by
`DRILL_STEP`."""

DRILL_STEP = 3
"""int: How many child block shells deep one block is."""

DRILL_TOP = DRILL_LEVELS[0]
"""int: The root's size: the first level whose block shell 0 reaches past the grid's outermost shell (65,535,
`spherical.MAX_ADDRESSABLE_SHELL`), so the whole galaxy sits in its three blocks."""

STAGE_LEVELS = (243, 81, 27, 9, 3, 1)
"""tuple: The levels the drill-down stages and map windows show (`galaxy_stage`, `window`): the ladder from the
level-243 blocks of the galaxy's own view down to a sector. The levels above are only for the map's nodes."""

STAGE_TOP = STAGE_LEVELS[0]
"""int: The size of the blocks the galaxy stage lists."""

DrillBlock = namedtuple("DrillBlock", ["m", "shell", "band", "wedge"])
"""One block of the ladder; at `m = 1`, a sector `(1, shell, band, wedge)`."""


def _js_round(x):
    """JavaScript's `Math.round`: halves round up, not to even."""
    return math.floor(x + 0.5)


def _odd_round(x):
    """The odd integer nearest `x`, a tie going up, and at least 1."""
    return max(1, 2 * math.floor(x / 2) + 1)


def _check_level(m):
    if m not in DRILL_LEVELS[:-1]:
        raise ValueError(f"block size must be one of {DRILL_LEVELS[:-1]}, got {m}")


def check_stage_level(m):
    """`m` when a stage or map window may show blocks of that size (`STAGE_LEVELS`, a sector included)."""
    if m not in STAGE_LEVELS:
        raise ValueError(f"stage block size must be one of {STAGE_LEVELS}, got {m}")
    return m


def block_radius(m, shell):
    """Centreline radius of block shell `shell` at size `m`, in edges."""
    return shell * m + (m - 1) / 2


@lru_cache(maxsize=65536)
def drill_band_count(m, shell):
    """
    Rows of blocks in block shell `shell` at size `m` (always odd):

    - the root (`DRILL_TOP`): the odd integer nearest `pi * R / m`, `R` the
      block shell's centreline radius in edges, so a row's arc is about `m`
      edges.
    - every other level: the parent block shell's (`shell // 3`) count times
      the odd integer nearest `T / parent`, `T` as above.

    Raises:
        ValueError: For a level that is not a block size or a negative shell.
    """
    _check_level(m)
    if shell < 0:
        raise ValueError(f"shell must be >= 0, got {shell}")
    target = math.pi * block_radius(m, shell) / m
    if m == DRILL_TOP:
        return _odd_round(target)
    parent = drill_band_count(m * DRILL_STEP, shell // DRILL_STEP)
    return parent * _odd_round(target / parent)


def _row_of(m, shell, band):
    """The north-to-south index of a signed row."""
    count = drill_band_count(m, shell)
    row = (count - 1) // 2 - band
    if not 0 <= row < count:
        raise ValueError(
            f"band {band} is not a row of block shell {shell} at size {m} (rows {-(count // 2)}..{count // 2})")
    return row


@lru_cache(maxsize=1 << 20)
def drill_wedge_count(m, shell, band):
    """
    Wedges in row `band` of block shell `shell` at size `m`:

    - the root: `max(1, round(2*pi*R*sin(theta) / m))`, `theta` the row's
      centre colatitude, so a wedge's arc at the row's centre is about `m`
      edges.
    - every other level: the count of the parent row holding this one times
      `max(1, round(T / parent))`, `T` as above.

    Raises:
        ValueError: For a level that is not a block size or a row that does
            not exist.
    """
    _check_level(m)
    count = drill_band_count(m, shell)
    row = _row_of(m, shell, band)
    theta = (row + 0.5) * math.pi / count
    target = 2 * math.pi * block_radius(m, shell) * math.sin(theta) / m
    if m == DRILL_TOP:
        return max(1, _js_round(target))
    parent_band = _parent_band(m, shell, band)
    parent = drill_wedge_count(m * DRILL_STEP, shell // DRILL_STEP, parent_band)
    return parent * max(1, _js_round(target / parent))


def _parent_band(m, shell, band):
    """The signed row of the parent block shell holding row `band`."""
    count = drill_band_count(m, shell)
    parent_count = drill_band_count(m * DRILL_STEP, shell // DRILL_STEP)
    row = (count - 1) // 2 - band
    parent_row = row // (count // parent_count)
    return (parent_count - 1) // 2 - parent_row


def _centre_row(m, shell, colatitude):
    """The row (signed) of block shell `shell` at size `m` holding a
    colatitude; the poles go to the end rows."""
    count = drill_band_count(m, shell)
    row = min(count - 1, max(0, math.floor(colatitude * count / math.pi)))
    return (count - 1) // 2 - row


def cell_colatitude(shell, band):
    """The colatitude a sector's centre lies at: the middle of its band, the
    plane for band 0, the axis for a cap, and `pi / 2` for the core."""
    if shell == 0:
        return math.pi / 2
    if spherical.is_cap(shell, band):
        return 0.0 if band > 0 else math.pi
    if band == 0:
        return math.pi / 2
    north, south = spherical.band_colatitude_bounds(shell, band)
    return (north + south) / 2


def _cell_wedge_in(cell_wedges, cell_wedge, block_wedges):
    """The block wedge holding a sector wedge's centre longitude, in integers;
    a centre exactly on a meridian goes to the higher wedge."""
    return ((2 * cell_wedge + 1) * block_wedges) // (2 * cell_wedges)


def drill_parent(block):
    """
    The block one size up that holds `block`, or `None` for a block of the
    root level.

    Args:
        block (DrillBlock or tuple): `(m, shell, band, wedge)`, or a sector
            `(1, shell, band, wedge)`.

    Returns:
        DrillBlock or None.
    """
    m, shell, band, wedge = block
    if m == 1:
        parent_shell = shell // DRILL_STEP
        parent_band = _centre_row(DRILL_STEP, parent_shell, cell_colatitude(shell, band))
        if shell == 0 or spherical.is_cap(shell, band):
            return DrillBlock(DRILL_STEP, parent_shell, parent_band, 0)
        parent_wedges = drill_wedge_count(DRILL_STEP, parent_shell, parent_band)
        return DrillBlock(DRILL_STEP, parent_shell, parent_band,
                          _cell_wedge_in(spherical.wedge_count(shell, band), wedge, parent_wedges))
    _check_level(m)
    if m == DRILL_TOP:
        return None
    parent_shell = shell // DRILL_STEP
    parent_band = _parent_band(m, shell, band)
    q = drill_wedge_count(m, shell, band) // drill_wedge_count(m * DRILL_STEP, parent_shell, parent_band)
    return DrillBlock(m * DRILL_STEP, parent_shell, parent_band, wedge // q)


def stage_parent(block):
    """`drill_parent`, but `None` for a level-243 block, the top of the stages (`STAGE_LEVELS`): the parent
    a drill-down stage names, `None` for the galaxy's own view."""
    return None if block[0] == STAGE_TOP else drill_parent(block)


def drill_chain_of(shell, band, wedge):
    """
    Sector `(shell, band, wedge)`'s blocks from the root: `[level 177147, 59049,
    ... 9, 3, the sector itself]`, one per `DRILL_LEVELS`.
    """
    chain = [DrillBlock(1, shell, band, wedge)]
    while chain[-1].m != DRILL_TOP:
        chain.append(drill_parent(chain[-1]))
    return chain[::-1]


def _wedge_cells(cell_wedges, wedge, block_wedges):
    """The sector wedges `[first, last]` of a band whose centres fall in block
    wedge `wedge` of a row of `block_wedges` (an empty range when none do)."""
    first = max(0, -((block_wedges - 2 * wedge * cell_wedges) // (2 * block_wedges)))
    last = min(cell_wedges - 1, -((block_wedges - 2 * (wedge + 1) * cell_wedges) // (2 * block_wedges)) - 1)
    return first, last


def _bands_in(shell, row, rows):
    """The sector bands of `shell` (not the core), north to south, whose
    centre colatitude lies in row `row` of `rows` equal steps; the poles belong
    to the end rows."""
    low, high = spherical.band_range(shell)
    bands = list(range(high, low - 1, -1))
    centres = [cell_colatitude(shell, b) for b in bands]
    start = 0 if row == 0 else bisect_left(centres, row * math.pi / rows)
    stop = len(bands) if row == rows - 1 else bisect_left(centres, (row + 1) * math.pi / rows)
    return bands[start:stop]


def drill_block_sectors(block):
    """
    The sectors of level-3 block `block`, as `DrillBlock(1, shell, band,
    wedge)`: in each of its three shells, the bands whose centre lies in its
    row and, in each, the wedges whose centre lies in its wedge.
    """
    m, shell, band, wedge = block
    if m != DRILL_STEP:
        raise ValueError(f"only level-3 blocks hold sectors, got m={m}")
    rows = drill_band_count(m, shell)
    row = _row_of(m, shell, band)
    block_wedges = drill_wedge_count(m, shell, band)
    sectors = []
    for k in range(shell * DRILL_STEP, shell * DRILL_STEP + DRILL_STEP):
        if k == 0:
            if band == 0 and wedge == 0:
                sectors.append(DrillBlock(1, 0, 0, 0))
            continue
        for b in _bands_in(k, row, rows):
            if spherical.is_cap(k, b):
                if wedge == 0:
                    sectors.append(DrillBlock(1, k, b, 0))
                continue
            first, last = _wedge_cells(spherical.wedge_count(k, b), wedge, block_wedges)
            sectors.extend(DrillBlock(1, k, b, w) for w in range(first, last + 1))
    return sectors


def drill_children(block):
    """
    `block`'s children, grouped by child shell: `[(shell, [DrillBlock, ...]),
    ...]`, innermost shell first. A level-3 block's children are its sectors,
    grouped by sector shell (shells it holds none of are left out).
    """
    m, shell, band, wedge = block
    _check_level(m)
    f = DRILL_STEP
    if m == DRILL_STEP:
        by_shell = {}
        for sector in drill_block_sectors(block):
            by_shell.setdefault(sector.shell, []).append(sector)
        return sorted(by_shell.items())
    child_m = m // f
    parent_rows = drill_band_count(m, shell)
    parent_wedges = drill_wedge_count(m, shell, band)
    row = _row_of(m, shell, band)
    groups = []
    for child_shell in range(shell * f, shell * f + f):
        child_rows = drill_band_count(child_m, child_shell)
        step = child_rows // parent_rows
        blocks = []
        for child_row in range(row * step, row * step + step):
            child_band = (child_rows - 1) // 2 - child_row
            q = drill_wedge_count(child_m, child_shell, child_band) // parent_wedges
            blocks.extend(DrillBlock(child_m, child_shell, child_band, w) for w in range(wedge * q, wedge * q + q))
        groups.append((child_shell, blocks))
    return groups


def block_addresses(block):
    """
    Every sector address `(shell, band, wedge)` inside `block` of any level
    (a sector's own, for level 1), shell by shell outward: a level-3 block's
    sectors (`drill_block_sectors`), and a bigger block's children's,
    recursively.
    """
    m = block.m if hasattr(block, "m") else block[0]
    if m == 1:
        yield tuple(block)[1:]
        return
    if m == DRILL_STEP:
        for sector in drill_block_sectors(block):
            yield (sector.shell, sector.band, sector.wedge)
        return
    for _shell, children in drill_children(block):
        for child in children:
            yield from block_addresses(child)


def drill_block_bounds(block, edge_pc):
    """
    A block's analytic bounds: `((inner radius, outer radius), (north
    colatitude, south colatitude), (start longitude, end longitude))`, parsecs
    and radians. Block shell 0's inner radius is 0.
    """
    m, shell, band, wedge = block
    if m == 1:
        inner, outer = spherical.shell_bounds_pc(shell, edge_pc)
        north, south = spherical.band_colatitude_bounds(shell, band)
        start, end = spherical.wedge_longitude_bounds(shell, band, wedge)
        return (inner, outer), (north, south), (start, end)
    _check_level(m)
    rows = drill_band_count(m, shell)
    row = _row_of(m, shell, band)
    inner = 0.0 if shell == 0 else (shell * m - 0.5) * edge_pc
    outer = ((shell + 1) * m - 0.5) * edge_pc
    step = 2 * math.pi / drill_wedge_count(m, shell, band)
    return (inner, outer), (row * math.pi / rows, (row + 1) * math.pi / rows), (wedge * step, (wedge + 1) * step)


def format_drill_key(block):
    """`"m.shell.band.wedge"`, the stage URL's `at` value."""
    return ".".join(str(int(v)) for v in block)


def parse_drill_key(key):
    """
    `format_drill_key`'s inverse, checked: a known size, a block shell, and a
    row and wedge that exist at that size.

    Raises:
        ValueError: On anything else.
    """
    parts = str(key or "").split(".")
    if len(parts) != 4:
        raise ValueError(f"not a block key: {key!r}")
    try:
        m, shell, band, wedge = (int(p) for p in parts)
    except ValueError:
        raise ValueError(f"not a block key: {key!r}") from None
    if m not in DRILL_LEVELS[:-1]:
        raise ValueError(f"block size must be one of {DRILL_LEVELS[:-1]}, got {m}")
    if shell < 0:
        raise ValueError(f"no such block: {key!r}")
    rows = drill_band_count(m, shell)
    if abs(band) > rows // 2 or not 0 <= wedge < drill_wedge_count(m, shell, band):
        raise ValueError(f"no such block: {key!r}")
    return DrillBlock(m, shell, band, wedge)


def parse_node_key(key):
    """
    A node key `m.shell.band.wedge` of any level of the ladder, a sector (`1.shell.band.wedge`) included: what the
    node read and cache take (MAP.260), checked like `parse_drill_key`.

    Raises:
        ValueError: On anything else.
    """
    parts = str(key or "").split(".")
    if len(parts) == 4 and parts[0] == "1":
        try:
            _m, shell, band, wedge = (int(p) for p in parts)
        except ValueError:
            raise ValueError(f"not a node key: {key!r}") from None
        if shell < 0 or not spherical.band_range(shell)[0] <= band <= spherical.band_range(shell)[1] \
                or not 0 <= wedge < spherical.wedge_count(shell, band):
            raise ValueError(f"no such node: {key!r}")
        return DrillBlock(1, shell, band, wedge)
    return parse_drill_key(key)


def parse_stage_key(key):
    """`parse_drill_key`, and the block must be of a stage level (`STAGE_LEVELS`): the `at` of a drill-down stage.

    Raises:
        ValueError: On anything else.
    """
    block = parse_drill_key(key)
    check_stage_level(block.m)
    return block
