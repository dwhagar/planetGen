"""
stellarObjects/galaxyDrill.py

The Galaxy Map drill-down's block ladder (docs/design/galaxy-drilldown-
navigation.md, section 3): blocks 243, 27, 3 and 1 sectors a side, nested
so that every block sits wholly inside one block of the next size up.
`static/galaxyprisms.js`'s `drill*` functions are the same rules for the
page; `tests/test_galaxydrill.py` checks that the two agree.

A block is `(m, ring, wedge, slab)`: `m` sectors a side, block ring `ring`
(sector rings `ring*m` to `ring*m + m - 1`), wedge `wedge` of
`drill_wedge_count(m, ring)` counterclockwise from +X, and slab `slab`
(sector layers `slab*m - (m-1)/2` to `slab*m + (m-1)/2`). At `m = 1` the
"block" is a sector: `(1, ring, slot, layer)`.

Rings and layers nest on their own (section 3.3). Wedges are made to nest
by the nested wedge rule (section 3.4): each child ring's wedge count is a
whole multiple of its parent ring's, so every parent wedge line is also a
child wedge line. Sectors join the level-3 block holding their center
(section 3.5), as `blockSlotRange` already does.

Everything here is integer math apart from the level-243 wedge count,
which copies `blockWedgeCount` step for step (including JavaScript's
round-half-up) so both sides pick the same count.
"""

import math
from collections import namedtuple
from functools import lru_cache

from stellarObjects.galaxyGeometry import ring_master_count, ring_sector_count

DRILL_LEVELS = (243, 27, 3, 1)
"""tuple: Block sizes from the galaxy down to a sector (Boss's "bigger
targets", 2026-10-01). Each step divides by 9, the last by 3."""

DRILL_TOP = DRILL_LEVELS[0]

BLOCK_WEDGE_TOLERANCE = 1.2
"""float: `galaxyprisms.js`'s own `BLOCK_WEDGE_TOLERANCE`, used by the
level-243 wedge count."""

DrillBlock = namedtuple("DrillBlock", ["m", "ring", "wedge", "slab"])
"""One block of the ladder; at `m = 1`, a sector `(1, ring, slot, layer)`."""


def _round_half_up(numerator, denominator):
    """`Math.round(numerator / denominator)` for integers, denominator > 0,
    without a float."""
    return (2 * numerator + denominator) // (2 * denominator)


def _js_round(x):
    """JavaScript's `Math.round`: halves round up, not to even."""
    return math.floor(x + 0.5)


def _target_wedges(ring):
    """`max(3, round(2*pi*(I + 1/2)))`, a block ring's aimed-for count."""
    return max(3, _js_round(2 * math.pi * (ring + 0.5)))


def _top_wedge_count(ring):
    """`blockWedgeCount(ring, 243)`: the master-aligned count nearest the
    target, if within `BLOCK_WEDGE_TOLERANCE`, else the target."""
    masters = ring_master_count(ring * DRILL_TOP)
    target = _target_wedges(ring)
    best = masters * max(1, _js_round(target / masters))
    d = 1
    while d < masters:
        for w in (d, 3 * d):
            if w < masters and abs(math.log(w / target)) < abs(math.log(best / target)):
                best = w
        d *= 2
    return best if abs(math.log(best / target)) <= math.log(BLOCK_WEDGE_TOLERANCE) else target


def _check_level(m):
    if m not in DRILL_LEVELS:
        raise ValueError(f"block size must be one of {DRILL_LEVELS}, got {m}")


def _step(m):
    """How many child blocks a side one block of size `m` holds."""
    return 3 if m == 3 else 9


@lru_cache(maxsize=65536)
def drill_wedge_count(m, ring):
    """
    Wedges in block ring `ring` at size `m` (the nested wedge rule):

    - 243: `blockWedgeCount(ring, 243)`, unchanged.
    - 27 and 3: the parent ring's count times
      `max(1, round(T / W_parent))`, `T` this ring's target
      `max(3, round(2*pi*(ring + 1/2)))`.
    - 1: the ring's slots, `ring_sector_count`.

    Args:
        m (int): One of `DRILL_LEVELS`.
        ring (int): Block ring `>= 0`.

    Returns:
        int: The wedge count.
    """
    _check_level(m)
    if ring < 0:
        raise ValueError(f"ring must be >= 0, got {ring}")
    if m == 1:
        return ring_sector_count(ring)
    if m == DRILL_TOP:
        return _top_wedge_count(ring)
    parent = drill_wedge_count(m * 9, ring // 9)
    return parent * max(1, _round_half_up(_target_wedges(ring), parent))


def drill_parent(block):
    """
    The block one size up that holds `block`, or `None` for a level-243
    block.

    Args:
        block (DrillBlock or tuple): `(m, ring, wedge, slab)`, or a sector
            `(1, ring, slot, layer)`.

    Returns:
        DrillBlock or None.
    """
    m, ring, wedge, slab = block
    _check_level(m)
    if m == DRILL_TOP:
        return None
    if m == 1:
        parent_ring = ring // 3
        # A sector joins the level-3 wedge holding its slot's center.
        parent_wedge = ((2 * wedge + 1) * drill_wedge_count(3, parent_ring)) // (2 * ring_sector_count(ring))
        return DrillBlock(3, parent_ring, parent_wedge, (slab + 1) // 3)
    parent_ring = ring // 9
    q = drill_wedge_count(m, ring) // drill_wedge_count(m * 9, parent_ring)
    return DrillBlock(m * 9, parent_ring, wedge // q, (slab + 4) // 9)


def drill_chain_of(ring, layer, slot):
    """
    Sector `(ring, layer, slot)`'s blocks from the top: `[level 243,
    level 27, level 3, the sector itself]`.
    """
    chain = [DrillBlock(1, ring, slot, layer)]
    while chain[-1].m != DRILL_TOP:
        chain.append(drill_parent(chain[-1]))
    return chain[::-1]


def drill_slabs(block, max_layer=None):
    """
    The child slabs of `block`, lowest first (child-level numbering; at
    `m = 3`, sector layers). For the galaxy (`block` is `None`), the
    level-243 slabs reaching `max_layer` (the outline's `galaxy_layer`),
    or just slab 0 without one.
    """
    if block is None:
        if max_layer is None:
            return [0]
        half = (DRILL_TOP - 1) // 2
        top = (max_layer + half) // DRILL_TOP
        return list(range(-top, top + 1))
    m, _ring, _wedge, slab = block
    _check_level(m)
    if m == 1:
        raise ValueError("a sector has no children")
    f = _step(m)
    half = (f - 1) // 2
    return list(range(slab * f - half, slab * f + half + 1))


def drill_block_sectors(block, layer):
    """
    The sectors of level-3 block `block` on sector layer `layer`, as
    `DrillBlock(1, ring, slot, layer)` -- the slots whose centers fall in
    its wedge, in each of its three rings (section 3.5).
    """
    m, ring, wedge, _slab = block
    if m != 3:
        raise ValueError(f"only level-3 blocks hold sectors, got m={m}")
    wedges = drill_wedge_count(3, ring)
    sectors = []
    for i in range(3 * ring, 3 * ring + 3):
        n = ring_sector_count(i)
        first = max(0, -((wedges - 2 * wedge * n) // (2 * wedges)))
        last = min(n - 1, -((wedges - 2 * (wedge + 1) * n) // (2 * wedges)) - 1)
        sectors.extend(DrillBlock(1, i, k, layer) for k in range(first, last + 1))
    return sectors


def drill_children(block, max_ring=None, max_layer=None):
    """
    `block`'s children, grouped by child slab: `[(slab, [DrillBlock,
    ...]), ...]`, lowest slab first. A level-3 block's children are its
    sectors, grouped by layer. For the galaxy (`block` is `None`), every
    level-243 block out to sector ring `max_ring` (required) on the slabs
    `drill_slabs(None, max_layer)` gives -- whether the stored outline
    allows any sector in them is the caller's test.
    """
    if block is None:
        if max_ring is None:
            raise ValueError("the galaxy's children need max_ring")
        rings = range(0, max_ring // DRILL_TOP + 1)
        return [
            (slab, [DrillBlock(DRILL_TOP, i, s, slab) for i in rings for s in range(drill_wedge_count(DRILL_TOP, i))])
            for slab in drill_slabs(None, max_layer)
        ]
    m, ring, wedge, _slab = block
    if m == 3:
        return [(layer, drill_block_sectors(block, layer)) for layer in drill_slabs(block)]
    f = _step(m)
    child_m = m // f
    parent_wedges = drill_wedge_count(m, ring)
    groups = []
    for slab in drill_slabs(block):
        blocks = []
        for i in range(ring * f, ring * f + f):
            q = drill_wedge_count(child_m, i) // parent_wedges
            blocks.extend(DrillBlock(child_m, i, s, slab) for s in range(wedge * q, wedge * q + q))
        groups.append((slab, blocks))
    return groups


def format_drill_key(block):
    """`"m.ring.wedge.slab"`, the stage URL's `at` value."""
    return ".".join(str(int(v)) for v in block)


def parse_drill_key(key):
    """
    `format_drill_key`'s inverse, checked: a known size, rings and wedges
    in range for that size. Slabs aren't checked against the outline.

    Raises:
        ValueError: On anything else.
    """
    parts = str(key or "").split(".")
    if len(parts) != 4:
        raise ValueError(f"not a block key: {key!r}")
    try:
        m, ring, wedge, slab = (int(p) for p in parts)
    except ValueError:
        raise ValueError(f"not a block key: {key!r}") from None
    if m not in DRILL_LEVELS[:-1]:
        raise ValueError(f"block size must be one of {DRILL_LEVELS[:-1]}, got {m}")
    if ring < 0 or not 0 <= wedge < drill_wedge_count(m, ring):
        raise ValueError(f"no such block: {key!r}")
    return DrillBlock(m, ring, wedge, slab)
