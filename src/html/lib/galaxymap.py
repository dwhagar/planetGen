# html/lib/galaxymap.py

"""
Quadrant/Ring classification shared by every page that groups sectors by
galaxy-frame position -- `galaxy.py`'s own Quadrant summary/drill-down
tables, and `sector.py`/`browse.py`'s "Quadrant N" links back into it.

This module used to also build the Galaxy Map's own visualization (a
flat, face-on SVG projection) -- superseded by a real 3D map
(`lib/galaxymap3d.py`/`static/galaxymap3d.js`), which `galaxy.py` now
renders directly instead. That SVG rendering code (and its own
`static/galaxymap.js`) has been removed; only the plain classification
math below survived, since `sector.py`/`browse.py` still need it
independent of how the map itself is drawn.

Terminology note: `star_systems.quadrant` (Roman numerals I-VIII,
`spaceSector.classify_octant`) is a different, older concept -- an octant
classification of a *system's* position within its own *sector*. That
column/attribute name is left alone for schema stability, but every place
`html/` displays it now says "Octant" (`sector.py`, `system.py`,
`static/sectormap.js`) specifically so it doesn't collide with *this*
module's Quadrant, which is the real, galaxy-scale, 4-region azimuthal
concept the word ordinarily means (and which is a genuine 2D mathematical
standard, unlike the 8-region octant workaround -- see `spaceSector.py`'s
own module docstring).
"""

import math

try:
    from stellarObjects.program_constants import DEFAULT_SECTOR_EDGE_LY
except ImportError:
    # The planetGen package isn't on the import path in this deployment --
    # duplicated fallback, matching every other lib/ module's identical
    # pattern.
    DEFAULT_SECTOR_EDGE_LY = 11.5

QUADRANT_LABELS = ("I", "II", "III", "IV")
"""tuple[str]: The four galaxy-scale Quadrants, in azimuthal order starting
from `+X` -- `sector_quadrant` indexes into this. Deliberately the plain
2D I-IV convention (see this module's own docstring) rather than the
sector-internal octant's Roman-numeral scheme."""

RING_TARGET_LY = 100.0
"""float: The approximate light-year thickness a Ring should aim for.
`RING_SHELL_WIDTH` is derived from this so Ring boundaries land close to
round light-year milestones (~100 ly, ~200 ly, ...) while still being a
*fixed* shell-count width per Ring end to end (the preferred approach --
simple, uniform, and tracks the addressing scheme's own `shell_index`
directly, unlike a milestone-exact width, which would have to vary per
Ring)."""

RING_SHELL_WIDTH = max(1, round(RING_TARGET_LY / DEFAULT_SECTOR_EDGE_LY))
"""int: How many consecutive `shell_index` values make up one Ring."""


def sector_quadrant(x_pc, y_pc):
    """
    Classifies a galaxy-frame `(x, y)` position into its Quadrant --
    `theta = atan2(y, x)` normalized to `[0, 2*pi)`, split into four
    90-degree bands starting at `+X` (matching
    `docs/design/galaxy-coordinate-system.md`'s own theta convention).
    Blind to `z`/height by design -- a Quadrant is purely azimuthal, the
    same way real astronomical galactic quadrants are.

    Args:
        x_pc (float): `sectors.center_x_pc`.
        y_pc (float): `sectors.center_y_pc`.

    Returns:
        str: One of `QUADRANT_LABELS` (`"I"`, `"II"`, `"III"`, `"IV"`).
    """
    theta = math.atan2(y_pc, x_pc) % (2 * math.pi)
    index = min(3, int(theta // (math.pi / 2)))
    return QUADRANT_LABELS[index]


def sector_ring(shell_index):
    """The Ring index (a group of `RING_SHELL_WIDTH` consecutive shells)
    that `shell_index` falls in."""
    return shell_index // RING_SHELL_WIDTH


def ring_bounds_ly(ring_index):
    """
    The `(inner_ly, outer_ly)` light-year bounds of Ring `ring_index` --
    exact conversion of its shell-index band's own boundary radius (design
    doc section 3: shell `k` spans `[k*edge_pc, (k+1)*edge_pc)`), not the
    Ring's nominal/average radius.

    Returns:
        tuple[float, float]: `(inner_ly, outer_ly)`.
    """
    inner_ly = ring_index * RING_SHELL_WIDTH * DEFAULT_SECTOR_EDGE_LY
    outer_ly = (ring_index + 1) * RING_SHELL_WIDTH * DEFAULT_SECTOR_EDGE_LY
    return inner_ly, outer_ly
