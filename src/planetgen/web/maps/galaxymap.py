# planetgen/web/maps/galaxymap.py

"""
Quadrant/Zone classification shared by every page that groups sectors by
galaxy-frame position -- `galaxy.py`'s own Quadrant summary/drill-down
tables, and `sector.py`/`browse.py`'s "Quadrant N" links back into it.

This module used to also build the Galaxy Map's own visualization (a
flat, face-on SVG projection) -- superseded by a real 3D map
(`planetgen/web/maps/galaxymap3d.py`/`static/galaxymap3d.js`), which `galaxy.py` now
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

from planetgen.galaxy import geometry
from planetgen.tuning import DEFAULT_SECTOR_EDGE_LY

QUADRANT_LABELS = ("I", "II", "III", "IV")
"""tuple[str]: The four galaxy-scale Quadrants, in azimuthal order starting
from `+X` -- `sector_quadrant` indexes into this. Deliberately the plain
2D I-IV convention (see this module's own docstring) rather than the
sector-internal octant's Roman-numeral scheme."""

ZONE_TARGET_LY = 100.0
"""float: The approximate light-year width a Zone should aim for.
`ZONE_RING_WIDTH` is derived from this so Zone boundaries land close to
round light-year milestones (~100 ly, ~200 ly, ...) while still being a
fixed number of the sector grid's own rings (`sectors.ring_index`, one
sector edge wide) end to end."""

ZONE_RING_WIDTH = max(1, round(ZONE_TARGET_LY / DEFAULT_SECTOR_EDGE_LY))
"""int: How many consecutive `ring_index` values make up one Zone."""


def sector_quadrant(x_pc, y_pc):
    """
    Labels a galaxy-frame `(x, y)` position with its Quadrant:
    `planetgen.galaxy.geometry.sector_quadrant`'s 1-4 (four 90-degree
    bands of `theta = atan2(y, x)` starting at `+X`, matching
    `docs/design/galaxy-coordinate-system.md`'s own theta convention) as
    `"I"`-`"IV"`. Blind to `z`/height by design -- a Quadrant is purely
    azimuthal, the same way real astronomical galactic quadrants are.

    Args:
        x_pc (float): `sectors.center_x_pc`.
        y_pc (float): `sectors.center_y_pc`.

    Returns:
        str: One of `QUADRANT_LABELS` (`"I"`, `"II"`, `"III"`, `"IV"`).
    """
    return QUADRANT_LABELS[geometry.sector_quadrant(x_pc, y_pc) - 1]


def sector_zone(ring_index):
    """The Zone index (a group of `ZONE_RING_WIDTH` consecutive rings)
    that `ring_index` falls in (`planetgen.galaxy.geometry.sector_zone`)."""
    return geometry.sector_zone(ring_index, DEFAULT_SECTOR_EDGE_LY, ZONE_TARGET_LY)


def zone_bounds_ly(zone_index):
    """
    The `(inner_ly, outer_ly)` cylindrical-radius bounds of Zone
    `zone_index` -- ring `i` spans `[i * edge, (i + 1) * edge)`, so these
    are exact ring boundaries, not the Zone's average radius.

    Returns:
        tuple[float, float]: `(inner_ly, outer_ly)`.
    """
    inner_ly = zone_index * ZONE_RING_WIDTH * DEFAULT_SECTOR_EDGE_LY
    outer_ly = (zone_index + 1) * ZONE_RING_WIDTH * DEFAULT_SECTOR_EDGE_LY
    return inner_ly, outer_ly
