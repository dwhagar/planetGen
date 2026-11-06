# planetgen/web/maps/galaxymap.py

"""
Zone classification shared by the pages that group sectors by galaxy-frame position (`sector_zone`, `zone_bounds_ly`).

Galaxy quadrants are gone from the site: pages name the regions where sectors are charted instead (`db/region_stats.py`).
`star_systems.quadrant` (Roman numerals I-VIII, `spaceSector.classify_octant`) is a different, older concept -- the octant of a
*system's* position within its own *sector* -- and every page that shows it says "Octant".
"""

from planetgen.galaxy import spherical
from planetgen.tuning import DEFAULT_SECTOR_EDGE_LY

ZONE_TARGET_LY = 100.0
"""float: The approximate light-year width a Zone should aim for.
`ZONE_SHELL_DEPTH` is derived from this so Zone boundaries land close to
round light-year milestones (~100 ly, ~200 ly, ...) while still being a
fixed number of the sector grid's own shells (`sectors.shell_index`, one
sector edge deep) end to end."""

ZONE_SHELL_DEPTH = max(1, round(ZONE_TARGET_LY / DEFAULT_SECTOR_EDGE_LY))
"""int: How many consecutive `shell_index` values make up one Zone."""


def sector_zone(shell_index):
    """The Zone index (a group of `ZONE_SHELL_DEPTH` consecutive shells)
    that `shell_index` falls in (`planetgen.galaxy.spherical.sector_zone`)."""
    return spherical.sector_zone(shell_index, DEFAULT_SECTOR_EDGE_LY, ZONE_TARGET_LY)


def zone_bounds_ly(zone_index):
    """
    The `(inner_ly, outer_ly)` distance-from-the-center bounds of Zone
    `zone_index` -- shell `k` spans `[(k - 1/2) * edge, (k + 1/2) * edge)`
    (shell 0 is the core ball), so these are exact shell boundaries, not
    the Zone's average radius.

    Returns:
        tuple[float, float]: `(inner_ly, outer_ly)`.
    """
    inner_ly = max(0.0, (zone_index * ZONE_SHELL_DEPTH - 0.5) * DEFAULT_SECTOR_EDGE_LY)
    outer_ly = ((zone_index + 1) * ZONE_SHELL_DEPTH - 0.5) * DEFAULT_SECTOR_EDGE_LY
    return inner_ly, outer_ly
