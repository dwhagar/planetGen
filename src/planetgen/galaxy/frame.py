# planetgen/galaxy/frame.py

"""
The Galaxy Frame and a Sector's Local Axes
=============================================

Distances from the galactic axis and center, and the local axes a sector
carries (`sector_orientation`), which the offsets of `star_systems.position_x/y/z_mpc` and the
Sector Map's scenes are expressed in. None of it depends on how the galaxy is
cut into sectors (`galaxy/spherical.py`). Every distance is in parsecs.
Pure geometry: no database, no I/O.
"""

import math

from planetgen.util import detmath
from planetgen.physics import unit_registry as units

_km_to_pc = units.converter("km", "pc")


def galactic_radius_pc(position):
    """
    `sqrt(x^2 + y^2 + z^2)` -- straight-line distance from the galactic
    center, matching `sectors.galactic_radius_pc`.
    """
    x, y, z = position
    return math.sqrt(x * x + y * y + z * z)


def sector_orientation(center_pc):
    """
    A sector's local axes in the galaxy frame, the frame
    `star_systems.position_x/y/z_mpc` offsets are expressed in: local
    `+X` points radially outward from the galactic axis, local `+Y` points
    toward increasing longitude (counterclockwise seen from galactic north),
    and local `+Z` is galactic north. The same convention for every
    sector, so a system's octant means the same thing everywhere.

    Args:
        center_pc (tuple): The sector's own `(x, y, z)` center, parsecs.
            Local `+X` is the direction of the center's longitude, so a
            center on the galactic axis (a polar cap, the core, or a
            hand-placed sector) gets the galaxy's own axes.

    Returns:
        tuple: `(local_x, local_y, local_z)`, each a unit `(x, y, z)`.
    """
    x, y, _z = center_pc
    r = detmath.hypot(x, y)
    if r < 1e-12:
        return (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)
    cos_t, sin_t = x / r, y / r
    return (cos_t, sin_t, 0.0), (-sin_t, cos_t, 0.0), (0.0, 0.0, 1.0)


def local_to_galaxy_pc(center_pc, offset_pc):
    """Converts a sector-local offset (parsecs, `sector_orientation`'s
    axes) into an absolute galaxy-frame point."""
    ax, ay, az = sector_orientation(center_pc)
    ox, oy, oz = offset_pc
    return tuple(center_pc[i] + ox * ax[i] + oy * ay[i] + oz * az[i] for i in range(3))


BINARY_PAIR_DIRECTION = (math.sqrt(0.5), -math.sqrt(0.5), 0.0)
"""tuple: The galaxy-frame unit vector a binary's secondary star lies along
from its primary, in every map (MAP.186) -- the orbit's own orientation
is not stored per pair, so all pairs share one direction."""



def binary_companion_offset_pc(separation_km):
    """
    Where a binary's secondary star sits from its primary, galaxy frame,
    parsecs: the stored `binary_separation_km` along `BINARY_PAIR_DIRECTION`
    (MAP.186), the one rule the Galaxy Map's nodes and the Sector Map's
    scene share. `(0, 0, 0)` when the separation is unknown.
    """
    scale = _km_to_pc(separation_km or 0.0)
    return tuple(scale * component for component in BINARY_PAIR_DIRECTION)


def galaxy_to_local_pc(center_pc, point_pc):
    """`local_to_galaxy_pc`'s inverse: an absolute galaxy-frame point as
    an offset (parsecs) from the sector's center along its local axes."""
    axes = sector_orientation(center_pc)
    delta = tuple(point_pc[i] - center_pc[i] for i in range(3))
    return tuple(sum(delta[i] * axis[i] for i in range(3)) for axis in axes)
