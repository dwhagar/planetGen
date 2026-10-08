# planetgen/galaxy/system_position.py

"""
Where everything in a star system is in the galaxy (GEN.74).

A star, planet, moon or comet holds a `SpatialPosition3D` in AU. Its
"system" frame is the offset from its primary: the star for a planet or
comet (a binary's barycenter for a circumbinary body, the secondary star for
a wide binary's own planets and comets) and the parent planet for a moon.
A system made before it has a place keeps those anchors at the origin;
`place_system` gives it its place, and every body's offset stays as it is
while its sector and galactic coordinates follow.
"""

from planetgen.physics import constants as physical_constants
from planetgen.physics.position import SpatialPosition3D

_AXES = "xyz"


def _star_offsets_au(system):
    """`[(star, (x, y, z))]`: each star's offset from the system's center, AU."""
    stars = list(getattr(system, "stars", None) or [system.primary_star])
    if len(stars) < 2:
        return [(stars[0], (0.0, 0.0, 0.0))]
    if getattr(system, "binary_type", None) == "wide" and getattr(system, "wide_binary", None) is not None:
        pair, template = system.wide_binary, "{which}_position_{axis}_au"
    else:
        pair, template = system.star, "binary_{which}_position_{axis}"
    offsets = []
    for star, which in zip(stars, ("primary", "secondary")):
        offsets.append((star, tuple(float(getattr(pair, template.format(which=which, axis=axis), 0.0) or 0.0)
                                    for axis in _AXES)))
    return offsets


def _add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def _carry(body, sector_au, primary_au):
    """Gives a planet, moon or comet its anchors, keeping its own offset."""
    if getattr(body, "spatial", None) is None:  # an asteroid belt is a ring, with no point
        return None
    body.spatial.carry_anchors(sector_au, primary_au)
    return body.spatial.get_coordinates("galactic", "cartesian")


def place_system(system, sector_center_ly, position_ly):
    """
    Puts every body of `system` where it is in the galaxy.

    Args:
        system (StarSystem): The system (a stand-in without bodies is left
            alone).
        sector_center_ly (tuple): The sector's center, galactic light-years.
        position_ly (tuple): The system's center, light-years from the
            sector's center.
    """
    if not hasattr(system, "primary_star"):
        return
    to_au = physical_constants.LY_TO_AU
    sector_au = tuple(c * to_au for c in sector_center_ly)
    center_au = tuple((c + p) * to_au for c, p in zip(sector_center_ly, position_ly))

    star_centers = {}
    for star, offset in _star_offsets_au(system):
        at = _add(center_au, offset)
        star_centers[id(star)] = at
        mass = getattr(star, "mass", None)
        star.spatial = SpatialPosition3D(
            at, sector_au, is_star=True, length_unit_m=physical_constants.AU_M,
            mass_kg=mass if isinstance(mass, (int, float)) and mass >= 0 else None)

    wide = getattr(system, "binary_type", None) == "wide"
    primary_at = star_centers[id(system.primary_star)]
    secondary = getattr(system, "secondary_star", None)
    secondary_at = star_centers.get(id(secondary), center_au) if secondary is not None else center_au
    groups = [
        (getattr(system, "planets", []) or [], getattr(system, "comets", []) or [], primary_at if wide else center_au),
        (getattr(system, "secondary_planets", []) or [], getattr(system, "secondary_comets", []) or [], secondary_at),
    ]
    for planets, comets, anchor in groups:
        for planet in planets:
            planet_at = _carry(planet, sector_au, anchor)
            for moon in getattr(planet, "moons", []) or []:
                _carry(moon, sector_au, planet_at if planet_at is not None else anchor)
        for comet in comets:
            _carry(comet, sector_au, anchor)
