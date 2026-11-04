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

import math

from planetgen.physics import system_tree
from planetgen.physics.position import SpatialPosition3D
from planetgen.util import draw
from planetgen.util import detmath
from planetgen.physics import unit_registry as units

_kg_to_msun = units.converter("kg", "Msun")

_ly_to_au = units.converter("ly", "au")

_AXES = "xyz"


def _star_offsets_au(system):
    """`[(star, (x, y, z))]`: each star's offset from the system's center, AU."""
    stars = list(getattr(system, "stars", None) or [system.primary_star])
    if getattr(system, "star_slots", None) is not None:    # GEN.630: more than a pair, placed by the whole tree
        dark = [(member.mass_solar, (member.parent_node_no, member.slot)) for member in getattr(system, "dark_members", ())]
        return list(zip(stars, system_tree.star_offsets_au(
            [_kg_to_msun(star.mass) for star in stars], system.star_slots, system.orbit_nodes,
            system.orbit_nodes[0].epoch_unix, dark)))
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


def galactic_velocity_ms(position, speed_kms):
    """
    The velocity, galactic axes, m/s, of something on the galaxy's circular
    rotation (`galaxy/galactic_orbit.py`): `speed_kms` along the tangent at
    `position` (any galactic unit), counterclockwise seen from galactic
    north, the way `db.store.advance_galactic_positions` turns it. At rest
    for no speed or on the axis.
    """
    if not isinstance(speed_kms, (int, float)) or isinstance(speed_kms, bool) or not speed_kms > 0.0:
        return (0.0, 0.0, 0.0)
    r = detmath.hypot(position[0], position[1])
    if r == 0.0:
        return (0.0, 0.0, 0.0)
    scale = speed_kms * 1000.0 / r
    return (-position[1] * scale, position[0] * scale, 0.0)


def random_unit_vector(rng=draw):
    """A direction uniform on the sphere, from `rng.random()` alone (the
    seeded generators replace that one function -- GEN.56); `rng` is the
    unit's draw stream unless a `draw.Stream` is passed."""
    z = 2.0 * rng.random() - 1.0
    phi = 2.0 * math.pi * rng.random()
    ring = math.sqrt(max(0.0, 1.0 - z * z))
    return (ring * math.cos(phi), ring * math.sin(phi), z)


def peculiar_velocity_ms(system):
    """
    The velocity, m/s, a runaway or hypervelocity system has beyond the
    galaxy's rotation: its `runaway_speed_kms` along its `runaway_direction`
    (galactic axes). At rest for an ordinary system.
    """
    speed = getattr(system, "runaway_speed_kms", None)
    direction = getattr(system, "runaway_direction", None)
    if (not isinstance(speed, (int, float)) or isinstance(speed, bool) or not speed > 0.0
            or not isinstance(direction, (tuple, list)) or len(direction) != 3):
        return (0.0, 0.0, 0.0)
    return tuple(speed * 1000.0 * float(c) for c in direction)


def system_velocity_ms(system, position):
    """The velocity, galactic axes, m/s, of `system` at `position` (any
    galactic unit): the rotation curve's tangent plus its runaway motion."""
    speed = getattr(getattr(system, "star", None), "galactic_orbital_speed_kms", None)
    return _add(galactic_velocity_ms(position, speed), peculiar_velocity_ms(system))


def set_system_epoch(system, epoch_unix):
    """Stamps every star, planet, moon and comet of `system` that holds a
    position with the time it holds at."""
    if not hasattr(system, "primary_star"):
        return
    bodies = list(getattr(system, "stars", None) or [system.primary_star])
    for planet in getattr(system, "all_planets", None) or (
            list(getattr(system, "planets", []) or []) + list(getattr(system, "secondary_planets", []) or [])):
        bodies.append(planet)
        bodies.extend(getattr(planet, "moons", []) or [])
    bodies.extend(getattr(system, "comets", []) or [])
    bodies.extend(getattr(system, "secondary_comets", []) or [])
    for body in bodies:
        spatial = getattr(body, "spatial", None)
        if spatial is not None:
            spatial.set_epoch_unix(epoch_unix)


def _carry(body, sector_au, primary_au, sector_edge_pc, primary_velocity_ms, report=False):
    """Gives a planet, moon or comet its anchors (and its primary's velocity),
    keeping its own offset and its velocity relative to the primary. Returns
    where it is and how it moves when `report` (a planet's moons are placed
    from it); the others' coordinates are worked out when first read."""
    if getattr(body, "spatial", None) is None:  # an asteroid belt is a ring, with no point
        return None
    body.spatial.set_sector_edge_pc(sector_edge_pc)
    body.spatial.carry_anchors(sector_au, primary_au)
    body.spatial.carry_star_velocity(primary_velocity_ms)
    if not report:
        return None
    return body.spatial.get_coordinates("galactic", "cartesian"), body.spatial.get_velocity_vector("galactic")


def place_system(system, sector_center_ly, position_ly, sector_edge_pc=None, velocity_ms=None):
    """
    Puts every body of `system` where it is in the galaxy.

    The system moves on the galaxy's rotation curve
    (`galactic_velocity_ms`: its `galactic_orbital_speed_kms` along the
    tangent at its center) plus its runaway motion
    (`peculiar_velocity_ms`); every star takes that one velocity, and every
    body carries its primary's velocity with the velocity it has relative to
    it.

    Args:
        system (StarSystem): The system (a stand-in without bodies is left
            alone).
        sector_center_ly (tuple): The sector's center, galactic light-years.
        position_ly (tuple): The system's center, light-years from the
            sector's center.
        sector_edge_pc (float or None): The sector grid's edge, parsecs, so
            every body knows its sector address; `None` when the sector is
            not a cell of the galaxy's grid.
        velocity_ms (tuple or None): The system center's velocity, galactic
            axes, m/s, when it is known (as stored); `None` works it out
            from the rotation curve and the system's runaway motion.

    Returns:
        tuple: The system center's velocity, galactic axes, m/s.
    """
    if not hasattr(system, "primary_star"):
        return None
    to_au = _ly_to_au(1.0)
    sector_au = tuple(c * to_au for c in sector_center_ly)
    center_au = tuple((c + p) * to_au for c, p in zip(sector_center_ly, position_ly))

    system_speed = getattr(getattr(system, "star", None), "galactic_orbital_speed_kms", None)
    circular_velocity = galactic_velocity_ms(center_au, system_speed)
    system_velocity = tuple(velocity_ms) if velocity_ms is not None else _add(circular_velocity,
                                                                              peculiar_velocity_ms(system))
    # Every star of a system moves with the system: the pair's bulk motion is one vector (a wide pair's two stars sit
    # up to a few hundredths of a parsec apart, where the rotation curve's tangent turns by about a millionth, 0.06 m/s
    # of 200 km/s; a star taking the tangent at its own place instead would load back, from the one stored vector,
    # differing by that much, GEN.527).
    star_centers, star_velocities = {}, {}
    for star, offset in _star_offsets_au(system):
        at = _add(center_au, offset)
        star_centers[id(star)] = at
        star_velocities[id(star)] = system_velocity
        mass = getattr(star, "mass", None)
        star.spatial = SpatialPosition3D(
            at, sector_au, is_star=True, length_unit="au", sector_edge_pc=sector_edge_pc,
            mass_kg=mass if isinstance(mass, (int, float)) and mass >= 0 else None,
            velocity_vector_cartesian=star_velocities[id(star)])

    wide = getattr(system, "binary_type", None) == "wide"
    # (a lone star with dark members about it is off the tree's barycentre, and its planets go round the star, GEN.1013)
    hosted = wide or (getattr(system, "binary_type", None) is None and bool(getattr(system, "dark_members", ())))
    primary_at = star_centers[id(system.primary_star)]
    secondary = getattr(system, "secondary_star", None)
    secondary_at = star_centers.get(id(secondary), center_au) if secondary is not None else center_au
    primary_v = star_velocities[id(system.primary_star)]
    secondary_v = star_velocities.get(id(secondary), system_velocity) if secondary is not None else system_velocity
    groups = [
        (getattr(system, "planets", []) or [], getattr(system, "comets", []) or [],
         primary_at if hosted else center_au, primary_v if hosted else system_velocity),
        (getattr(system, "secondary_planets", []) or [], getattr(system, "secondary_comets", []) or [],
         secondary_at, secondary_v),
    ]
    for member, bodies in system.member_hosts() if hasattr(system, "member_hosts") else ():    # GEN.332
        groups.append((bodies, [], star_centers.get(id(member), center_au), star_velocities.get(id(member), system_velocity)))
    for planets, comets, anchor, anchor_v in groups:
        for planet in planets:
            placed = _carry(planet, sector_au, anchor, sector_edge_pc, anchor_v, report=True)
            planet_at, planet_v = placed if placed is not None else (anchor, anchor_v)
            for moon in getattr(planet, "moons", []) or []:
                _carry(moon, sector_au, planet_at, sector_edge_pc, planet_v)
        for comet in comets:
            _carry(comet, sector_au, anchor, sector_edge_pc, anchor_v)
    return system_velocity
