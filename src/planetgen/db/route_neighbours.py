# planetgen/db/route_neighbours.py

"""
The systems round a point, for the graph round a route (NAV.65)
===============================================================

`nearby_systems` reads the placed systems nearest a point from the cell-key index (`db/region.py`), the same
covering index `near.objects_within` and `near.nearest` use, so its cost follows the systems near the point and never
the size of the galaxy. Only the id and the position are read: the name and sector of a neighbour are not needed to
draw a link. Radii climb the ladder `LADDER_PC` (the one `near.nearest` uses) and stop at the first sphere that holds
enough systems; a sphere over `near.SPHERE_LIMIT` systems (the dense core) is not read, the previous one answers.
"""

from planetgen.db import region
from planetgen.physics.units import ly_to_pc, pc_to_ly
from planetgen.util import detmath

LADDER_PC = (2.0, 4.0, 8.0, 16.0, 32.0, 50.0)
"""tuple: The radii tried, smallest first (`near.LADDER_PC`, then `near.MAX_DISTANCE_PC`)."""

SPHERE_LIMIT = 100_000
"""int: A sphere holding more systems than this is not read (`near.SPHERE_LIMIT`)."""

GUARDED_FROM_PC = 16.0
"""float: Spheres of this radius or more are counted before they are read. The densest place there is, the core at 15
systems per cubic parsec, holds about 500 systems in 2 pc, 4,000 in 4 pc and 32,000 in 8 pc, all under `SPHERE_LIMIT`, so
the smaller rungs need no count."""

_POSITION = ("p.galaxy_x_pc", "p.galaxy_y_pc", "p.galaxy_z_pc")


def nearby_systems(conn, point_ly, count, exclude_id=None):
    """
    The `count` placed systems nearest a galaxy-frame point, nearest first.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        point_ly (tuple): The `(x, y, z)` point, light-years, galaxy frame.
        count (int): How many systems.
        exclude_id (int, optional): A system left out (the stop's own).

    Returns:
        tuple: `(systems, reach_ly)`: `systems` is `{id: (x, y, z)}` in light-years, nearest first, at most `count`;
        `reach_ly` the distance within which the list is complete (the farthest listed system's when `count` were
        found, else the radius searched).
    """
    center = tuple(ly_to_pc(c) for c in point_ly)
    found, radius = [], 0.0
    for rung in LADDER_PC:
        if rung >= GUARDED_FROM_PC and region.capped_count(
                conn, "star_systems", _POSITION, center, rung, SPHERE_LIMIT + 1) > SPHERE_LIMIT:
            break
        rows = region.region_rows(conn, "star_systems", "p.id, p.galaxy_x_pc, p.galaxy_y_pc, p.galaxy_z_pc",
                                  _POSITION, center, rung)
        found, radius = [], rung
        for row in rows:
            if row["id"] == exclude_id:
                continue
            point = (row["galaxy_x_pc"], row["galaxy_y_pc"], row["galaxy_z_pc"])
            away = detmath.dist(center, point)
            if away <= rung:
                found.append((away, row["id"], point))
        if len(found) >= count:
            break
    found.sort(key=lambda entry: entry[:2])
    if len(found) >= count:
        found = found[:count]
        reach_pc = found[-1][0]
    else:
        reach_pc = radius
    return ({system_id: tuple(pc_to_ly(c) for c in point) for _away, system_id, point in found}, pc_to_ly(reach_pc))


def system_points(conn, system_ids):
    """
    The index position of each system given: `{id: (x, y, z)}` in light-years, galaxy frame (the columns
    `nearby_systems` reads, so a stop and its neighbours are measured from one source). A system with no galaxy
    position is left out.
    """
    points = {}
    ids = sorted(set(system_ids))
    for start in range(0, len(ids), 500):
        chunk = ids[start:start + 500]
        rows = conn.execute("SELECT id, galaxy_x_pc, galaxy_y_pc, galaxy_z_pc FROM star_systems WHERE id IN ("
                            + ", ".join("?" * len(chunk)) + ") AND galaxy_x_pc IS NOT NULL", tuple(chunk)).fetchall()
        for row in rows:
            points[row["id"]] = tuple(pc_to_ly(row[axis]) for axis in ("galaxy_x_pc", "galaxy_y_pc", "galaxy_z_pc"))
    return points
