# planetgen/db/corridor.py

"""
What lies along a line (NAV.10)
===============================

Routing and course-steering need the objects near a straight segment of the
galaxy, not near a point. Both queries here cover the segment with boxes read
from `sectors` through its center index (`idx_sectors_center`), then read the
systems of those sectors by their indexed sector id and measure each exact
distance to the segment in the galaxy frame:

- `positions_near_segment`: the placed systems within a half-width of a
  segment, as the route search's node positions (light-years). It is what
  `query.nav_between` loads instead of every placed system in the galaxy.
- `objects_near_segment`: every system, star and phenomenon within a
  distance of a segment, in order along it (parsecs). NAV.6's gravity-well
  steering uses it.

A long segment is cut into pieces and each piece gets its own box, so the
read grows with the length of the line and not with the volume of the box
around it. Only generated objects count; nothing is generated here.
"""

import math

from planetgen import tuning
from planetgen.db.store import get_galaxy_shape
from planetgen.galaxy import objectref as object_ref
from planetgen.galaxy.geometry import sectors_along_segment
from planetgen.physics.units import ly_to_pc, mpc_to_pc, pc_to_ly

_CHUNK = 500
"""int: Ids per `IN (...)` list."""

PIECE_RATIO = 4.0
"""float: A piece of the segment is at most this many reaches long, so its
box wastes little volume."""

MAX_PIECES = 4000
"""int: Pieces a segment is cut into at most (a longer one gets longer
pieces)."""


def segment_distance(point, start, end):
    """`(distance, along)` from `point` to the segment `start`-`end`: the
    shortest distance, and how far along the segment (0 at `start`, in the
    segment's own units) the nearest point on it lies."""
    ax, ay, az = start
    dx, dy, dz = end[0] - ax, end[1] - ay, end[2] - az
    length_squared = dx * dx + dy * dy + dz * dz
    if length_squared == 0.0:
        return math.dist(point, start), 0.0
    t = ((point[0] - ax) * dx + (point[1] - ay) * dy + (point[2] - az) * dz) / length_squared
    t = max(0.0, min(1.0, t))
    nearest = (ax + t * dx, ay + t * dy, az + t * dz)
    return math.dist(point, nearest), t * math.sqrt(length_squared)


def _pieces(start, end, reach):
    """The segment cut into `(midpoint, half_length)` pieces."""
    length = math.dist(start, end)
    count = max(1, min(MAX_PIECES, math.ceil(length / (PIECE_RATIO * max(reach, 1e-9)))))
    pieces = []
    for index in range(count):
        low, high = index / count, (index + 1) / count
        middle = tuple(start[axis] + (end[axis] - start[axis]) * (low + high) / 2.0 for axis in range(3))
        pieces.append((middle, length * (high - low) / 2.0))
    return pieces


def _sector_edge_pc(conn):
    shape = get_galaxy_shape(conn)
    return float(shape.edge_pc) if shape is not None else float(tuning.DEFAULT_SECTOR_EDGE_PC)


def generated_cells(conn, cells):
    """The subset of `(ring_index, layer_index, ring_slot_index)` cells that
    hold a generated sector."""
    by_band = {}
    for ring, layer, slot in cells:
        by_band.setdefault((ring, layer), []).append(slot)
    found = set()
    for (ring, layer), slots in by_band.items():
        for start in range(0, len(slots), _CHUNK):
            chunk = slots[start:start + _CHUNK]
            marks = ",".join("?" * len(chunk))
            rows = conn.execute(
                f"SELECT ring_slot_index FROM sectors WHERE ring_index = ? AND layer_index = ? "
                f"AND ring_slot_index IN ({marks})", (ring, layer, *chunk)).fetchall()
            found.update((ring, layer, row["ring_slot_index"]) for row in rows)
    return found


def unknown_space_flags(conn, points_ly):
    """
    NAV.12: for each hop of a route (`points_ly`, the stops' galaxy-frame
    light-year positions in order), whether its straight line crosses a sector
    that has not been generated -- "a jump through unknown space". The cells a
    line crosses come from `sectors_along_segment` (NAV.38); a cell is known
    when `sectors` holds a row at its address.

    Returns:
        list[bool]: One flag per hop (`len(points_ly) - 1`).
    """
    edge = _sector_edge_pc(conn)
    per_hop = [
        sectors_along_segment(tuple(ly_to_pc(c) for c in a), tuple(ly_to_pc(c) for c in b), edge)
        for a, b in zip(points_ly, points_ly[1:])
    ]
    known = generated_cells(conn, {cell for cells in per_hop for cell in cells})
    return [any(cell not in known for cell in cells) for cells in per_hop]


def _sectors_near(conn, start, end, distance, edge):
    """`{sector_id: (center, name)}` of every placed sector whose cell a
    system within `distance` of the segment can lie in."""
    reach = distance + edge * math.sqrt(3.0) / 2.0 + edge * 0.5
    found = {}
    for middle, half in _pieces(start, end, reach):
        box = half + reach
        rows = conn.execute(
            "SELECT id, name, center_x_pc, center_y_pc, center_z_pc FROM sectors"
            " WHERE center_x_pc BETWEEN ? AND ? AND center_y_pc BETWEEN ? AND ? AND center_z_pc BETWEEN ? AND ?",
            (middle[0] - box, middle[0] + box, middle[1] - box, middle[1] + box,
             middle[2] - box, middle[2] + box)).fetchall()
        for row in rows:
            center = (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])
            if segment_distance(center, start, end)[0] <= reach:
                found[row["id"]] = (center, row["name"])
    return found


def _systems_of(conn, sectors, start, end, distance):
    """The placed systems in `sectors` within `distance` of the segment:
    `{id: {"name", "sector_id", "point", "distance", "along"}}` (parsecs)."""
    found = {}
    ids = list(sectors)
    for offset in range(0, len(ids), _CHUNK):
        chunk = ids[offset:offset + _CHUNK]
        rows = conn.execute(
            "SELECT id, name, sector_id, position_x_mpc, position_y_mpc, position_z_mpc FROM star_systems"
            f" WHERE sector_id IN ({', '.join('?' * len(chunk))}) AND position_x_mpc IS NOT NULL", tuple(chunk)
        ).fetchall()
        for row in rows:
            center = sectors[row["sector_id"]][0]
            point = (center[0] + mpc_to_pc(row["position_x_mpc"]), center[1] + mpc_to_pc(row["position_y_mpc"]),
                     center[2] + mpc_to_pc(row["position_z_mpc"]))
            away, along = segment_distance(point, start, end)
            if away <= distance:
                found[row["id"]] = {"name": row["name"], "sector_id": row["sector_id"], "point": point,
                                    "distance": away, "along": along}
    return found


def positions_near_segment(conn, start_ly, end_ly, half_width_ly):
    """
    The placed systems within `half_width_ly` of the segment `start_ly` to
    `end_ly` (galaxy frame): `{star_systems.id: (x, y, z)}` in light-years,
    the shape `query.nav_between` builds its route graph from.
    """
    start, end = (tuple(ly_to_pc(value) for value in point) for point in (start_ly, end_ly))
    distance = ly_to_pc(half_width_ly)
    sectors = _sectors_near(conn, start, end, distance, _sector_edge_pc(conn))
    systems = _systems_of(conn, sectors, start, end, distance)
    return {system_id: tuple(pc_to_ly(value) for value in system["point"]) for system_id, system in systems.items()}


def objects_near_segment(conn, start_pc, end_pc, distance_pc, kinds=None):
    """
    Every generated system, star and phenomenon within `distance_pc` of the
    segment `start_pc` to `end_pc` (galaxy frame, parsecs), in order along it
    (ties: nearer to the line, then kind and id).

    Args:
        kinds (iterable, optional): Any of `"system"`, `"star"` and
            `planetgen.galaxy.objectref.PHENOMENON_KINDS`; all by default.

    Returns:
        list[dict]: `ref`, `kind`, `id`, `name`, `system_id` (the system it
            is or belongs to, else `None`), `distance_pc` (to the line),
            `along_pc` (where along the line its nearest point lies) and
            `point_pc`.

    Raises:
        ValueError: For a negative distance or a kind that cannot be used.
    """
    if distance_pc < 0:
        raise ValueError("the distance is 0 or more")
    allowed = {"system", "star", *object_ref.PHENOMENON_KINDS}
    wanted = set(allowed if kinds is None else kinds)
    if not wanted <= allowed:
        raise ValueError(f"unknown kinds: {', '.join(sorted(wanted - allowed))}")
    start, end = tuple(start_pc), tuple(end_pc)
    found = []
    if wanted & {"system", "star"}:
        sectors = _sectors_near(conn, start, end, distance_pc, _sector_edge_pc(conn))
        systems = _systems_of(conn, sectors, start, end, distance_pc)
        for system_id, system in systems.items():
            if "system" in wanted:
                found.append({"ref": f"system:{system_id}", "kind": "system", "id": system_id,
                              "name": system["name"], "system_id": system_id,
                              "distance_pc": system["distance"], "along_pc": system["along"],
                              "point_pc": system["point"]})
        if "star" in wanted:
            ids = list(systems)
            for offset in range(0, len(ids), _CHUNK):
                chunk = ids[offset:offset + _CHUNK]
                rows = conn.execute(
                    "SELECT id, name, star_system_id FROM stars"
                    f" WHERE star_system_id IN ({', '.join('?' * len(chunk))})", tuple(chunk)).fetchall()
                for row in rows:
                    system = systems[row["star_system_id"]]
                    found.append({"ref": f"star:{row['id']}", "kind": "star", "id": row["id"], "name": row["name"],
                                  "system_id": row["star_system_id"], "distance_pc": system["distance"],
                                  "along_pc": system["along"], "point_pc": system["point"]})
    phenomena = wanted & set(object_ref.PHENOMENON_KINDS)
    if phenomena:
        from planetgen.db import query  # not at the top: query imports this module
        widest = max(query._widest_placed_phenomenon_radius_ly(conn), 0.0)
        margin = ly_to_pc(widest) if widest > 0 else 0.0
        reach = distance_pc + margin
        seen = set()
        for middle, half in _pieces(start, end, reach):
            box = (middle[0], middle[1], middle[2], half + reach)
            for row in query._placed_phenomenon_rows(conn, bbox=box):
                key = (row["type"], row["id"])
                if row["type"] not in phenomena or key in seen:
                    continue
                point = (row["x"], row["y"], row["z"])
                away, along = segment_distance(point, start, end)
                if away <= distance_pc:
                    seen.add(key)
                    found.append({"ref": f"{row['type']}:{row['id']}", "kind": row["type"], "id": row["id"],
                                  "name": row["name"], "system_id": None, "distance_pc": away,
                                  "along_pc": along, "point_pc": point})
    order = {"system": 0, "star": 1}
    found.sort(key=lambda row: (row["along_pc"], row["distance_pc"], order.get(row["kind"], 2), row["id"]))
    return found
