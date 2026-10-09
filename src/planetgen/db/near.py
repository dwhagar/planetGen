# planetgen/db/near.py

"""
What is near a place (NAV.43)
=============================

`objects_within` answers "what is within N parsecs of this place?" for every
kind of object: star systems with their stars, planets, moons, asteroid belts
and comets (placed at their system's position), facilities, and the
standalone phenomena. A place is an object reference (`planetgen.galaxy.
objectref`, a bare number being a system) or a point in the galaxy frame
(parsecs).

Only generated objects count. The sectors the sphere reaches that have not
been generated are counted, never filled: the search reads, it does not
generate.

Distances are measured to an object's position -- the center of a nebula or
supernova remnant, not its edge. Rows come nearest first (ties: system, its
bodies, facilities, then the phenomena, each by id), and the result pages
like every list (`limit`, `offset`).

How it stays cheap: the sectors the sphere can reach come from one bounding
box read of `sectors` (`idx_sectors_center`); systems, rogue objects and
facilities are read by their indexed sector or system ids, phenomena by
their own indexed center columns, and exact distances are measured in the
galaxy frame. A system's bodies are never loaded for the whole sphere: only
a per-system count is, and the bodies of the systems that fall on the
requested page are read afterwards.
"""

import math

from planetgen import tuning
from planetgen.db import query
from planetgen.db.store import get_galaxy_shape
from planetgen.galaxy import objectref as object_ref
from planetgen.galaxy.geometry import enumerate_sectors_within_radius
from planetgen.physics.units import ly_to_pc, mpc_to_pc

MAX_DISTANCE_PC = 50.0
"""float: The largest distance a search takes (Boss, NAV.43: about 8,000
sectors at 4 pc a side)."""

DEFAULT_LIMIT = 50
MAX_LIMIT = 200

FACILITY_KIND = "facility"

BODY_KINDS = ("star", "planet", "moon", "belt", "comet")
"""tuple: The kinds read from a system's own tables, in the order they sort
under their system."""

SEARCH_KINDS = ("system",) + BODY_KINDS + (FACILITY_KIND,) + object_ref.PHENOMENON_KINDS
"""tuple: Every kind the search can return (and filter on)."""

_CHUNK = 500
"""int: Ids per `IN (...)` list."""


class NearError(ValueError):
    """A place or distance the search cannot use."""


def _chunks(items, size=_CHUNK):
    items = list(items)
    for start in range(0, len(items), size):
        yield items[start:start + size]


def _marks(count):
    return ", ".join(["?"] * count)


def place_from_reference(conn, reference):
    """
    Resolves an object reference to a place.

    Returns:
        dict: `ref`, `kind`, `id`, `name` and `point_pc` (galaxy frame,
            parsecs).

    Raises:
        NearError: For a reference that is malformed, names nothing, or
            names something not placed in the galaxy.
    """
    try:
        kind, object_id = object_ref.parse(reference)
        found = query.resolve_object(conn, kind, object_id)
    except ValueError as exc:
        raise NearError(str(exc)) from exc
    point = found["positions"]["galaxy_pc"]
    if point is None:
        raise NearError(f"{found['ref']} is not placed in the galaxy, so there is nothing to measure from")
    return {"ref": found["ref"], "kind": found["kind"], "id": object_id, "name": found["name"],
            "point_pc": tuple(float(v) for v in point)}


def place_from_point(point):
    """A bare galaxy-frame point `(x, y, z)` in parsecs as a place.

    Raises:
        NearError: For anything but three finite numbers."""
    try:
        values = tuple(float(v) for v in point)
    except (TypeError, ValueError) as exc:
        raise NearError("a point is three numbers: x, y, z in parsecs") from exc
    if len(values) != 3 or not all(math.isfinite(v) for v in values):
        raise NearError("a point is three finite numbers: x, y, z in parsecs")
    return {"ref": None, "kind": "point", "id": None, "name": "(%.2f, %.2f, %.2f) pc" % values, "point_pc": values}


def _check_distance(distance_pc):
    try:
        distance = float(distance_pc)
    except (TypeError, ValueError) as exc:
        raise NearError("the distance must be a number of parsecs") from exc
    if not math.isfinite(distance) or distance <= 0.0:
        raise NearError("the distance must be a number of parsecs greater than 0")
    if distance > MAX_DISTANCE_PC:
        raise NearError(f"the largest distance is {MAX_DISTANCE_PC:g} pc")
    return distance


def _check_kinds(kinds):
    if kinds is None:
        return set(SEARCH_KINDS)
    wanted = set(kinds)
    unknown = sorted(wanted - set(SEARCH_KINDS))
    if unknown:
        raise NearError(f"unknown kind {unknown[0]!r}; the kinds are {', '.join(SEARCH_KINDS)}")
    return wanted


def _ref(kind, object_id):
    return f"{kind}:{object_id}"


def _systems_in_reach(conn, center, distance, margin):
    """Every placed system whose position is within `distance` of `center`:
    `{id: {...}}` with `name`, `sector_id`, `sector_name`, `point`, `distance`."""
    reach = distance + margin
    rows = conn.execute(
        """
        SELECT ss.id, ss.name, ss.sector_id, sec.name AS sector_name,
               sec.center_x_pc, sec.center_y_pc, sec.center_z_pc,
               ss.position_x_mpc, ss.position_y_mpc, ss.position_z_mpc
        FROM sectors sec JOIN star_systems ss ON ss.sector_id = sec.id
        WHERE sec.center_x_pc BETWEEN ? AND ? AND sec.center_y_pc BETWEEN ? AND ?
          AND sec.center_z_pc BETWEEN ? AND ? AND ss.position_x_mpc IS NOT NULL
        """,
        (center[0] - reach, center[0] + reach, center[1] - reach, center[1] + reach,
         center[2] - reach, center[2] + reach),
    ).fetchall()
    found = {}
    for row in rows:
        point = (row["center_x_pc"] + mpc_to_pc(row["position_x_mpc"]),
                 row["center_y_pc"] + mpc_to_pc(row["position_y_mpc"]),
                 row["center_z_pc"] + mpc_to_pc(row["position_z_mpc"]))
        away = math.dist(center, point)
        if away <= distance:
            found[row["id"]] = {"name": row["name"], "sector_id": row["sector_id"], "sector_name": row["sector_name"],
                                "point": point, "distance": away}
    return found


def _sector_counts(conn, center, distance, margin):
    """`(reached, generated)`: how many grid sectors the sphere touches and
    how many of them have been generated."""
    shape = get_galaxy_shape(conn)
    if shape is None:
        return 0, 0
    bounds = query.get_galaxy_bounds(conn)
    reached = {(ring, layer, slot)
               for ring, layer, slot, *_rest in enumerate_sectors_within_radius(center, distance, float(shape.edge_pc))
               if bounds.contains(ring, layer)}
    reach = distance + margin
    rows = conn.execute(
        """
        SELECT ring_index, layer_index, ring_slot_index FROM sectors
        WHERE center_x_pc BETWEEN ? AND ? AND center_y_pc BETWEEN ? AND ? AND center_z_pc BETWEEN ? AND ?
        """,
        (center[0] - reach, center[0] + reach, center[1] - reach, center[1] + reach,
         center[2] - reach, center[2] + reach)).fetchall()
    made = {(row["ring_index"], row["layer_index"], row["ring_slot_index"]) for row in rows}
    return len(reached), len(reached & made)


def _body_counts(conn, system_ids, wanted):
    """`{system_id: {kind: count}}` for the wanted body kinds."""
    tables = {"star": "stars", "planet": "planets", "moon": "moons", "belt": "asteroid_belts", "comet": "comets"}
    counts = {}
    for kind in BODY_KINDS:
        if kind not in wanted:
            continue
        for chunk in _chunks(system_ids):
            rows = conn.execute(
                f"SELECT star_system_id AS sid, COUNT(*) AS n FROM {tables[kind]} "
                f"WHERE star_system_id IN ({_marks(len(chunk))}) GROUP BY star_system_id", tuple(chunk)).fetchall()
            for row in rows:
                counts.setdefault(row["sid"], {})[kind] = row["n"]
    return counts


def _bodies_of(conn, system_id, system, wanted):
    """The wanted bodies of one system, in `BODY_KINDS` order, then id."""
    rows = []
    parent = {"ref": _ref("system", system_id), "kind": "system", "name": system["name"]}
    if "star" in wanted:
        for row in conn.execute("SELECT id, name FROM stars WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchall():
            rows.append(("star", row["id"], row["name"], parent))
    planet_names = {}
    if "planet" in wanted or "moon" in wanted:
        for row in conn.execute("SELECT id, name FROM planets WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchall():
            planet_names[row["id"]] = row["name"]
            if "planet" in wanted:
                rows.append(("planet", row["id"], row["name"], parent))
    if "moon" in wanted:
        for row in conn.execute("SELECT id, name, planet_id FROM moons WHERE star_system_id = ? ORDER BY id",
                                (system_id,)).fetchall():
            rows.append(("moon", row["id"], row["name"],
                         {"ref": _ref("planet", row["planet_id"]), "kind": "planet",
                          "name": planet_names.get(row["planet_id"])}))
    if "belt" in wanted:
        for row in conn.execute("SELECT id, orbital_index FROM asteroid_belts WHERE star_system_id = ? ORDER BY id",
                                (system_id,)).fetchall():
            rows.append(("belt", row["id"], f"Asteroid belt {row['orbital_index'] + 1}", parent))
    if "comet" in wanted:
        for row in conn.execute("SELECT id, name FROM comets WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchall():
            rows.append(("comet", row["id"], row["name"], parent))
    return rows


def _facilities(conn, center, distance, margin, systems):
    """Facilities within reach: `(distance, id, name, parent, point)`."""
    found = []
    reach = distance + margin
    for chunk in _chunks(systems):
        rows = conn.execute(
            f"SELECT id, name, star_system_id FROM facilities WHERE star_system_id IN ({_marks(len(chunk))})",
            tuple(chunk)).fetchall()
        for row in rows:
            system = systems[row["star_system_id"]]
            found.append((system["distance"], row["id"], row["name"],
                          {"ref": _ref("system", row["star_system_id"]), "kind": "system", "name": system["name"]},
                          system["point"]))
    rows = conn.execute(
        """
        SELECT f.id, f.name, f.center_x_pc, f.center_y_pc, f.center_z_pc, f.sector_id, sec.name AS sector_name
        FROM facilities f LEFT JOIN sectors sec ON sec.id = f.sector_id
        WHERE f.star_system_id IS NULL AND f.center_x_pc BETWEEN ? AND ? AND f.center_y_pc BETWEEN ? AND ?
          AND f.center_z_pc BETWEEN ? AND ?
        """,
        (center[0] - reach, center[0] + reach, center[1] - reach, center[1] + reach,
         center[2] - reach, center[2] + reach)).fetchall()
    for row in rows:
        point = (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])
        away = math.dist(center, point)
        if away <= distance:
            parent = ({"ref": _ref("sector", row["sector_id"]), "kind": "sector", "name": row["sector_name"]}
                      if row["sector_id"] is not None else None)
            found.append((away, row["id"], row["name"], parent, point))
    return found


def _phenomena(conn, center, distance, kinds):
    """Placed phenomena within `distance` of `center`: `(distance, kind, id, name, parent, point)`."""
    found = []
    widest = max(query._widest_placed_phenomenon_radius_ly(conn), 0.0)
    margin = 0.0 if widest <= 0 else ly_to_pc(widest)  #  an extended object's center may lie outside
    reach = distance + margin
    box = (center[0], center[1], center[2], reach)
    for row in query._placed_phenomenon_rows(conn, bbox=box):
        if row["type"] not in kinds:
            continue
        point = (row["x"], row["y"], row["z"])
        away = math.dist(center, point)
        if away <= distance:
            found.append((away, row["type"], row["id"], row["name"], None, point))
    return found


def objects_within(conn, place, distance_pc, kinds=None, limit=DEFAULT_LIMIT, offset=0):
    """
    Everything generated within `distance_pc` of a place, nearest first.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        place (dict): `place_from_reference` or `place_from_point`.
        distance_pc (float): The search distance, 0 up to `MAX_DISTANCE_PC`.
        kinds (iterable, optional): Only these `SEARCH_KINDS`; all by default.
        limit (int): Rows per page, at most `MAX_LIMIT`.
        offset (int): Rows to skip.

    Returns:
        dict: `place` (`ref`, `kind`, `name`, `point_pc`), `distance_pc`,
            `total` (rows over all pages), `by_kind` (`{kind: count}`),
            `limit`, `offset`, `sectors_in_range` and `sectors_generated`
            (the sectors the sphere reaches and how many of them hold
            generated stars), and `rows`: each `ref`, `kind`, `id`, `name`,
            `parent` (`ref`, `kind`, `name`, or `None`) and `distance_pc`.
            The place itself is left out when it is a system or a phenomenon (a
            planet or moon is not removed from its own system's bodies).

    Raises:
        NearError: For a distance, kind, or paging value that cannot be used.
    """
    distance = _check_distance(distance_pc)
    wanted = _check_kinds(kinds)
    limit, offset = int(limit), int(offset)
    if limit < 1 or limit > MAX_LIMIT or offset < 0:
        raise NearError(f"limit is 1 to {MAX_LIMIT} and offset is 0 or more")
    center = place["point_pc"]
    shape = get_galaxy_shape(conn)
    edge = float(shape.edge_pc) if shape is not None else float(tuning.DEFAULT_SECTOR_EDGE_PC)
    margin = edge * 1.5  # a system lies inside its sector's cell; allow for the cell's slack

    systems = _systems_in_reach(conn, center, distance, margin)
    self_system = place["id"] if place["kind"] == "system" else None
    body_wanted = wanted & set(BODY_KINDS)
    counts = _body_counts(conn, list(systems), body_wanted) if body_wanted else {}

    # Everything that sorts on its own: (distance, rank, id, ...). A system's
    # bodies sort as one group right after the system.
    entries = []
    if "system" in wanted:
        for system_id, system in systems.items():
            if system_id != self_system:
                entries.append((system["distance"], 0, system_id, "system", system_id))
    for system_id in counts:
        entries.append((systems[system_id]["distance"], 1, system_id, "bodies", system_id))
    if FACILITY_KIND in wanted:
        for away, facility_id, name, parent, point in _facilities(conn, center, distance, margin, systems):
            entries.append((away, 2, facility_id, FACILITY_KIND, (name, parent)))
    phenomena = _phenomena(conn, center, distance, wanted)
    for away, kind, object_id, name, _parent, _point in phenomena:
        if place["kind"] == kind and place["id"] == object_id:
            continue
        entries.append((away, 3 + object_ref.PHENOMENON_KINDS.index(kind), object_id, "phenomenon",
                        (kind, name)))
    entries.sort(key=lambda entry: entry[:3])

    by_kind = {}
    sizes = []
    for entry in entries:
        if entry[3] == "bodies":
            here = counts[entry[2]]
            size = sum(here.values())
            for kind, count in here.items():
                by_kind[kind] = by_kind.get(kind, 0) + count
        else:
            size = 1
            kind = entry[3] if entry[3] != "phenomenon" else entry[4][0]
            by_kind[kind] = by_kind.get(kind, 0) + 1
        sizes.append(size)
    total = sum(sizes)

    rows = []
    position = 0
    for entry, size in zip(entries, sizes):
        if position + size <= offset:
            position += size
            continue
        if position >= offset + limit:
            break
        away, _rank, object_id, what, extra = entry
        if what == "system":
            system = systems[object_id]
            made = [("system", object_id, system["name"],
                     {"ref": _ref("sector", system["sector_id"]), "kind": "sector", "name": system["sector_name"]})]
        elif what == "bodies":
            made = _bodies_of(conn, object_id, systems[object_id], body_wanted)
        elif what == FACILITY_KIND:
            made = [(FACILITY_KIND, object_id, extra[0], extra[1])]
        else:
            made = [(extra[0], object_id, extra[1], None)]
        for kind, row_id, name, parent in made:
            if offset <= position < offset + limit:
                rows.append({"ref": _ref(kind, row_id), "kind": kind, "id": row_id, "name": name,
                             "parent": parent, "distance_pc": away})
            position += 1
    reached, generated = _sector_counts(conn, center, distance, margin)
    return {
        "place": {"ref": place["ref"], "kind": place["kind"], "name": place["name"],
                  "point_pc": list(center)},
        "distance_pc": distance, "total": total, "by_kind": by_kind, "limit": limit, "offset": offset,
        "sectors_in_range": reached, "sectors_generated": generated, "rows": rows,
    }
