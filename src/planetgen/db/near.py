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

How it stays cheap (NAV.60): every positioned table is read by its cell key (`db/region.py`): the sphere becomes at
most 256 runs of keys (`galaxy/cellkey.py`), joined to the table's covering index `(cell_key, x, y, z)` with the exact
distance test in the same statement, so the cost follows the rows near the place, never the size of the galaxy, and
an object with no sector is found like any other. The `sectors` table is not read for the search. A system's bodies
are never loaded for the whole sphere: only a per-system count is, and the bodies of the systems that fall on the
requested page are read afterwards.

What has not been generated yet is counted, never listed and never filled: the pre-placed bright stars and the
scatter objects (`bright_stars` and `phenomenon_scatter` rows whose system or object is not built).

The dense core (NAV.61, design note section 6.5): a sphere's cost follows the rows inside it, and the core holds 15
systems per cubic parsec, so a 50 pc sphere there holds millions. Before any row is read, `_rows_in_reach` counts what
the sphere holds with an early stop (`region.capped_count`, 11 ms at most however full the sphere is). A sphere over
`SPHERE_LIMIT` is refused with the text that says to shrink the distance; one over `COUNT_CAP` systems is reported as
"`COUNT_CAP`+" in `systems_in_range`. `nearest` asks for the nearest K rows through a ladder of radii (`LADDER_PC`),
stopping at the first sphere that holds K, so it is cheap wherever the place is.
"""

import math

from planetgen.db import query, region
from planetgen.galaxy import objectref as object_ref
from planetgen.util import detmath

MAX_DISTANCE_PC = 50.0
"""float: The largest distance a search takes (Boss, NAV.43: about 8,000
sectors at 4 pc a side)."""

DEFAULT_LIMIT = 50
MAX_LIMIT = 200

SPHERE_LIMIT = 100_000
"""int: The most systems (and unbuilt objects) a sphere may hold for a search to run. About 70 ms of scan on the design
note's measurements; a larger sphere is refused (NAV.61, Core's call)."""

COUNT_CAP = 10_000
"""int: A count of the systems in a sphere is shown as `COUNT_CAP`+ above this (the early-stop count reads one past it)."""

LADDER_PC = (2.0, 4.0, 8.0, 16.0, 32.0)
"""tuple: The radii `nearest` tries, smallest first; the last rung is always the largest distance asked for."""

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


def _systems_in_reach(conn, center, distance):
    """Every placed system whose position is within `distance` of `center`:
    `{id: {...}}` with `point` and `distance`. Only the columns of the covering cell index are read (PERF.217): the name and
    sector of a system are not in it, and a sphere of 60,000 systems read them from the table row by row, three times the cost
    of the id list; `_system_labels` reads them for the systems of one page."""
    columns = "p.id, p.galaxy_x_pc, p.galaxy_y_pc, p.galaxy_z_pc"
    position = ("p.galaxy_x_pc", "p.galaxy_y_pc", "p.galaxy_z_pc")
    found = {}
    for row in region.region_rows(conn, "star_systems", columns, position, center, distance):
        point = (row["galaxy_x_pc"], row["galaxy_y_pc"], row["galaxy_z_pc"])
        away = detmath.dist(center, point)
        if away <= distance:
            found[row["id"]] = {"point": point, "distance": away}
    return found


def _system_labels(conn, system_ids):
    """`{id: (name, sector_id)}` of the systems of one page (a read by id, never a search)."""
    labels = {}
    for chunk in _chunks(sorted(system_ids)):
        for row in conn.execute(f"SELECT id, name, sector_id FROM star_systems WHERE id IN ({_marks(len(chunk))})",
                                tuple(chunk)).fetchall():
            labels[row["id"]] = (row["name"], row["sector_id"])
    return labels


def _star_counts(conn, system_ids):
    """`{id: stars}` of the systems given, the rows of the `stars` table each holds (a system with none is absent)."""
    counts = {}
    for chunk in _chunks(sorted(system_ids)):
        for row in conn.execute(f"SELECT star_system_id AS sid, COUNT(*) AS n FROM stars WHERE star_system_id IN ({_marks(len(chunk))}) "
                                "GROUP BY star_system_id", tuple(chunk)).fetchall():
            counts[row["sid"]] = row["n"]
    return counts


def _systems_with_stars(conn, systems, stars):
    """`systems` cut to those that hold `stars` (`(least, most)` of `util.star_count.parse_star_count`, `most` `None` for "or more")."""
    least, most = stars
    counts = _star_counts(conn, systems)
    return {system_id: system for system_id, system in systems.items()
            if counts.get(system_id, 0) >= least and (most is None or counts[system_id] <= most)}


def _sector_names(conn, sector_ids):
    """`{id: name}` of the sectors that label the rows of one page (a read of the sectors by id, never a search)."""
    names = {}
    for chunk in _chunks(sorted({sector_id for sector_id in sector_ids if sector_id is not None})):
        for row in conn.execute(f"SELECT id, name FROM sectors WHERE id IN ({_marks(len(chunk))})",
                                tuple(chunk)).fetchall():
            names[row["id"]] = row["name"]
    return names


def _unbuilt(conn, center, distance):
    """`{kind: count}` of what lies in the sphere and is not generated yet: pre-placed bright stars (`bright_star`)
    and the scatter objects not built (their phenomenon kind), read by the same statement on the scatter tables."""
    position = ("p.position_x_mpc / 1000.0", "p.position_y_mpc / 1000.0", "p.position_z_mpc / 1000.0")
    counts = {}
    runs = None
    for table, kind, built in (("bright_stars", None, "p.star_system_id IS NULL"),
                               ("phenomenon_scatter", "p.kind", "p.built_at IS NULL")):
        columns = ("'bright_star'" if kind is None else kind) + " AS kind, " + ", ".join(
            f"{axis} AS {name}" for axis, name in zip(position, ("x", "y", "z")))
        sql, runs = region.region_statement(table, columns, position, center, distance, built, runs)
        for row in conn.execute(sql).fetchall():
            if detmath.dist(center, (row["x"], row["y"], row["z"])) <= distance:
                counts[row["kind"]] = counts.get(row["kind"], 0) + 1
    return counts


def _rows_in_reach(conn, center, distance):
    """
    What the sphere holds, counted with an early stop: `(systems, unbuilt)`, each no more than `SPHERE_LIMIT + 1` (a
    count of `SPHERE_LIMIT + 1` means "more than the limit"). `unbuilt` is the bright stars and scatter objects not built.
    """
    cap = SPHERE_LIMIT + 1
    systems = region.capped_count(conn, "star_systems", ("p.galaxy_x_pc", "p.galaxy_y_pc", "p.galaxy_z_pc"),
                                  center, distance, cap)
    position = ("p.position_x_mpc / 1000.0", "p.position_y_mpc / 1000.0", "p.position_z_mpc / 1000.0")
    unbuilt = (region.capped_count(conn, "bright_stars", position, center, distance, cap,
                                   where="p.star_system_id IS NULL")
               + region.capped_count(conn, "phenomenon_scatter", position, center, distance, cap,
                                     where="p.built_at IS NULL"))
    return systems, unbuilt


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


def _bodies_of(conn, system_id, system_name, wanted):
    """The wanted bodies of one system, in `BODY_KINDS` order, then id."""
    rows = []
    parent = {"ref": _ref("system", system_id), "kind": "system", "name": system_name}
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


def _facilities(conn, center, distance, systems, standalone=True):
    """Facilities within reach: `(distance, id, name, parent, point)`; a standalone one's `parent` carries its
    `sector_id` until the page names it, and one in a system has the system's name left `None` until the page is made
    (`_system_labels`)."""
    found = []
    for chunk in _chunks(systems):
        rows = conn.execute(
            f"SELECT id, name, star_system_id FROM facilities WHERE star_system_id IN ({_marks(len(chunk))})",
            tuple(chunk)).fetchall()
        for row in rows:
            system = systems[row["star_system_id"]]
            found.append((system["distance"], row["id"], row["name"],
                          {"ref": _ref("system", row["star_system_id"]), "kind": "system", "name": None},
                          system["point"]))
    if not standalone:
        return found
    position = ("p.center_x_pc", "p.center_y_pc", "p.center_z_pc")
    for row in region.region_rows(conn, "facilities", "p.id, p.name, p.sector_id, p.center_x_pc, p.center_y_pc, p.center_z_pc",
                                  position, center, distance, where="p.star_system_id IS NULL"):
        point = (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])
        away = detmath.dist(center, point)
        if away <= distance:
            parent = {"ref": _ref("sector", row["sector_id"]), "kind": "sector", "name": None} \
                if row["sector_id"] is not None else None
            found.append((away, row["id"], row["name"], parent, point))
    return found


def _phenomena(conn, center, distance, kinds):
    """Placed phenomena within `distance` of `center`: `(distance, kind, id, name, parent, point)`. Measured to the
    object's position (a nebula's center), as NAV.43 states; the key finds it by that position."""
    found = []
    runs = None
    position = ("p.center_x_pc", "p.center_y_pc", "p.center_z_pc")
    for table, type_label, _descriptor, _radius in query._PHENOMENON_TABLES:
        if type_label not in kinds:
            continue
        sql, runs = region.region_statement(
            table, "p.id, p.name, p.center_x_pc, p.center_y_pc, p.center_z_pc", position, center, distance, runs=runs)
        for row in conn.execute(sql).fetchall():
            point = (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])
            away = detmath.dist(center, point)
            if away <= distance:
                found.append((away, type_label, row["id"], row["name"], None, point))
    return found


def objects_within(conn, place, distance_pc, kinds=None, limit=DEFAULT_LIMIT, offset=0, stars=None):
    """
    Everything generated within `distance_pc` of a place, nearest first.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        place (dict): `place_from_reference` or `place_from_point`.
        distance_pc (float): The search distance, 0 up to `MAX_DISTANCE_PC`.
        kinds (iterable, optional): Only these `SEARCH_KINDS`; all by default.
        limit (int): Rows per page, at most `MAX_LIMIT`.
        offset (int): Rows to skip.
        stars (tuple, optional): The stars filter `(least, most)` of `util.star_count.parse_star_count` (UX.239): only the systems that hold that
            many stars, and what lives in them. A phenomenon and a facility outside a system are not in one, so they are left out.

    Returns:
        dict: `place` (`ref`, `kind`, `name`, `point_pc`), `distance_pc`,
            `total` (rows over all pages), `by_kind` (`{kind: count}`),
            `limit`, `offset`, `not_generated` (`{kind: count}` of the pre-placed
            bright stars, `bright_star`, and scatter objects in the sphere that are not
            built yet; counted, never listed), `not_generated_total`, and `rows`: each `ref`, `kind`, `id`, `name`,
            `parent` (`ref`, `kind`, `name`, or `None`), `system_id` (the system
            it is in or is, else `None`), `distance_pc` and, on a `system` row, `star_count`.
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

    in_reach, unbuilt_in_reach = _rows_in_reach(conn, center, distance)
    if max(in_reach, unbuilt_in_reach) > SPHERE_LIMIT:
        raise NearError(
            f"{distance:g} pc around {place['name']} holds more than {SPHERE_LIMIT:,} systems, too many to list; "
            f"shrink the distance (a smaller sphere is much cheaper) or ask for the nearest few with nearest")

    systems = _systems_in_reach(conn, center, distance)
    if stars is not None:
        systems = _systems_with_stars(conn, systems, stars)
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
        for away, facility_id, name, parent, point in _facilities(conn, center, distance, systems, standalone=stars is None):
            entries.append((away, 2, facility_id, FACILITY_KIND, (name, parent)))
    phenomena = _phenomena(conn, center, distance, wanted) if stars is None else []
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

    on_page = set()                     # the systems the page names: its own, its bodies' and its facilities' parents
    position = 0
    for entry, size in zip(entries, sizes):
        if position + size > offset and position < offset + limit:
            if entry[3] in ("system", "bodies"):
                on_page.add(entry[2])
            elif entry[3] == FACILITY_KIND and entry[4][1] and entry[4][1]["kind"] == "system":
                on_page.add(object_ref.parse(entry[4][1]["ref"])[1])
        position += size
    labels = _system_labels(conn, on_page)
    star_counts = _star_counts(conn, on_page)

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
            name, sector_id = labels.get(object_id, (None, None))
            made = [("system", object_id, name,
                     {"ref": _ref("sector", sector_id), "kind": "sector", "name": None} if sector_id is not None else None)]
            home = object_id
        elif what == "bodies":
            made = _bodies_of(conn, object_id, labels.get(object_id, (None, None))[0], body_wanted)
            home = object_id
        elif what == FACILITY_KIND:
            parent = extra[1]
            if parent and parent["kind"] == "system":
                parent = dict(parent, name=labels.get(object_ref.parse(parent["ref"])[1], (None, None))[0])
            made = [(FACILITY_KIND, object_id, extra[0], parent)]
            home = (object_ref.parse(extra[1]["ref"])[1] if extra[1] and extra[1]["kind"] == "system" else None)
        else:
            made = [(extra[0], object_id, extra[1], None)]
            home = None
        for kind, row_id, name, parent in made:
            if offset <= position < offset + limit:
                row = {"ref": _ref(kind, row_id), "kind": kind, "id": row_id, "name": name,
                       "parent": parent, "system_id": home, "distance_pc": away}
                if kind == "system":
                    row["star_count"] = star_counts.get(row_id, 0)
                rows.append(row)
            position += 1
    in_sector = [row for row in rows if row["parent"] and row["parent"]["kind"] == "sector"]
    names = _sector_names(conn, [object_ref.parse(row["parent"]["ref"])[1] for row in in_sector])
    for row in in_sector:
        row["parent"]["name"] = names.get(object_ref.parse(row["parent"]["ref"])[1])
    unbuilt = _unbuilt(conn, center, distance)
    return {
        "place": {"ref": place["ref"], "kind": place["kind"], "name": place["name"],
                  "point_pc": list(center)},
        "distance_pc": distance, "total": total, "by_kind": by_kind, "limit": limit, "offset": offset,
        "not_generated": unbuilt, "not_generated_total": sum(unbuilt.values()), "rows": rows,
        "systems_in_range": min(in_reach, COUNT_CAP), "systems_in_range_capped": in_reach > COUNT_CAP,
    }


def nearest(conn, place, count, kinds=None, max_distance_pc=MAX_DISTANCE_PC, stars=None):
    """
    The nearest `count` rows to a place, found through a ladder of radii (NAV.61).

    Asks for `LADDER_PC` in turn (only the rungs under `max_distance_pc`, then `max_distance_pc` itself) and stops at the
    first sphere that holds `count` rows, so the cost follows the neighbourhood and not the distance to the next
    sparse star: at the core 200 systems lie within a few parsecs, at the Sun within about 8. Everything outside a
    sphere is farther than everything inside it, so the rows of the first sphere that holds `count` are the nearest.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        place (dict): `place_from_reference` or `place_from_point`.
        count (int): How many rows, 1 up to `MAX_LIMIT`.
        kinds (iterable, optional): Only these `SEARCH_KINDS`; all by default.
        max_distance_pc (float): Do not look farther than this, up to `MAX_DISTANCE_PC`.
        stars (tuple, optional): The stars filter of `objects_within`.

    Returns:
        dict: The result of `objects_within` for the sphere that answered, with at most `count` rows; its `distance_pc`
        is the rung that answered, and `total` may be below `count` when even the last rung holds fewer.

    Raises:
        NearError: For a count, distance or kind that cannot be used, or a rung the dense-core guard refuses.
    """
    count = int(count)
    if count < 1 or count > MAX_LIMIT:
        raise NearError(f"the number of nearest rows is 1 to {MAX_LIMIT}")
    farthest = _check_distance(max_distance_pc)
    rungs = [rung for rung in LADDER_PC if rung < farthest] + [farthest]
    for rung in rungs:
        result = objects_within(conn, place, rung, kinds=kinds, limit=count, stars=stars)
        if result["total"] >= count:
            break
    return result
