# planetgen/db/children.py

"""
What belongs to an object, one keyset page at a time (UX.169): the feed of the summary trees on the sector page, the
system page and every object's own page, and `GET /api/objects/<ref>/children`.

The tree is the parent links the kind tables carry (`parent_type`, `parent_id`, DB.73) and the system they sit in:

- a sector holds its systems, then the phenomena placed in it;
- a system holds its stars, then the planets, belts and comets that orbit the whole system (a close pair's
  circumbinary bodies; the head of a single star's or of a wide pair's own bodies is the star);
- a star holds its planets and belts in their order out from it (the order of the system's list), then its comets;
- a planet holds its moons, and so does a free (rogue) planet;
- a moon, belt, comet and the rest hold nothing.

A page never uses `OFFSET`: the children are read as an ordered run of streams (each a union of tables ordered by
`orbital_index`, then the table, then the row id) and the next page starts after the last row's `(stream, order,
table, id)`, the opaque token `next` carries. Each item has a short `summary` of what its row says about it, `counts`
of everything below it by kind (a system's planets, a star's moons) and `child_count`, how many it holds directly, so a
row can say "3 planets, 5 moons" and show its expand arrow without a second request.

Reads only; row ids throughout (the API turns them into Object IDs, `api/idpaths.py`).
"""

from planetgen.db import query
from planetgen.galaxy import objectref

DEFAULT_LIMIT = 50
"""int: Children a page holds when none is asked for."""

MAX_LIMIT = 200
"""int: The most children a page holds."""

class PageError(ValueError):
    """A `limit` or `after` that is not one the page can use (the API answers 400, not 404)."""


PARENT_KINDS = ("sector", "system", "star", "planet", "rogue_planet")
"""tuple: The kinds that can hold children (any other kind answers an empty list)."""

_TABLES = {"system": "star_systems", "star": "stars", "planet": "planets", "moon": "moons", "belt": "asteroid_belts",
           "comet": "comets", "sector": "sectors", "rogue_planet": "rogue_planets"}

_SUMMARY = {
    "system": ("name", "is_binary", "binary_configuration", "quadrant"),
    "star": ("name", "role", "star_type", "yerkes_class", "mass_kg", "radius_km", "temperature_k", "luminosity_w"),
    "planet": ("name", "planet_class", "body_type", "distance_km", "period_years", "radius_km", "mass_kg", "gravity_g",
               *query._tier_columns()),
    "moon": ("name", "planet_class", "body_type", "distance_km", "period_years", "radius_km", "mass_kg", "gravity_g",
             *query._tier_columns()),
    "belt": ("orbital_index", "density", "distance_km", "lower_limit_km", "upper_limit_km"),
    "comet": ("name", "orbit_type", "period_class", "perihelion_distance_km", "orbital_period_years", "is_active"),
}
"""dict: The columns of its own row that a kind's `summary` carries."""

_PHENOMENA = tuple((kind, table) for kind, table in sorted(query._PHENOMENON_TYPE_TO_TABLE.items())
                   if kind in objectref.PHENOMENON_KINDS)


class _Member:
    """One table of a stream: the rows `where` selects, ordered by `order` (a column or constant), then by `sub`."""

    def __init__(self, kind, table, where, params, order="0", sub=0):
        self.kind, self.table, self.where, self.params, self.order, self.sub = kind, table, where, tuple(params), order, sub


def _parent(kind, table, parent_type, parent_id, order="orbital_index", sub=0):
    return _Member(kind, table, "parent_type = ? AND parent_id = ?", (parent_type, parent_id), order, sub)


def _streams(kind, object_id):
    """The ordered streams of a parent: each a list of `_Member`s whose rows are merged by `(order, sub, id)`."""
    if kind == "sector":
        return [[_Member("system", "star_systems", "sector_id = ?", (object_id,))],
                [_Member(phenomenon, table, "sector_id = ?", (object_id,), order=str(number))
                 for number, (phenomenon, table) in enumerate(_PHENOMENA)]]
    if kind == "system":
        return [[_Member("star", "stars", "star_system_id = ?", (object_id,), order="COALESCE(slot, 0)")],
                [_parent("planet", "planets", "system", object_id), _parent("belt", "asteroid_belts", "system", object_id,
                                                                            sub=1)],
                [_parent("comet", "comets", "system", object_id, order="0")]]
    if kind == "star":
        return [[_parent("planet", "planets", "star", object_id), _parent("belt", "asteroid_belts", "star", object_id, sub=1)],
                [_parent("comet", "comets", "star", object_id, order="0")]]
    if kind in ("planet", "rogue_planet"):
        return [[_parent("moon", "moons", kind, object_id)]]
    return []


def _token(stream, row):
    return f"{stream}.{row['ord']}.{row['sub']}.{row['id']}"


def _parse_token(after):
    try:
        stream, order, sub, row_id = (int(part) for part in str(after).split("."))
    except ValueError:
        raise PageError(f"not a children page token: {after!r}") from None
    if min(stream, sub, row_id) < 0:
        raise PageError(f"not a children page token: {after!r}")
    return stream, (order, sub, row_id)


def _stream_page(conn, members, low, count):
    selects = " UNION ALL ".join(
        f"SELECT '{member.kind}' AS kind, id, {member.order} AS ord, {member.sub} AS sub FROM {member.table} "
        f"WHERE {member.where}" for member in members)
    params = [value for member in members for value in member.params]
    bound = ""
    if low is not None:
        bound = " WHERE (ord, sub, id) > (?, ?, ?)"
        params += list(low)
    return conn.execute(
        f"SELECT kind, id, ord, sub FROM ({selects}) AS members{bound} ORDER BY ord, sub, id LIMIT {int(count)}",
        params).fetchall()


def _count_rows(conn, sql, ids):
    marks = ", ".join("?" for _ in ids)
    return conn.execute(sql.format(marks=marks), list(ids)).fetchall()


def _direct_counts(conn, kind, ids):
    """`{id: how many children it holds}` for rows of a parent `kind`."""
    counts = {row_id: 0 for row_id in ids}
    if not ids or kind not in PARENT_KINDS:
        return counts
    if kind == "sector":
        for table, key in (("star_systems", "sector_id"), *((table, "sector_id") for _kind, table in _PHENOMENA)):
            for row in _count_rows(conn, f"SELECT {key} AS owner, COUNT(*) AS n FROM {table} WHERE {key} IN ({{marks}}) "
                                   f"GROUP BY {key}", ids):
                counts[row["owner"]] += row["n"]
        return counts
    if kind == "system":
        for row in _count_rows(conn, "SELECT star_system_id AS owner, COUNT(*) AS n FROM stars "
                                     "WHERE star_system_id IN ({marks}) GROUP BY star_system_id", ids):
            counts[row["owner"]] += row["n"]
    for table in {"system": ("planets", "asteroid_belts", "comets"), "star": ("planets", "asteroid_belts", "comets"),
                  "planet": ("moons",), "rogue_planet": ("moons",)}[kind]:
        for row in _count_rows(conn, f"SELECT parent_id AS owner, COUNT(*) AS n FROM {table} WHERE parent_type = '{kind}' "
                               f"AND parent_id IN ({{marks}}) GROUP BY parent_id", ids):
            counts[row["owner"]] += row["n"]
    return counts


def _below_counts(conn, kind, ids):
    """`{id: {kind: how many of it lie anywhere below}}` for systems, stars and planets (a sector's would be a whole
    scan; it has none)."""
    counts = {row_id: {} for row_id in ids}
    if not ids or kind not in ("system", "star", "planet"):
        return counts

    def add(sql, label):
        for row in _count_rows(conn, sql, ids):
            if row["n"]:
                counts[row["owner"]][label] = row["n"]

    if kind == "system":
        for label, table in (("star", "stars"), ("planet", "planets"), ("moon", "moons"), ("belt", "asteroid_belts"),
                             ("comet", "comets")):
            add(f"SELECT star_system_id AS owner, COUNT(*) AS n FROM {table} WHERE star_system_id IN ({{marks}}) "
                f"GROUP BY star_system_id", label)
    elif kind == "star":
        for label, table in (("planet", "planets"), ("belt", "asteroid_belts"), ("comet", "comets")):
            add(f"SELECT parent_id AS owner, COUNT(*) AS n FROM {table} WHERE parent_type = 'star' "
                f"AND parent_id IN ({{marks}}) GROUP BY parent_id", label)
        add("SELECT p.parent_id AS owner, COUNT(*) AS n FROM moons m JOIN planets p ON p.id = m.planet_id "
            "WHERE p.parent_type = 'star' AND p.parent_id IN ({marks}) GROUP BY p.parent_id", "moon")
    else:
        add("SELECT parent_id AS owner, COUNT(*) AS n FROM moons WHERE parent_type = 'planet' "
            "AND parent_id IN ({marks}) GROUP BY parent_id", "moon")
    return counts


def _fetch(conn, kind, ids):
    """`{id: the row's summary columns}` for rows of one child kind."""
    if kind in _SUMMARY:
        columns = ", ".join(("id", *_SUMMARY[kind]))
        rows = _count_rows(conn, f"SELECT {columns} FROM {_TABLES[kind]} WHERE id IN ({{marks}})", ids)
    else:
        table = query._PHENOMENON_TYPE_TO_TABLE[kind]
        rows = _count_rows(conn, f"SELECT id, name FROM {table} WHERE id IN ({{marks}})", ids)
    return {row["id"]: dict(row) for row in rows}


def _star_types(conn, system_ids):
    found = {row_id: [] for row_id in system_ids}
    for row in _count_rows(conn, "SELECT star_system_id, star_type FROM stars WHERE star_system_id IN ({marks}) "
                                 "ORDER BY COALESCE(slot, 0), id", system_ids):
        found[row["star_system_id"]].append(row["star_type"])
    return found


def _items(conn, picked):
    """The items of a page: `picked` is `[(kind, id)]` in order."""
    by_kind = {}
    for kind, row_id in picked:
        by_kind.setdefault(kind, []).append(row_id)
    rows, direct, below, star_types = {}, {}, {}, {}
    for kind, ids in by_kind.items():
        rows[kind] = _fetch(conn, kind, ids)
        direct[kind] = _direct_counts(conn, kind, ids)
        below[kind] = _below_counts(conn, kind, ids)
    if "system" in by_kind:
        star_types = _star_types(conn, by_kind["system"])
    items = []
    for kind, row_id in picked:
        row = rows[kind][row_id]
        name = f"Asteroid belt {row['orbital_index'] + 1}" if kind == "belt" else row.pop("name", None)
        summary = {key: value for key, value in row.items() if key not in ("id", "name")}
        if kind == "system":
            summary["is_binary"] = bool(summary["is_binary"])
            summary["star_types"] = star_types[row_id]
        if kind == "comet":
            summary["is_active"] = bool(summary["is_active"])
        items.append({"ref": objectref.format(kind, row_id), "kind": kind, "id": row_id, "name": name, "summary": summary,
                      "counts": below[kind][row_id], "child_count": direct[kind][row_id]})
    return items


def children(conn, kind, object_id, limit=DEFAULT_LIMIT, after=None):
    """
    One page of what belongs to an object (the module's docstring has the shape of the tree).

    Args:
        conn: An open, read-only connection.
        kind (str): The parent's kind.
        object_id (int): The parent's row id.
        limit (int): Children in the page, 1 to `MAX_LIMIT`.
        after (str, optional): The `next` of the previous page.

    Returns:
        dict: `ref`, `kind`, `total` (every child it holds), `items` (each `ref`, `kind`, `id`, `name`, `summary`,
            `counts`, `child_count`) and `next` (the token for the page after, `None` at the end).

    Raises:
        PageError: For a `limit` out of range or a bad `after`.
        ValueError: For an unknown kind or a missing parent row.
    """
    if kind not in objectref.KINDS:
        raise ValueError(f"unknown object kind: {kind!r}")
    if not 1 <= int(limit) <= MAX_LIMIT:
        raise PageError(f"limit must be 1 to {MAX_LIMIT}, got {limit!r}")
    table = _TABLES.get(kind) or query._PHENOMENON_TYPE_TO_TABLE.get(kind)
    if table is None or conn.execute(f"SELECT 1 FROM {table} WHERE id = ?", (object_id,)).fetchone() is None:
        raise ValueError(f"no {kind} with id {object_id}")
    result = {"ref": objectref.format(kind, object_id), "kind": kind, "total": 0, "items": [], "next": None}
    if kind not in PARENT_KINDS:
        return result
    result["total"] = _direct_counts(conn, kind, [object_id])[object_id]
    streams = _streams(kind, object_id)
    first, low = _parse_token(after) if after else (0, None)
    if first >= len(streams):
        raise PageError(f"not a children page token: {after!r}")
    picked = []
    for index in range(first, len(streams)):
        want = int(limit) + 1 - len(picked)
        if want <= 0:
            break
        picked += [(index, row) for row in _stream_page(conn, streams[index], low if index == first else None, want)]
    page = picked[:int(limit)]
    result["items"] = _items(conn, [(row["kind"], row["id"]) for _index, row in page])
    if len(picked) > len(page):
        result["next"] = _token(*page[-1])
    return result
