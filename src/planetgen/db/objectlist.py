# planetgen.db.objectlist

"""
The list of every kind of object, read from the object index (UX.141, DB.79, PERF.140).

One keyset page of `object_index` for a kind (a system, star, planet, moon, belt, comet, facility or any phenomenon
type): sorted by name, class, mass or radius with the object ID as the tie-break, filtered by a name prefix, a class and
a sector. Every sort and filter is served by an index of the table (`idx_object_index_name`, `_class`, `_mass`,
`_radius`, `_sector`), a page never uses `OFFSET`, and the next page starts after the last row's object ID. The index
holds what the builder last wrote (`planetgen warm object-index`); a kind it has not been built for lists nothing.
"""

from planetgen.db import objectindex
from planetgen.db import query
from planetgen.galaxy import object_uid

SORTS = {
    "name": ("name",),
    "class": ("object_class", "name"),
    "mass": ("mass_kg",),
    "radius": ("radius_km",),
}
"""dict: Sort key -> the index columns it orders by (the object ID follows them): `idx_object_index_name`, `_class`,
`_mass` and `_radius`. A list sorted by class, mass or radius holds the objects that have one."""

KINDS = objectindex.KINDS
"""tuple: The kinds the lists hold."""

MAX_LIMIT = 200
"""int: The most rows a page holds."""

DEFAULT_LIMIT = 50
"""int: The rows of a page when none is asked for."""

_COLUMNS = ("kind, uid, source_id, name, sector_id, system_id, parent_uid, object_class, mass_kg, radius_km, placed")


def _printed(raw):
    return None if raw is None else object_uid.format_id(object_uid.from_bytes(raw))


def _item(row):
    return {"kind": row["kind"], "id": _printed(row["uid"]), "source_id": row["source_id"], "name": row["name"],
            "sector_id": row["sector_id"], "system_id": row["system_id"], "parent_id": _printed(row["parent_uid"]),
            "object_class": row["object_class"], "mass_kg": row["mass_kg"], "radius_km": row["radius_km"],
            "placed": bool(row["placed"])}


HABITABILITY_KINDS = {"planet": "planets", "moon": "moons"}
"""dict: The kinds whose items carry their PHI-6 rating (UX.100), and the table that stores it."""

PHYSICS_COLUMNS = {
    "star": ("temperature_k", "luminosity_w", "age_gy"),
    "planet": ("distance_km", "eccentricity", "surface_temperature_k"),
    "moon": ("distance_km", "eccentricity", "surface_temperature_k"),
}
"""dict: The stored physics values each kind's items carry beside the index's own (UX.187: a star's temperature, luminosity
and age; a planet's or moon's orbit size, eccentricity and surface temperature), read from the kind's table. The lists show
them; they are not sorts (the index holds no column to order by)."""

_TABLES = {"star": "stars", "planet": "planets", "moon": "moons"}

HABITABILITY_COLUMNS = query._tier_columns()
"""tuple: The stored overall score and the six domain tiers (0 Blue to 3 Red) a planet or moon item carries."""


def _with_table_values(conn, kind, items):
    """Adds the kind's table values to each item (`None` for a row the table lacks), read by `source_id` in one query for
    the page: a planet's or moon's `HABITABILITY_COLUMNS` (UX.100, `None` for a body that was not scored) and the
    `PHYSICS_COLUMNS` of a star, planet or moon (UX.187)."""
    table = _TABLES.get(kind)
    columns = (*(HABITABILITY_COLUMNS if kind in HABITABILITY_KINDS else ()), *PHYSICS_COLUMNS.get(kind, ()))
    if table is None or not columns or not items:
        return items
    ids = [item["source_id"] for item in items]
    rows = conn.execute(f"SELECT id, {', '.join(columns)} FROM {table} WHERE id IN ({','.join('?' * len(ids))})", ids).fetchall()
    found = {row["id"]: row for row in rows}
    for item in items:
        row = found.get(item["source_id"])
        item.update({column: row[column] if row is not None else None for column in columns})
    return items


def _like_prefix(text):
    """`text` as the pattern of a `LIKE` that matches names starting with it."""
    return text.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_") + "%"


def _where(kind, sort, name, object_class, sector_id, start=None, descending=False):
    where, params = ["kind = ?"], [kind]
    if start:
        where.append("name <= ?" if descending else "name >= ?")
        params.append(start)
    if name:
        where.append("name LIKE ?")
        params.append(_like_prefix(name))
    if object_class:
        where.append("object_class = ?")
        params.append(object_class)
    if sector_id is not None:
        where.append("sector_id = ?")
        params.append(int(sector_id))
    for column in SORTS[sort]:
        if column in ("object_class", "mass_kg", "radius_km"):
            where.append(f"{column} IS NOT NULL")
    return where, params


def count_objects(conn, kind, name=None, object_class=None, sector_id=None, sort="name"):
    """How many objects of `kind` the filters leave, counted up to `query.COUNT_FALLBACK_CAP` (a `CappedCount` when
    there are more)."""
    where, params = _where(kind, sort, name, object_class, sector_id)
    row = conn.execute(f"SELECT COUNT(*) AS n FROM (SELECT 1 FROM object_index WHERE {' AND '.join(where)} LIMIT ?) c",
                       (*params, query.COUNT_FALLBACK_CAP)).fetchone()
    return query._capped(int(row["n"]))


def list_objects(conn, kind, limit=DEFAULT_LIMIT, sort="name", descending=False, after=None, name=None,
                 object_class=None, sector_id=None, start=None):
    """
    One keyset page of the objects of `kind`.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        kind (str): One of `KINDS`.
        limit (int): The rows of the page, at most `MAX_LIMIT`.
        sort (str): A key of `SORTS`; ties fall back to the object ID, in the same direction.
        descending (bool): Reverse the sort.
        after (str, optional): The printed object ID of the last row of the previous page: the page starts after it.
            An object that is gone (or is not of this kind) starts the list from the top.
        name (str, optional): Only objects whose name starts with this (case is not significant).
        object_class (str, optional): Only objects of exactly this class (`object_class`: a star's type, a planet's
            class, a phenomenon's descriptor).
        sector_id (int, optional): Only objects in this sector.
        start (str, optional): Begin the list at the first name from this one on in the name order (the jump box of
            the list pages); only with `sort="name"`. When no name lies at or past it the last page is answered with
            `past_end` true. `total` still counts the whole list.

    Returns:
        dict: `items` (as `_item`; a planet's or moon's also with `phi4` and the six `tier_<domain>` ratings, UX.100; a star's with its `temperature_k`, `luminosity_w` and `age_gy`, a planet's or moon's with its `distance_km`, `eccentricity` and `surface_temperature_k`, UX.187), `next` (the printed ID to pass as `after`, `None` at the end), `total` (a
        `CappedCount` past `query.COUNT_FALLBACK_CAP`) and, with a `start`, `past_end`.

    Raises:
        ValueError: An unknown kind or sort, a `limit` under 1, an `after` that is not an object ID, or a `start`
            under another sort.
    """
    if kind not in KINDS:
        raise ValueError(f"kind must be one of: {', '.join(KINDS)}")
    if sort not in SORTS:
        raise ValueError(f"sort must be one of: {', '.join(SORTS)}")
    if int(limit) < 1:
        raise ValueError("limit must be at least 1")
    if start and sort != "name":
        raise ValueError("a start needs sort=name")
    limit = min(int(limit), MAX_LIMIT)
    columns = list(SORTS[sort])
    where, params = _where(kind, sort, name, object_class, sector_id, start, descending)
    direction = "DESC" if descending else "ASC"
    if after:
        uid = object_uid.to_bytes(object_uid.parse_id(after))
        previous = conn.execute(f"SELECT {', '.join(columns)} FROM object_index WHERE kind = ? AND uid = ?",
                                (kind, uid)).fetchone()
        if previous is not None:
            condition, keyset = query._keyset_condition(columns, [previous[c] for c in columns], descending, "uid", uid)
            where.append(condition)
            params += keyset
    order = ", ".join(f"{column} {direction}" for column in [*columns, "uid"])
    rows = conn.execute(f"SELECT {_COLUMNS} FROM object_index WHERE {' AND '.join(where)} ORDER BY {order} LIMIT ?",
                        (*params, limit + 1)).fetchall()
    total = count_objects(conn, kind, name, object_class, sector_id, sort)
    if start and not rows and not after:
        last = list_objects(conn, kind, limit, sort, not descending, None, name, object_class, sector_id)
        last["items"].reverse()
        return {**last, "next": None, "total": total, "past_end": True}
    more = len(rows) > limit
    items = _with_table_values(conn, kind, [_item(row) for row in rows[:limit]])
    page = {"items": items, "next": items[-1]["id"] if more and items else None, "total": total}
    if start:
        page["past_end"] = False
    return page
