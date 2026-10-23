# planetgen/db/containers.py

"""
The galaxy's containers (DB.85): the stored globular clusters and superbubbles, which of them hold a point, and what a
container holds.

`galaxy_containers` has one row per container (`generation.containers` builds them from the galaxy seed at a plan). The
reverse of `cloudcontents`: that lists what sits in a nebula or supernova remnant by the `inside_*` columns, while a
container has no column on the objects it holds, so containment is geometric: an object is in a container when its
position lies within the container's extent of its centre. The UX.144 'Inside:' line asks `holding_object` and a
container's page asks `contents`.

Positions are read by the cell key (`db/region.py`), so a query costs the rows near the container and not the size of
the galaxy. The SQL test of a sphere is a billionth loose (`region.SLACK`), so a count of objects lying exactly on a
container's surface can differ by that.
"""

from planetgen.db import cloudcontents, region
from planetgen.galaxy import object_uid
from planetgen.util import detmath

KINDS = ("globular_cluster", "superbubble")
"""tuple: The kinds of container (the table's CHECK)."""

COLUMNS = ("uid", "kind", "name", "component_id", "bubble_col", "bubble_row", "bubble_n", "center_x_pc", "center_y_pc",
           "center_z_pc", "extent_pc", "shell_index", "band_index", "wedge_index")
"""tuple: The columns `replace` writes and `get` reads (the table's, less `id` and `created_at`)."""

DEFAULT_LIMIT = cloudcontents.DEFAULT_LIMIT

_FREE_POSITION = ("p.center_x_pc", "p.center_y_pc", "p.center_z_pc")
_SYSTEM_POSITION = ("p.galaxy_x_pc", "p.galaxy_y_pc", "p.galaxy_z_pc")

HELD_TABLES = tuple((kind, table, extra, _SYSTEM_POSITION if table == "star_systems" else _FREE_POSITION)
                    for kind, table, extra in cloudcontents.CONTENT_TABLES if table != "nebulae") + (
    ("nebula", "nebulae", "", _FREE_POSITION),)
"""tuple: `(kind, table, extra condition, position expressions)` for each table whose rows a container can hold, in
the order `cloudcontents.CONTENT_TABLES` lists them (a nebula last: it is measured by its centre like the rest)."""


def replace(conn, containers):
    """
    Replaces every stored container with `containers` (a plan draws the whole set again from the seed). Does not
    commit.

    Args:
        conn (Connection): The connection to write through.
        containers (iterable[generation.containers.Container]): The rows to write.
    """
    conn.execute("DELETE FROM galaxy_containers")
    rows = []
    for c in containers:
        bubble = c.bubble or (None, None, None)
        rows.append((c.uid, c.kind, c.name, c.component_id, *bubble, *c.center_pc, c.extent_pc, *c.birth))
    if rows:
        conn.executemany(
            f"INSERT INTO galaxy_containers ({', '.join(COLUMNS)}) VALUES ({', '.join('?' * len(COLUMNS))})", rows)


def _public(row):
    """A row as the JSON-able dict every function here returns."""
    out = {"id": row["id"], "uid": object_uid.format_id(object_uid.from_bytes(bytes(row["uid"]))), "kind": row["kind"],
           "name": row["name"], "component_id": row["component_id"],
           "center_pc": [row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]],
           "extent_pc": row["extent_pc"],
           "region": {"shell": row["shell_index"], "band": row["band_index"], "wedge": row["wedge_index"]}}
    if row["bubble_col"] is not None:
        out["bubble"] = [row["bubble_col"], row["bubble_row"], row["bubble_n"]]
    return out


def get(conn, container_id=None, uid=None):
    """
    One container by row id or by printed object ID (`FE81000A2B-0000005-000`), or `None`.

    Raises:
        ValueError: For a printed ID that is not an object ID, or when neither is given.
    """
    if container_id is not None:
        row = conn.execute("SELECT * FROM galaxy_containers WHERE id = ?", (container_id,)).fetchone()
    elif uid is not None:
        value = uid if isinstance(uid, int) else object_uid.parse_id(uid)
        row = conn.execute("SELECT * FROM galaxy_containers WHERE uid = ?", (object_uid.to_bytes(value),)).fetchone()
    else:
        raise ValueError("a container is asked for by its row id or its object ID")
    return None if row is None else _public(row)


def listing(conn, kind=None):
    """Every container, or those of one `kind`, by kind then id."""
    if kind is not None and kind not in KINDS:
        raise ValueError(f"a container is one of {', '.join(KINDS)}, not {kind!r}")
    sql = "SELECT * FROM galaxy_containers" + (" WHERE kind = ?" if kind else "") + " ORDER BY kind, id"
    return [_public(row) for row in conn.execute(sql, (kind,) if kind else ()).fetchall()]


def holding(conn, point_pc, kinds=None):
    """
    The containers whose sphere holds a point, deepest first: the one the point is nearest the middle of, as a share
    of the container's extent, comes first (superbubbles overlap, so there can be several).

    Args:
        conn (Connection): An open connection.
        point_pc (sequence): `(x, y, z)`, galaxy-frame parsecs.
        kinds (iterable, optional): Only these `KINDS`.

    Returns:
        list[dict]: Each container as `get` gives it, plus `depth` (the distance from its centre over its extent,
            under 1 inside) and `distance_pc`.
    """
    wanted = KINDS if kinds is None else tuple(kinds)
    unknown = [kind for kind in wanted if kind not in KINDS]
    if unknown:
        raise ValueError(f"a container is one of {', '.join(KINDS)}, not {unknown[0]!r}")
    x, y, z = (float(c) for c in point_pc)
    marks = ", ".join("?" * len(wanted))
    rows = conn.execute(
        f"SELECT * FROM galaxy_containers WHERE kind IN ({marks})"
        " AND ABS(center_x_pc - ?) <= extent_pc AND ABS(center_y_pc - ?) <= extent_pc"
        " AND ABS(center_z_pc - ?) <= extent_pc ORDER BY id", (*wanted, x, y, z)).fetchall()
    found = []
    for row in rows:
        away = detmath.dist((row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]), (x, y, z))
        if away <= row["extent_pc"]:
            found.append({**_public(row), "depth": away / row["extent_pc"], "distance_pc": away})
    found.sort(key=lambda item: (item["depth"], item["id"]))
    return found


def position_of(conn, kind, object_id):
    """
    The galaxy-frame position, parsecs, of a held object (`HELD_TABLES` kind and row id), or `None` when there is no
    such row or it has no position.

    Raises:
        ValueError: For a kind a container cannot hold.
    """
    for held, table, _extra, position in HELD_TABLES:
        if held == kind:
            columns = ", ".join(f"{expression} AS c{index}" for index, expression in enumerate(position))
            row = conn.execute(f"SELECT {columns} FROM {table} p WHERE p.id = ?", (object_id,)).fetchone()
            if row is None or row["c0"] is None:
                return None
            return (row["c0"], row["c1"], row["c2"])
    raise ValueError(f"a container holds {', '.join(h[0] for h in HELD_TABLES)}, not {kind!r}")


def holding_object(conn, kind, object_id, kinds=None):
    """
    The containers an object is inside: `holding` at its position (a system's, a phenomenon's centre). Empty for an
    object with no row or no position.
    """
    point = position_of(conn, kind, object_id)
    return [] if point is None else holding(conn, point, kinds)


def contents(conn, container_id, limit=DEFAULT_LIMIT):
    """
    What sits inside a container: every system, nebula, free-standing black hole, neutron star, asteroid field, rogue
    planet and interstellar comet whose position lies within its extent.

    Args:
        conn (Connection): An open connection.
        container_id (int): The container's row id.
        limit (int): The most objects to list (by kind, then name, then id); the counts are complete.

    Returns:
        dict: `container` (as `get` gives it), `counts` (`{kind: how many}` for each kind with any), `items`
            (`[{"kind", "id", "name"}]`, at most `limit`) and `total`.

    Raises:
        ValueError: If there is no such container.
    """
    row = conn.execute("SELECT * FROM galaxy_containers WHERE id = ?", (container_id,)).fetchone()
    if row is None:
        raise ValueError(f"no container {container_id}")
    center = (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])
    counts = {}
    items = []
    for kind, table, extra, position in HELD_TABLES:
        where = extra.removeprefix(" AND ") if extra else ""
        sql, runs = region.region_statement(table, "p.id, p.name", position, center, row["extent_pc"], where=where)
        found = conn.execute(f"SELECT COUNT(*) AS n FROM ({sql}) AS held").fetchone()["n"]
        if not found:
            continue
        counts[kind] = found
        room = limit - len(items)
        if room > 0:
            listed = conn.execute(f"SELECT id, name FROM ({sql}) AS held ORDER BY name, id LIMIT {int(room)}").fetchall()
            items.extend({"kind": kind, "id": item["id"], "name": item["name"]} for item in listed)
    return {"container": _public(row), "counts": counts, "items": items, "total": sum(counts.values())}
