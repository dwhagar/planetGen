# planetgen/db/uncharted.py

"""
The objects no sector has stored yet, read where they are (MAP.400).

The object index (`db/objectindex.py`, PERF.140) holds one row per object that has an Object ID in a kind's own table. Two sets of
objects the map draws are not in it, and are not put there:

- the scattered phenomena (`phenomenon_scatter`, 9.7 million rows on the test galaxy): black holes, neutron stars, supernova remnants,
  planetary nebulae, hypervelocity stars and the nucleus quasar whose sector has not been built. The row's `built_at` is empty until a
  sector stores the object, and then the stored object is the one the index holds (the same Object ID, GEN.1254);
- the field clouds (`nebula_field`, 110 thousand rows): the molecular and dark clouds the galaxy's seed places. A cloud has no Object ID
  until a sector stores it (`is_stored`), so it is named by the key of its row, `field-<i>.<j>.<k>.<n>` (`db.field_clouds.cloud_id`).

An index of the scatter would take about 300 bytes a row, 3 GB, which is three times the cache budget (`cache.max_mb`, kept first by the
index) for rows nothing ever changes. So the lists, the search and Nearby read these tables beside the index: the scatter table by the
row id, the field by its key, the counts from `phenomenon_scatter_classes` (exact, kept by the scatter) and one count of the field.

A list shows the charted objects first, in the order asked for, and the uncharted ones after them, field clouds then scattered rows, each
in its table's order (`page`). A scattered row is named "Uncharted <kind> <shell>.<band>.<wedge>" like the Phenomena table names it and
links to the page of its scatter row (`/uncharted-phenomenon/<id>`); a field cloud has no page.
"""

from planetgen import tuning
from planetgen.db import query, region
from planetgen.galaxy import object_uid
from planetgen.physics import constants
from planetgen.util import detmath

TAIL_KINDS = ("nebula", "black_hole", "neutron_star", "supernova_remnant", "quasar", "hypervelocity_star")
"""tuple: The kinds that have uncharted objects. A hypervelocity star is a system once its sector is built, so its kind has no row in the
index and its list is this module's alone (`LIST_ONLY_KINDS`)."""

LIST_ONLY_KINDS = ("hypervelocity_star",)
"""tuple: The list kinds that have no row in `object_index`."""

CLOUD_PREFIX = "field-"
"""str: What a field cloud's printed ID starts with."""

FIELD_FAMILIES = {code: tuning.NEBULA_CLASSES[code]["family"] for code in ("M", "N")}
"""dict: The classes of `nebula_field` and the family each lists as (`object_class`)."""

_FIELD_COLUMNS = "cell_i, cell_j, cell_k, n, nebula_class, radius_ly, center_x_pc, center_y_pc, center_z_pc"

SOLAR_MASS_KG = constants.SOLAR_MASS_TO_KG
KM_PER_LY = constants.LY_TO_AU * (constants.AU_M / 1000.0)


def scatter_type(kind):
    """The list kind of a `phenomenon_scatter` row's `kind` (`black-hole` is `black_hole`; a hypervelocity star is its own list kind)."""
    return query._SCATTER_TYPES.get(kind, kind)


def cloud_id(row):
    """The printed ID of a field cloud from the key of its row."""
    return f"{CLOUD_PREFIX}{int(row['cell_i'])}.{int(row['cell_j'])}.{int(row['cell_k'])}.{int(row['n'])}"


def parse_cloud_id(text):
    """`(i, j, k, n)` of a printed field cloud ID, or `None` when `text` is not one."""
    if not isinstance(text, str) or not text.startswith(CLOUD_PREFIX):
        return None
    parts = text[len(CLOUD_PREFIX):].split(".")
    if len(parts) != 4:
        return None
    try:
        return tuple(int(part) for part in parts)
    except ValueError:
        return None


def _field_where(object_class=None, codes=None):
    """`(sql, params)` for the field clouds a class filter leaves (`object_class`, a family as the lists name it; `codes`, class letters as the
    search page names them), or `None` when none can match."""
    wanted = sorted(FIELD_FAMILIES)
    if object_class:
        wanted = [code for code in wanted if FIELD_FAMILIES[code] == object_class]
    if codes is not None:
        wanted = [code for code in wanted if code in codes]
    if not wanted:
        return None
    if len(wanted) == len(FIELD_FAMILIES):
        return "is_stored = 0", []
    return f"is_stored = 0 AND nebula_class IN ({', '.join('?' for _ in wanted)})", wanted


def _classes(conn, kind, object_class=None):
    """The scattered classes of `kind` that still have unbuilt rows, narrowed to one descriptor when `object_class` is given."""
    return query._scattered_classes(conn, types=(kind,), descriptors=(object_class,) if object_class else ())


def count(conn, kind, object_class=None):
    """How many uncharted objects of `kind` there are (those of `object_class` when it is given): zero for a kind with none."""
    if kind not in TAIL_KINDS:
        return 0
    total = sum(entry["unbuilt"] for entry in _classes(conn, kind, object_class))
    if kind == "nebula":
        field = _field_where(object_class)
        if field is not None:
            total += conn.execute(f"SELECT COUNT(*) AS n FROM nebula_field WHERE {field[0]}", field[1]).fetchone()["n"]
    return int(total)


def counts(conn):
    """`{kind: how many}` for every kind in `TAIL_KINDS`, in one pass over the class totals and one count of the field."""
    found = {kind: 0 for kind in TAIL_KINDS}
    for entry in query._scattered_classes(conn):
        if entry["type"] in found:
            found[entry["type"]] += entry["unbuilt"]
    found["nebula"] += int(conn.execute("SELECT COUNT(*) AS n FROM nebula_field WHERE is_stored = 0").fetchone()["n"])
    return found


def cursor(conn, kind, after):
    """
    Where a list that is past the index continues from, for the printed ID `after` of a row of the tail: `("field", (i, j, k, n))` or
    `("scatter", id)`; `None` when `after` is no uncharted object of `kind`.
    """
    if kind not in TAIL_KINDS or not after:
        return None
    key = parse_cloud_id(after)
    if key is not None:
        return ("field", key) if kind == "nebula" else None
    try:
        uid = object_uid.to_bytes(object_uid.parse_id(after))
    except ValueError:
        return None
    classes = query._scattered_classes(conn, types=(kind,))
    if not classes:
        return None
    row = conn.execute("SELECT id FROM phenomenon_scatter WHERE uid = ? AND built_at IS NULL", (uid,)).fetchone()
    return ("scatter", row["id"]) if row is not None else None


def _printed(raw):
    return None if raw is None else object_uid.format_id(object_uid.from_bytes(raw))


def _cloud_item(row):
    radius_ly = float(row["radius_ly"])
    return {"kind": "nebula", "id": cloud_id(row), "source_id": None, "name": f"Field cloud {cloud_id(row)[len(CLOUD_PREFIX):]}",
            "sector_id": None, "system_id": None, "parent_id": None, "object_class": FIELD_FAMILIES[row["nebula_class"]],
            "mass_kg": None, "radius_km": radius_ly * KM_PER_LY, "placed": True, "uncharted": True, "scatter_id": None}


def _scatter_item(row, kind, descriptor):
    mass = row["mass_solar"]
    return {"kind": kind, "id": _printed(row["uid"]), "source_id": None,
            "name": f"Uncharted {row['kind'].replace('-', ' ')} {row['shell_index']}.{row['band_index']}.{row['wedge_index']}",
            "sector_id": None, "system_id": None, "parent_id": None, "object_class": descriptor,
            "mass_kg": None if mass is None else float(mass) * SOLAR_MASS_KG, "radius_km": None, "placed": True,
            "uncharted": True, "scatter_id": row["id"]}


def page(conn, kind, limit, after=None, object_class=None):
    """
    Up to `limit` uncharted objects of `kind` that follow the cursor `after` (what `cursor` returned; `None` for the start): the field
    clouds in the order of their key, then the scattered rows in the order of their id. An item has the keys of `objectlist._item`, and
    `uncharted` (True) and `scatter_id` (the scatter row, `None` for a field cloud).
    """
    if kind not in TAIL_KINDS or limit < 1:
        return []
    items = []
    if kind == "nebula" and (after is None or after[0] == "field"):
        field = _field_where(object_class)
        if field is not None:
            where, params = field
            if after is not None:
                where += " AND (cell_i, cell_j, cell_k, n) > (?, ?, ?, ?)"
                params = [*params, *after[1]]
            rows = conn.execute(f"SELECT {_FIELD_COLUMNS} FROM nebula_field WHERE {where} "
                                "ORDER BY cell_i, cell_j, cell_k, n LIMIT ?", [*params, limit]).fetchall()
            items += [_cloud_item(row) for row in rows]
    room = limit - len(items)
    if room > 0:
        items += _scatter_items(conn, kind, room, after[1] if after is not None and after[0] == "scatter" else None, object_class)
    return items


def _scatter_items(conn, kind, limit, after_id, object_class):
    classes = _classes(conn, kind, object_class)
    if not classes:
        return []
    by_class = {(entry["kind"], entry["subtype"]): entry["descriptor"] for entry in classes}
    clauses, params = [], []
    for entry in classes:
        clauses.append("(kind = ? AND COALESCE(subtype, '') = ?)")
        params.extend([entry["kind"], entry["subtype"]])
    later = ""
    if after_id is not None:
        later, params = "id > ? AND ", [after_id, *params]
    rows = conn.execute(
        "SELECT id, uid, kind, COALESCE(subtype, '') AS subtype, shell_index, band_index, wedge_index, mass_solar"
        f" FROM phenomenon_scatter WHERE built_at IS NULL AND uid IS NOT NULL AND {later}({' OR '.join(clauses)})"
        " ORDER BY id LIMIT ?", [*params, limit]).fetchall()
    return [_scatter_item(row, kind, by_class[(row["kind"], row["subtype"])]) for row in rows]


def within(conn, center, distance, kinds):
    """
    The uncharted objects of the `kinds` within `distance` parsecs of `center`, each `(distance, kind, id, name, point, scatter_id)`:
    the scattered rows by their position and the field clouds by their center, both found by the cell key (`db.region`). `id` is the printed
    ID of the row (a field cloud's `field-i.j.k.n`); `scatter_id` is the row of `phenomenon_scatter` and `None` for a field cloud.
    """
    wanted = set(kinds) & set(TAIL_KINDS)
    found = []
    if not wanted:
        return found
    runs = None
    position = ("p.position_x_mpc / 1000.0", "p.position_y_mpc / 1000.0", "p.position_z_mpc / 1000.0")
    sql, runs = region.region_statement(
        "phenomenon_scatter", "p.id, p.uid, p.kind, p.shell_index, p.band_index, p.wedge_index, " + ", ".join(
            f"{axis} AS {name}" for axis, name in zip(position, ("x", "y", "z"))),
        position, center, distance, "p.built_at IS NULL AND p.uid IS NOT NULL", runs)
    for row in conn.execute(sql).fetchall():
        kind = scatter_type(row["kind"])
        if kind not in wanted:
            continue
        point = (row["x"], row["y"], row["z"])
        away = detmath.dist(center, point)
        if away <= distance:
            name = f"Uncharted {row['kind'].replace('-', ' ')} {row['shell_index']}.{row['band_index']}.{row['wedge_index']}"
            found.append((away, kind, _printed(row["uid"]), name, point, row["id"]))
    if "nebula" in wanted:
        position = ("p.center_x_pc", "p.center_y_pc", "p.center_z_pc")
        sql, runs = region.region_statement(
            "nebula_field", "p.cell_i, p.cell_j, p.cell_k, p.n, p.center_x_pc, p.center_y_pc, p.center_z_pc", position, center,
            distance, "p.is_stored = 0", runs)
        for row in conn.execute(sql).fetchall():
            point = (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])
            away = detmath.dist(center, point)
            if away <= distance:
                found.append((away, "nebula", cloud_id(row), f"Field cloud {cloud_id(row)[len(CLOUD_PREFIX):]}", point, None))
    return found


def search_rows(conn, types, limit, offset, field_codes=None, scatter_types=None):
    """
    A page of the uncharted objects the search page's Phenomena panel lists (MAP.400), by offset like the panel: the field clouds first (in
    the order of their key), then the scattered rows (by id).

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        types (iterable): The kinds asked for (those in `TAIL_KINDS` count).
        limit (int): Rows on the page.
        offset (int): Rows to skip.
        field_codes (iterable, optional): Only the field clouds of these classes (`M`, `N`); all by default.
        scatter_types (iterable, optional): Only the scattered rows of these kinds (a class tag leaves out the kinds it narrows: a scattered
            row has no class letter); all of `types` by default.

    Returns:
        tuple: `(rows, total)`: rows `type`, `id` (the scatter row's id, or a field cloud's printed ID), `name`, `phenomenon_class`,
        `sector_id` (`None`), `scatter_id`, `uncharted` (True); and how many there are in all.
    """
    types = [kind for kind in TAIL_KINDS if kind in set(types)]
    scatter_wanted = types if scatter_types is None else [kind for kind in types if kind in set(scatter_types)]
    field = _field_where(codes=field_codes) if "nebula" in types else None
    field_total = 0
    rows = []
    if field is not None:
        field_total = conn.execute(f"SELECT COUNT(*) AS n FROM nebula_field WHERE {field[0]}", field[1]).fetchone()["n"]
        if offset < field_total and limit > 0:
            for row in conn.execute(f"SELECT {_FIELD_COLUMNS} FROM nebula_field WHERE {field[0]} "
                                    "ORDER BY cell_i, cell_j, cell_k, n LIMIT ? OFFSET ?", [*field[1], limit, offset]).fetchall():
                rows.append({"type": "nebula", "id": cloud_id(row), "name": f"Field cloud {cloud_id(row)[len(CLOUD_PREFIX):]}",
                             "phenomenon_class": row["nebula_class"], "sector_id": None, "scatter_id": None, "uncharted": True})
    classes = [entry for entry in query._scattered_classes(conn) if entry["type"] in scatter_wanted]
    scatter_total = sum(entry["unbuilt"] for entry in classes)
    room = limit - len(rows)
    if classes and room > 0:
        clauses, params = [], []
        for entry in classes:
            clauses.append("(kind = ? AND COALESCE(subtype, '') = ?)")
            params.extend([entry["kind"], entry["subtype"]])
        by_class = {(entry["kind"], entry["subtype"]): entry["type"] for entry in classes}
        for row in conn.execute(
                "SELECT id, kind, COALESCE(subtype, '') AS subtype, shell_index, band_index, wedge_index FROM phenomenon_scatter"
                f" WHERE built_at IS NULL AND ({' OR '.join(clauses)}) ORDER BY id LIMIT ? OFFSET ?",
                [*params, room, max(0, offset - field_total)]).fetchall():
            rows.append({"type": by_class[(row["kind"], row["subtype"])], "id": row["id"],
                         "name": f"Uncharted {row['kind'].replace('-', ' ')} {row['shell_index']}.{row['band_index']}.{row['wedge_index']}",
                         "phenomenon_class": row["subtype"] or None, "sector_id": None, "scatter_id": row["id"], "uncharted": True})
    return rows, int(field_total + scatter_total)
