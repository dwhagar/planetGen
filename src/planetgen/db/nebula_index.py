# planetgen/db/nebula_index.py

"""
The nebula index (GEN.202): one read that gives every nebula of a box, whatever table holds it, with no object twice.

A nebula lives in up to three places before and after its sector is generated: the stored nebulae and supernova remnants
(`query.galaxy_clouds_in_box`), the field clouds the galaxy's seed places (`field_clouds.field_clouds_in_box`, listed only until a sector stores
one) and the scattered planetary nebulae and remnants of `phenomenon_scatter` (`query.galaxy_scattered_nebulae_in_box`, listed only until a sector
builds one). Each reader leaves out what is stored already by its own flag, so a node's lists do not repeat. `nebulae_in_box` reads the three and
also makes the object ID the key: a cloud, a field row or a scatter row that holds the Object ID of a stored one is dropped whatever its flag says
(a flag a redrawn plan has not set again, a row built before its mark was written), so the one answer lists each object once. The stored one is
kept, then the field cloud, then the scattered mark.

An entry is the reader's dict (`type`, `id`, `name`, `descriptor`, `x`/`y`/`z`, and for a cloud `class` and `radius_pc`) with `source`
(`"stored"`, `"field"` or `"scattered"`) and `object_id` (the Object ID as hex, `None` for an object that has none yet).
"""

from planetgen.db import field_clouds, query

STORED, FIELD, SCATTERED = "stored", "field", "scattered"
"""str: Where an entry of `nebulae_in_box` was read from."""

_STORED_TABLES = {"nebula": "nebulae", "supernova_remnant": "supernova_remnants"}
"""dict: The table of each stored cloud's `type`."""


def _hex(value):
    return value.hex() if isinstance(value, (bytes, bytearray)) else None


def _stored_object_ids(conn, clouds):
    """`{(type, id): object ID hex}` of the stored clouds that have one."""
    found = {}
    for type_label, table in _STORED_TABLES.items():
        ids = [cloud["id"] for cloud in clouds if cloud["type"] == type_label]
        if not ids:
            continue
        marks = ", ".join("?" for _ in ids)
        for row in conn.execute(f"SELECT id, uid FROM {table} WHERE id IN ({marks})", ids).fetchall():
            if row["uid"] is not None:
                found[(type_label, row["id"])] = _hex(row["uid"])
    return found


def _field_object_ids(conn, clouds):
    """`{cloud id: object ID hex}` of the field clouds that carry one (a stored cloud's row does)."""
    found = {}
    for cloud in clouds:
        cell_i, cell_j, cell_k, number = (int(part) for part in cloud["id"].split("-", 1)[1].split("."))
        row = conn.execute("SELECT uid FROM nebula_field WHERE cell_i = ? AND cell_j = ? AND cell_k = ? AND n = ?",
                           (cell_i, cell_j, cell_k, number)).fetchone()
        if row is not None and row["uid"] is not None:
            found[cloud["id"]] = _hex(row["uid"])
    return found


def _scattered_object_ids(conn, points):
    """`{mark id: object ID hex}` of the scattered marks whose row has one."""
    found = {}
    rows = {int(point["id"][1:]): point["id"] for point in points}
    if rows:
        marks = ", ".join("?" for _ in rows)
        for row in conn.execute(f"SELECT id, uid FROM phenomenon_scatter WHERE id IN ({marks})", list(rows)).fetchall():
            if row["uid"] is not None:
                found[rows[row["id"]]] = _hex(row["uid"])
    return found


def nebulae_in_box(conn, lo, hi, edge_pc, limit):
    """
    Every nebula and supernova remnant of the box `[lo, hi)`: the stored ones whose sphere reaches into it, the field clouds no sector has stored
    and the scattered planetary nebulae and remnants no sector has built, each object once by its Object ID.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        lo (tuple): `(x, y, z)` lower corner, parsecs.
        hi (tuple): `(x, y, z)` upper corner, parsecs.
        edge_pc (float): The galaxy's sector edge, parsecs (the scatter is filed by sector).
        limit (int): The most to read from each source (the answer holds at most three times this).

    Returns:
        list[dict]: The stored clouds largest first, then the field clouds largest first, then the scattered marks; each with `source` and `object_id`.
    """
    stored = query.galaxy_clouds_in_box(conn, lo, hi, max_clouds=limit)
    field = field_clouds.field_clouds_in_box(conn, lo, hi, limit)
    scattered = query.galaxy_scattered_nebulae_in_box(conn, lo, hi, edge_pc, limit)
    stored_ids = _stored_object_ids(conn, stored)
    field_ids = _field_object_ids(conn, field)
    scattered_ids = _scattered_object_ids(conn, scattered)
    seen, found = set(), []
    for source, entries, ids, key in ((STORED, stored, stored_ids, lambda entry: (entry["type"], entry["id"])),
                                      (FIELD, field, field_ids, lambda entry: entry["id"]),
                                      (SCATTERED, scattered, scattered_ids, lambda entry: entry["id"])):
        for entry in entries:
            object_id = ids.get(key(entry))
            if object_id is not None and object_id in seen:
                continue
            if object_id is not None:
                seen.add(object_id)
            found.append({**entry, "source": source, "object_id": object_id})
    return found
