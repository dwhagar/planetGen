# planetgen/db/cloudcontents.py

"""
What a nebula or supernova remnant contains (UX.144): every row whose `inside_nebula_id` / `inside_remnant_id` names it,
across the tables that carry those columns. The reverse of `store.surrounding_cloud`, which finds the cloud a row is in.
"""

CONTAINER_COLUMNS = {"nebula": "inside_nebula_id", "supernova_remnant": "inside_remnant_id"}
"""dict: The kinds that can contain something, and the column that names them."""

CONTENT_TABLES = (
    ("system", "star_systems", ""),
    ("nebula", "nebulae", ""),
    ("asteroid_field", "asteroid_fields", ""),
    ("rogue_planet", "rogue_planets", ""),
    ("interstellar_comet", "interstellar_comets", ""),
    # A system's own compact remnant is part of that system, not a free-standing object.
    ("black_hole", "black_holes", " AND star_id IS NULL"),
    ("neutron_star", "neutron_stars", " AND star_id IS NULL"),
)
"""tuple: `(kind, table, extra condition)` for each table that can sit inside a cloud, in the order they are listed."""

DEFAULT_LIMIT = 100
"""int: The most contained objects `cloud_contents` lists (the counts are always complete)."""


def cloud_contents(conn, cloud_kind, cloud_id, limit=DEFAULT_LIMIT):
    """
    What sits inside one nebula or supernova remnant.

    Args:
        conn (Connection): An open connection.
        cloud_kind (str): `"nebula"` or `"supernova_remnant"`.
        cloud_id (int): The cloud's row id.
        limit (int): The most objects to list (nearest to a stable order: by kind, then name, then id).

    Returns:
        dict: `counts` (`{kind: how many}` for each kind with any) and `items` (`[{"kind", "id", "name"}]`, at most
            `limit`; `id` is the row id), and `total`.

    Raises:
        ValueError: If `cloud_kind` can't contain anything.
    """
    column = CONTAINER_COLUMNS.get(cloud_kind)
    if column is None:
        raise ValueError(f"a {cloud_kind} contains nothing: use one of {', '.join(CONTAINER_COLUMNS)}")
    counts = {}
    items = []
    for kind, table, extra in CONTENT_TABLES:
        # A nebula never lists itself.
        skip = " AND id <> ?" if table == "nebulae" and cloud_kind == "nebula" else ""
        params = (cloud_id, cloud_id) if skip else (cloud_id,)
        found = conn.execute(f"SELECT COUNT(*) AS n FROM {table} WHERE {column} = ?{extra}{skip}", params).fetchone()["n"]
        if not found:
            continue
        counts[kind] = found
        room = limit - len(items)
        if room > 0:
            rows = conn.execute(
                f"SELECT id, name FROM {table} WHERE {column} = ?{extra}{skip} ORDER BY name, id LIMIT {int(room)}",
                params).fetchall()
            items.extend({"kind": kind, "id": row["id"], "name": row["name"]} for row in rows)
    return {"counts": counts, "items": items, "total": sum(counts.values())}
