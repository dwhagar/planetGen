# planetgen/describe/host_star.py
"""The stored star a body goes round, for the derived rows of the sheet: the star of its `star_id`, else the first star of its
system (a body stored before stars had ids of their own, or one whose star row is gone)."""


def row(conn, star_id, star_system_id, columns):
    """The `stars` row (`columns`, a SQL list) of the star `star_id`, else the first star of `star_system_id`, else `None`."""
    if star_id:
        star = conn.execute(f"SELECT {columns} FROM stars WHERE id = ?", (star_id,)).fetchone()
        if star is not None:
            return star
    if star_system_id:
        return conn.execute(f"SELECT {columns} FROM stars WHERE star_system_id = ? ORDER BY id LIMIT 1", (star_system_id,)).fetchone()
    return None


def of(conn, body, columns):
    """`row` for a stored planet row `body` (its `star_id` and `star_system_id` columns)."""
    return row(conn, body["star_id"], body["star_system_id"], columns)
