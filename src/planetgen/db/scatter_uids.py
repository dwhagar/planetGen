"""
planetgen/db/scatter_uids.py

The Object ID of every scatter row (GEN.1254; docs/design/architecture/object-ids.md, "Scatter rows"). A star of `bright_stars` and a
scattered phenomenon of `phenomenon_scatter` are born in the sector they sit in, so they have an ID before the sector is charted:
the sector, a serial and body number 0. The scatter writers (`store.insert_bright_stars`, `store.insert_phenomenon_scatter`) call `assign`
when they finish, and the rows they just wrote that have no `uid` get one.

The serial is the row's rank in its sector among the rows written in the same call, brightest first (`luminosity_w` for a star,
`map_luminosity_sol` for a phenomenon), with the star's own `seed` and then its row number for a tie (the seed is a draw of 63 bits, so the
row number is a last resort and the order does not depend on how the rows were inserted). The rank continues from the sector's counter in
`id_counters`, so a later call for the same sector (a backfill that adds dimmer stars) continues the numbers and never repeats one. The same
seed written in the same order gives the same IDs.

The two tables share a sector's serials without sharing a counter: stars take serials from 0, phenomena from `PHENOMENON_SERIAL_BASE` (half the
serial space, 33 million each). `_UidIssuer` (`store.py`) starts the serials of the sector's own fill past the stars' counter, and a charted
system round a scatter star, or a built scatter phenomenon, keeps the ID its row has.

The work is set-based in the database (a window function and three joined updates), so a batch of a million rows costs a few statements.
"""
from planetgen.galaxy import object_uid

STAR_COUNTER = "scatstar"
PHENOMENON_COUNTER = "scatphen"
"""str: The `id_counters.kind`s (the column is 8 characters): serials given out so far per birth sector, stars and phenomena."""

PHENOMENON_SERIAL_BASE = 1 << (object_uid.LAYOUT.serial_count_bits - 1)
"""int: The first serial of a scattered phenomenon: the stars have the lower half of the serial space and the phenomena the upper."""

_TABLES = {
    "bright_stars": (STAR_COUNTER, 0, "luminosity_w DESC, seed, id"),
    "phenomenon_scatter": (PHENOMENON_COUNTER, PHENOMENON_SERIAL_BASE, "map_luminosity_sol DESC, kind, subtype, seed, id"),
}
"""dict: Per scatter table, its counter, its first serial and the order of the rank."""

COUNTERS = (STAR_COUNTER, PHENOMENON_COUNTER)

_LAYOUT = object_uid.LAYOUT
_SECTOR_SQL = (f"(((CAST(shell_index AS UNSIGNED) << {_LAYOUT.band_bits}) | CAST(band_index + {_LAYOUT.band_bias} AS UNSIGNED))"
               f" << {_LAYOUT.wedge_bits}) | CAST(wedge_index AS UNSIGNED)")
"""str: The birth sector of a scatter row as the integer `object_uid.pack` makes of its address."""
_SECTOR_HEX_DIGITS = _LAYOUT.sector_bits // 4
_BELOW_SECTOR_HEX_DIGITS = (_LAYOUT.serial_bits + _LAYOUT.body_bits) // 4


def _create_work_tables(conn):
    conn.execute(
        "CREATE TEMPORARY TABLE IF NOT EXISTS scatter_uid_rank (id BIGINT UNSIGNED NOT NULL PRIMARY KEY,"
        " rn BIGINT UNSIGNED NOT NULL, scope BINARY(7) NOT NULL, KEY (scope))")
    conn.execute(
        "CREATE TEMPORARY TABLE IF NOT EXISTS scatter_uid_sector (scope BINARY(7) NOT NULL PRIMARY KEY,"
        " cnt BIGINT UNSIGNED NOT NULL, base BIGINT UNSIGNED NOT NULL DEFAULT 0)")
    conn.execute("DELETE FROM scatter_uid_rank")
    conn.execute("DELETE FROM scatter_uid_sector")


def assign(conn, table):
    """
    Gives every row of `table` that has no `uid` its Object ID (module note) and moves each sector's counter on. The caller's transaction.

    Args:
        conn (Connection): Open.
        table (str): `bright_stars` or `phenomenon_scatter`.

    Returns:
        int: Rows given an ID.

    Raises:
        ValueError: A table that is not a scatter table, or a sector whose serials the table's half of the space cannot hold.
    """
    if table not in _TABLES:
        raise ValueError(f"not a scatter table: {table!r}")
    counter, first, order = _TABLES[table]
    _create_work_tables(conn)
    conn.execute(
        "INSERT INTO scatter_uid_rank (id, rn, scope)"
        f" SELECT id, ROW_NUMBER() OVER (PARTITION BY shell_index, band_index, wedge_index ORDER BY {order}) - 1,"
        f" UNHEX(LPAD(HEX({_SECTOR_SQL}), {_SECTOR_HEX_DIGITS}, '0'))"
        f" FROM {table} WHERE uid IS NULL")
    given = conn.execute("SELECT COUNT(*) AS n FROM scatter_uid_rank").fetchone()["n"]
    if not given:
        return 0
    conn.execute("INSERT INTO scatter_uid_sector (scope, cnt) SELECT scope, COUNT(*) FROM scatter_uid_rank GROUP BY scope")
    conn.execute(
        "INSERT INTO id_counters (kind, scope, next_value) SELECT ?, scope, 0 FROM scatter_uid_sector"
        " ON DUPLICATE KEY UPDATE next_value = next_value", (counter,))
    conn.execute("UPDATE scatter_uid_sector s JOIN id_counters c ON c.kind = ? AND c.scope = s.scope SET s.base = c.next_value", (counter,))
    room = (1 << (_LAYOUT.serial_count_bits - (1 if table == "bright_stars" else 0))) - first
    over = conn.execute("SELECT scope FROM scatter_uid_sector WHERE base + cnt > ? LIMIT 1", (room,)).fetchone()
    if over is not None:
        raise ValueError(f"a sector's {table} rows have used all of the serials set aside for them ({room})")
    conn.execute("UPDATE id_counters c JOIN scatter_uid_sector s ON c.kind = ? AND c.scope = s.scope SET c.next_value = s.base + s.cnt", (counter,))
    conn.execute(
        f"UPDATE {table} t JOIN scatter_uid_rank r ON r.id = t.id JOIN scatter_uid_sector s ON s.scope = r.scope"
        f" SET t.uid = UNHEX(CONCAT(HEX(r.scope), LPAD(HEX((? + s.base + r.rn) << {_LAYOUT.body_bits}), {_BELOW_SECTOR_HEX_DIGITS}, '0')))",
        (first,))
    return given


def next_star_serial(conn, shell, band, wedge):
    """The first serial of the sector at `(shell, band, wedge)` that no scatter star has: where the sector's own fill starts (0 for a sector
    with no scatter stars)."""
    row = conn.execute("SELECT next_value FROM id_counters WHERE kind = ? AND scope = ?",
                       (STAR_COUNTER, object_uid.sector_bytes(shell, band, wedge))).fetchone()
    return 0 if row is None else int(row["next_value"])


def forget(conn, counter=None):
    """Deletes the counters of one scatter table (`STAR_COUNTER`, `PHENOMENON_COUNTER`) or both, when its rows are cleared: a scatter made again
    starts again from serial 0."""
    kinds = COUNTERS if counter is None else (counter,)
    for kind in kinds:
        conn.execute("DELETE FROM id_counters WHERE kind = ?", (kind,))
