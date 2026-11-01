# planetgen/db/chart_ledger.py

"""
The chart ledger (v197, GEN.1127, docs/design/systems/region-object-store.md sections 5, 6 and 10)
=================================================================================================

Which sector cells have been charted. A cell is a whole sector, the galaxy's unit of randomness, so it is charted or it is not: there is no
partial state and no stored boundary, and the uncharted remainder of a region is every cell with no row here. A row records the run and the
generator epoch that charted the cell, when, and when the rogue top-up of its region run filled it (`topped_up_at`, NULL until it does; a cell
keeps one top-up, GEN.1120). The row is keyed by the cell's address, not by `sectors.id`, so it survives a sector row that is replaced.

Written by the save of a sector (`store.save_sector` calls `record`) and cleared when a sector is deleted (`store.forget_sector_fill` calls
`forget`); read by `store.filled_sector_addresses` (which cells a run leaves alone) and by `galaxy/epoch_record.py` through `state`.

    record(conn, (shell, band, wedge), run_id, epoch)     # the save of a sector
    forget(conn, address)                                 # a sector was deleted
    state(conn, address, current_epoch)                   # "charted", "stale" (charted under an older epoch) or "uncharted"
    charted_addresses(conn)                               # {(shell, band, wedge)}
    mark_topped_up(conn, address)                         # the top-up filled the cell; True when this call did it
"""

UNCHARTED = "uncharted"
CHARTED = "charted"
STALE = "stale"
STATES = (UNCHARTED, CHARTED, STALE)
"""tuple: The charted states of a cell: no row; a row made under the current epoch; a row made under an older one."""


def record(conn, address, run_id, epoch):
    """
    Charts the cell `address` (`(shell, band, wedge)`): a new row, or the old one made again (a regenerated cell is charted afresh, so its
    time moves and its top-up is forgotten: the top-up measured a shortfall that is gone).
    """
    conn.execute(
        "INSERT INTO chart_ledger (shell_index, band_index, wedge_index, run_id, generator_epoch, charted_at, topped_up_at)"
        " VALUES (?, ?, ?, ?, ?, CURRENT_TIMESTAMP(3), NULL)"
        " ON DUPLICATE KEY UPDATE run_id = VALUES(run_id), generator_epoch = VALUES(generator_epoch),"
        " charted_at = CURRENT_TIMESTAMP(3), topped_up_at = NULL",
        (*address, run_id, epoch))


def forget(conn, address):
    """Uncharts the cell `address` (its sector was deleted). Nothing for a cell that was never charted."""
    conn.execute("DELETE FROM chart_ledger WHERE shell_index = ? AND band_index = ? AND wedge_index = ?", tuple(address))


def row(conn, address):
    """The ledger row of a cell as a dict (`run_id`, `generator_epoch`, `charted_at`, `topped_up_at`), or `None` when it is uncharted."""
    return conn.execute(
        "SELECT run_id, generator_epoch, charted_at, topped_up_at FROM chart_ledger"
        " WHERE shell_index = ? AND band_index = ? AND wedge_index = ?", tuple(address)).fetchone()


def state(conn, address, current_epoch):
    """
    `UNCHARTED`, `CHARTED` or `STALE` for the cell `address`: a cell is stale when an epoch since it was charted changed how the fill
    makes a sector (`galaxy/epoch_record.needs`); a row with no recorded epoch counts as the first epoch.
    """
    found = row(conn, address)
    if found is None:
        return UNCHARTED
    from planetgen.galaxy import epoch_record
    return STALE if epoch_record.needs("fill", found["generator_epoch"], current_epoch) else CHARTED


def charted_addresses(conn):
    """The `(shell, band, wedge)` of every charted cell."""
    return {(r["shell_index"], r["band_index"], r["wedge_index"])
            for r in conn.execute("SELECT shell_index, band_index, wedge_index FROM chart_ledger").fetchall()}


def count(conn):
    """How many cells are charted."""
    return conn.execute("SELECT COUNT(*) AS n FROM chart_ledger").fetchone()["n"]


def mark_topped_up(conn, address):
    """
    Notes that the rogue top-up filled the cell `address` (GEN.1120): a cell keeps one top-up, so a cell already topped up, or not
    charted, is left alone.

    Returns:
        bool: Whether this call set the time.
    """
    return conn.execute(
        "UPDATE chart_ledger SET topped_up_at = CURRENT_TIMESTAMP(3) WHERE shell_index = ? AND band_index = ? AND wedge_index = ?"
        " AND topped_up_at IS NULL", tuple(address)).rowcount > 0
