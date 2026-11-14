# planetgen/db/chart_ledger.py

"""
The chart ledger (v201, GEN.1127, docs/design/systems/region-object-store.md sections 5, 6 and 10)
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
        " ON DUPLICATE KEY UPDATE run_id = ?, generator_epoch = ?,"
        " charted_at = CURRENT_TIMESTAMP(3), topped_up_at = NULL",
        (*address, run_id, epoch, run_id, epoch))


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


def charted_among(conn, addresses):
    """The cells of `addresses` that are charted (VIEW.25 reads a few hundred at a time, by the primary key)."""
    found = set()
    keys = sorted({(int(a[0]), int(a[1]), int(a[2])) for a in addresses})
    for start in range(0, len(keys), 500):
        chunk = keys[start:start + 500]
        marks = ", ".join("(?, ?, ?)" for _ in chunk)
        for row in conn.execute(f"SELECT shell_index, band_index, wedge_index FROM chart_ledger WHERE (shell_index, band_index, wedge_index) IN ({marks})",
                                tuple(value for key in chunk for value in key)).fetchall():
            found.add((row["shell_index"], row["band_index"], row["wedge_index"]))
    return found


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


_NO_ROW = ("NOT EXISTS (SELECT 1 FROM chart_ledger c WHERE c.shell_index = s.shell_index AND c.band_index = s.band_index AND c.wedge_index = s.wedge_index)")
FILLED_SELECTS = (
    "SELECT DISTINCT s.shell_index, s.band_index, s.wedge_index, s.run_id, s.generator_epoch FROM star_systems ss JOIN sectors s ON s.id = ss.sector_id "
    f"WHERE s.shell_index IS NOT NULL AND {_NO_ROW}",
    "SELECT DISTINCT s.shell_index, s.band_index, s.wedge_index, s.run_id, s.generator_epoch FROM sector_stats st JOIN sectors s ON s.shell_index = st.shell_index "
    f"AND s.band_index = st.band_index AND s.wedge_index = st.wedge_index WHERE st.bright_level_sol = 0 AND {_NO_ROW}",
)
"""tuple: The placed sectors that were charted though the ledger holds no row for their cell (`repair`): one that holds a star system, and one whose fill was recorded (`sector_stats.bright_level_sol` 0). They start
from the small tables, not from every sector row, since the placeholders of the scatters are most of those. A rogue planet alone does not make a sector charted: one drifting into a cell nothing was generated in makes an
uncharted placeholder row (`store.uncharted_sector_id`). OPS.180."""


def diagnosis(conn):
    """
    The counts that say whether the ledger agrees with the sectors (OPS.180): `placed` (sector rows with a grid address), `ledger` (cells charted), `charted` (placed sectors whose cell the ledger holds) and `missing`
    (placed sectors that hold a system or a recorded fill but whose cell it does not).
    """
    placed = conn.execute("SELECT COUNT(*) AS n FROM sectors WHERE shell_index IS NOT NULL").fetchone()["n"]
    charted = conn.execute("SELECT COUNT(*) AS n FROM sectors s WHERE s.shell_index IS NOT NULL AND EXISTS (SELECT 1 FROM chart_ledger c WHERE c.shell_index = s.shell_index"
                           " AND c.band_index = s.band_index AND c.wedge_index = s.wedge_index)").fetchone()["n"]
    missing = len({(r["shell_index"], r["band_index"], r["wedge_index"]) for select in FILLED_SELECTS for r in conn.execute(select).fetchall()})
    return {"placed": int(placed), "ledger": count(conn), "charted": int(charted), "missing": missing}


def repair(conn):
    """
    Writes the ledger row of every placed sector that holds a system or a recorded fill and has none (OPS.180), as its save would have: the run and epoch of the sector's row, the time now. A sector
    in a cell the ledger does not hold reads as an uncharted placeholder everywhere (the Sectors list, the Galaxy page and the Generate page's counts), so a charted sector left out of the ledger vanished
    from them. A second run finds nothing. Returns the counts of `diagnosis` taken before, with `repaired`, the rows written.
    """
    found = diagnosis(conn)
    repaired = 0
    for select in FILLED_SELECTS:
        repaired += conn.execute(
            "INSERT IGNORE INTO chart_ledger (shell_index, band_index, wedge_index, run_id, generator_epoch, charted_at, topped_up_at) "
            "SELECT f.shell_index, f.band_index, f.wedge_index, f.run_id, f.generator_epoch, CURRENT_TIMESTAMP(3), NULL FROM (" + select + ") f").rowcount
    conn.commit()
    return {**found, "repaired": repaired}
