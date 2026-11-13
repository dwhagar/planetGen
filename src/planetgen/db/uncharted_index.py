# planetgen/db/uncharted_index.py

"""
The placeholder sector rows of uncharted cells (MAP.402, ended by OPS.154).

MAP.402 made a `sectors` row for every cell a scatter left a waiting bright star or an unbuilt scattered phenomenon in, so the Sectors list, the search and the
sector pages could reach what the long-range scans detected. At full size that is up to 1.5 billion rows of 408 bytes, about 600 GB, and everything a row holds
is a function of the cell's address. OPS.154 ended it: the cells are read from the block light parts (`db.block_lights`, `db.query.virtual_cells`) and the stage
that made the rows is gone. What stays is the delete: `clear_unoccupied` removes the placeholder rows nothing lives in, for a new plan (whose scatters are cleared
with the old outline) and for the migration (0247) that sheds the rows an older run made. A cell a rogue planet drifted into keeps its row
(`store.uncharted_sector_id`): that row is the planet's address.
"""

_UNCHARTED = (
    "(sec.shell_index IS NOT NULL AND NOT EXISTS (SELECT 1 FROM chart_ledger cl WHERE cl.shell_index = sec.shell_index "
    "AND cl.band_index = sec.band_index AND cl.wedge_index = sec.wedge_index))")
_UNOCCUPIED = ("NOT EXISTS (SELECT 1 FROM rogue_planets rp WHERE rp.sector_id = sec.id) "
               "AND NOT EXISTS (SELECT 1 FROM star_systems ss WHERE ss.sector_id = sec.id)")
CLEAR_CHUNK = 5000


def clear_unoccupied(conn, on_phase=None, on_progress=None):
    """
    Deletes the placeholder rows nothing lives in (no rogue planet filed in them, no system), for a new plan: the scatters they were made for are
    cleared with the old outline, so they would be empty rows for cells that may not exist in the new one. A row something references stays (the
    delete is skipped for it by the database's own keys). `on_progress(deleted)` is called after each batch. Returns the rows deleted.
    """
    deleted, last = 0, 0
    while True:
        ids = [r["id"] for r in conn.execute(
            f"SELECT sec.id FROM sectors sec WHERE sec.id > ? AND {_UNCHARTED} AND {_UNOCCUPIED} ORDER BY sec.id LIMIT ?", (last, CLEAR_CHUNK)).fetchall()]
        if not ids:
            break
        last = ids[-1]
        deleted += conn.execute(f"DELETE IGNORE FROM sectors WHERE id IN ({', '.join('?' for _ in ids)})", ids).rowcount
        conn.commit()
        if on_progress is not None:
            on_progress(deleted)
    if on_phase is not None and deleted:
        on_phase(f"Cleared {deleted:,} uncharted sector rows nothing lives in.")
    return deleted
