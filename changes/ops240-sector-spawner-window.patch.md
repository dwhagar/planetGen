**Headline:** a galaxy run makes its sectors one at a time at any size, and the engine writes them to its ledger a window at a time.

### Changed
- OPS.240: the sector stages (`generate-start-sector`, `generate-neighborhood`, `generate-sectors`) have one unit for each sector however many sectors the galaxy has (the chunking of more than 4,000 sectors a stage into units of several, `SECTOR_UNITS_MAX`, is gone). The cut writes the first window of unit rows (4 for each of the stage's workers, between 200 and 2,000) and a `spawn` row that holds the cursor; the launcher's turn writes the next window when the unit rows written and not finished fall under half of one (`planrows.spawn_window`). The stage still counts every sector (bar, ETA, the follower's total) and a run that carries on from another, or a held run that is retried, picks the cursor up from its row.
- The sectors of the windows after the first wait in a packed, checksummed file in the state folder (`queue/itemfile.py`, `stage-items/<run>-<stage>.items`, 12 bytes a sector), which the stage's end deletes; a file that is missing or torn fails the spawn row and holds the run with the reason, never a stage with fewer sectors.

generation-output: unchanged (no golden, no epoch, no lock item)
