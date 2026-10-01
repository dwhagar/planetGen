### Changed

- Bright stars now come after the sectors (GEN.30). A new galaxy generates its first sectors and then scatters the bright stars galaxy-wide, leaving those sectors out (`generate.py galaxy --then-scatter`).
- The bright-star backfill runs once, after a run has generated every sector it was asked for. By default it backfills only around the requested sector: the random start, the center sector or the slot address, or for ring, column, shell and block runs the generated sector nearest the middle. `--backfill-from all` (on the Generate page, "Backfill from every generated sector (farthest out)") backfills around every generated sector instead, and `--backfill-from none` skips the backfill.
- The scatter always leaves filled sectors out. The "Leave filled sectors out" checkbox is gone, and `--force` is accepted but no longer needed. A generated sector never gets scattered or backfilled stars.
