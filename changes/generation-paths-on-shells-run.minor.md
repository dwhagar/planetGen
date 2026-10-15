### Changed
- `planetgen galaxy` and `plan` run on shells to the end: single sectors, bands, shells, spans, courses and random starts fill by `(shell, band, wedge)`; the plan stages, the backfill and the sector paths read and write the shell columns; the galactic nucleus is the one core sector. The bright-star backfill and the sector-stats lock read `sector_stats` by its primary key, so concurrent backfills lock rows in address order (they could deadlock on the level index).
- The block-totals step of a `galaxy` run and `plan --fill-block-stats` are gone with the `galaxy_block_stats` table; the Galaxy Map's ladder brings its totals back with MAP.205.
- `SectorCell.contains` refuses a non-finite point (a cap or the core accepted NaN).

TODO items resolved: GEN.227 / Added: none
