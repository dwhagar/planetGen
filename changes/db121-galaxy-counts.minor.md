**Headline:** The galaxy's counts are kept as they change, so a page no longer counts millions of systems again every time a sector is added.

### Changed

- DB.121: two new tables, `sector_facets` (for each generated sector, its systems, stars, planets, moons, belts and habitable worlds, and its count of each class, body type, life class, spectral and Yerkes class, binary state and octant) and `galaxy_totals` (one row, their sum plus the galaxy's mass and luminosity). They are written in the same transaction as the sector's save, edit and deletion (`db/galaxy_counts.py`, called from `record_sector_stats`, `refresh_edited_sectors` and `forget_sector_fill`), by the difference between the sector's old rows and its new ones. The migration makes them from the stored bodies at once, and `planetgen fill-block-totals` makes them again.
- The Systems list's unfiltered count and filter menu, and the search page's type, spectral class, Yerkes class, planet and moon class, body type and life class counts, read these rows (plus a small count of what no sector holds) while they are whole, and count the tables otherwise. The stored count (`countcache`) is still remade when a sector is added, but the remaking is now one read.
- `check-counts` compares the stored counts with a full count of the tables when it runs over the whole galaxy (a warning), and says to run `planetgen fill-block-totals` when they differ.
- Nothing stored about a system, planet or star changes.

generation-output: unchanged
