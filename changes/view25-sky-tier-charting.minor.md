**Headline:** Charting can now build the bright and massive stars a system's sky needs out to the reach of each tier, and remembers which spheres it has drawn.

### Added

- `planetgen galaxy`, `sectors` and `backfill` take `--sky-limit MAG` (VIEW.25): after the bright-star backfill, each 50 parsec observer cell the new sectors lie in gets the tiers of its sky table (`generation/sky_tiers.py`) drawn over spheres of their reach, the luminosity bands from the brightest down and then the mass rings. Each draw is a unit of the backfill command, costed from the table's block visits, committed in batches and resumed after a stop, and one sky draw runs at a time in the galaxy (`GET_LOCK`). The stars are the ones the cell backfill would give for the same cells (`generation/sky_charting.py`).
- The table `sky_tier_draws` (schema v205, `db/sky_tier_ledger.py`) records each sphere drawn. A cell inside a complete sphere is complete to the band's floor without a `sector_stats` row; the bright-star bands, the mass passes and the fill settle the ledger's depth into the rows they read before they draw, so a later backfill of a cell does not draw the stars a sphere already did, whichever came first.

### Changed

- `db/sky_sphere.cell_levels` and so the sky's `filling` and `built_to_pc` read the ledger, so a sphere too large to list is complete down to the deepest sphere that holds all of it.

generation-output: unchanged
