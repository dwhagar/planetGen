**Headline:** Every sector cell that the star and object scatters put something in now has an entry in the Sectors list as an uncharted sector, so the list, the filters, the search and the sector pages can reach what long-range scans detected there.

### Changed
- New `db/uncharted_index.py`: `index` makes the placeholder `sectors` row (MAP.399) for every cell that holds a waiting bright star or an unbuilt scattered phenomenon and has none, in batches of 500 cells per transaction through `contention.run_with_retries` (gate `uncharted-sectors`), a shell at a time; a second run adds nothing. It runs when a scatter ends (beside the block light totals, in `refresh_block_lights` and `refresh_block_lights_unit`), in `planetgen fill-block-totals`, and in migration 0226 for a galaxy that was scattered before. A cell whose objects are later built or cleared keeps its row as an empty uncharted sector (MAP.403 leaves those out of the list by default).
- A new plan clears the placeholder rows nothing lives in (`uncharted_index.clear_unoccupied`) when it clears the scatters of the old outline, and `store.has_galaxy_sectors` (the guard on a different seed) no longer counts a placeholder as a charted sector.
- The Galaxy Map's placed-sector queries (`galaxy_placed_sectors`, `galaxy_sectors_in_view`, `galaxy_sectors_in_box`) leave the placeholders out, so they do not draw as generated sectors.
- The totals that count generated sectors leave the placeholders out (`query.charted_sectors_sql`): the block totals (`block_stats.rebuild`), the regions (`region_stats`), the shell spec counts, the stage sums and the sector count of a box. Migration 0168, which runs the block rebuild on a schema older than the chart ledger, reads every row (`charted=False`). Test fixtures that insert sector rows by hand record the chart ledger row a saved sector has (`tests/chart_support.py`).
- Migration 0221 changes no schema; `helpers.StoreConnection` results gain `rowcount` and `lastrowid`.
- Measured: 30,000 cells in 2.7 seconds on a lane machine (about 11,000 cells a second, so a million cells is about 90 seconds); a second run over them takes 0.1 seconds.
- New tests: `test_uncharted_index.py` (the rows, a second run, listing, empty cells, the migration, the seed guard, the new-plan clear).

generation-output: unchanged
