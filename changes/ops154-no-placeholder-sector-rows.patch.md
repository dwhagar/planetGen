**Headline:** The galaxy no longer stores a sector row for every uncharted cell that holds a waiting star or phenomenon: the block light stage makes none, a migration deletes the ones already made (up to 600 GB at full size), and the Sectors list, the counts, the sector page and the search find those cells from the block light parts (part 2 of OPS.154).

### Changed
- **The block light stage is the roll-up alone (OPS.154).** `run_plan.make_block_lights` makes the totals and no sector rows; `refresh-block-lights` is one unit (the runs of shells of the uncharted index, OPS.153, are gone), and `fill-block-totals` repairs the totals and says so in its closing line without a count of rows. `uncharted_index.index`, `shell_runs`, `shell_part_counts` and the other helpers of the index are removed; `clear_unoccupied` stays (a new plan uses it) and reports each batch.
- **Migration 0254 deletes the rows already made.** It removes the placeholder rows nothing lives in (no rogue planet, no system; not in the chart ledger) 5,000 at a time, each batch committed, so a stopped run is started again and the database never holds one transaction of that size. A row a rogue planet or a system refers to stays: it is that body's address. No schema change, no epoch (the rows are derived and are not in the golden). Migration 0227 does nothing now, as there is nothing for it to make.
- **Search finds a cell with no row by its designation.** With uncharted sectors asked for, typing a cell's whole designation in the Sectors box returns the cell (`virtual: true`), as it would a row.
- **While part 1's guard lasts:** nothing changes for a galaxy that still has its bulk rows until the migration has run; after it the list reads the cells from the parts.

generation-output: unchanged
generation-affects: none
