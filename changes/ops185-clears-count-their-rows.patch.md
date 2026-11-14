**Headline:** The long clears of a plan and a scatter (the bright stars, the phenomena, the nebula field and the block light totals) now go a batch of rows at a time and say how many they have done, so a slow machine never shows them as quiet.

### Changed

- OPS.185: the units of OPS.170 that can pass the 300 s quiet floor on a slow host now work in steps that move the unit's heartbeat and say a real count, as OPS.184 did for the generator battery's digest. New `db/batched.py` (`repeat`: a `DELETE` or an `UPDATE` run a batch of rows at a time, each batch committed, `unitbeat.moving()` and a `Ticker` line with the rows so far) is used for the sectors' bright star levels, the spheres the star draw kept (`sky_tier_draws`), the bright stars of one pass, the Galaxy Map's block light totals and parts, and the members and nodes of the scattered stars; `scatter_systems.delete_member_trees_scanning` finds the rows of a pass a chunk of ids at a time and removes their member trees, where one `DELETE` over a join was minutes; `uncharted_index.clear_unoccupied` and the nebula field's inserts tell the unit after each chunk; the nebula field's slab says the column of cells it has walked (`nebula_field.table_rows(on_column)`); the phenomenon scatter's indexes come off one statement each, counted. No `waiting` tag is used. The `plan-skeleton`, `clear-bright-stars`, `empty-phenomenon-scatter`, `drop-phenomenon-indexes`, `clear-scatter` and `draw-nebula-field` units go through these.
- The phases that cannot pass the floor (a few rows written, an outline of tens of thousands of cheap calls, one range of 200,000 ids, a block light roll that already counts) are left alone; the reasons are in the pull request's audit table.

generation-output: unchanged
generation-affects: none
