**Headline:** The Galaxy Map's block light totals are now kept as the scatters write, so the last stage of a scatter rolls them up from small stored parts instead of reading every star and phenomenon again.

### Added

- `galaxy_block_light_parts` table (migration 0236): the light sums of each sector cell that holds an unbuilt bright star or scattered phenomenon, added to by `store.insert_bright_stars` and `store.insert_phenomenon_scatter` in the transaction that writes the rows.
- `db.block_lights.roll`, `parts_current`, `sync_cells`, `forget_parts`: the block light stage rolls the parts up (the same rows a rebuild from the objects gives) and checks them against the scatter first, making them again from the objects when they differ (a galaxy scattered before v236, a prune).

### Changed

- The `refresh-block-lights` stage no longer reads the star and phenomenon tables up to four times: the totals and the uncharted sector index both come from the parts. `fill-map-luminosity` puts right the parts of the cells it fills, and `planetgen fill-block-totals` still makes everything again from the objects, as the repair and the proof.

generation-output: unchanged
generation-affects: none
