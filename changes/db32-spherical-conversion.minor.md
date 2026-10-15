### Changed
- Schema v91: the galaxy is converted from cylinder ring/layer/slot sectors to spherical shell/band/wedge sectors (DB.32, migration 0091). Every object keeps its position, velocity and epoch and is refiled under its new sector; planets and moons follow their star. `sectors`, `bright_stars`, `sector_stats` and `phenomenon_scatter` use `shell_index`, `band_index`, `wedge_index`; `galaxy_shell` replaces `galaxy_layer` and `galaxy_column`; `galaxy_shape.outer_ring_index` is `outer_shell_index`.
- Object IDs are renumbered from the new birth sector (body numbers are kept). Printed IDs, bookmarks, URLs and cached API v3 IDs change.
- The Galaxy Map queries of the block ladder raise until they are rebuilt on shells (MAP.205 and MAP.206).

### Removed
- `galaxy_block_stats` and `db/block_stats.py`, and the block totals refresh at the end of a `galaxy` run.
