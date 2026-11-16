**Headline:** The Sectors list no longer stores a Quadrant for each sector: the column and its index are dropped, with the `position` sort and the Quadrant filter that had no page left to use them.

### Removed

- DB.112 (schema v264, migration 0264): `sector_system_counts.quadrant` and `idx_sector_system_counts_quadrant` are dropped (the index online, the column instantly; a stopped run can be started again). The Quadrants were taken off the site by UX.272, so nothing filtered by them; `query.list_sectors`, `count_sectors` and `virtual_cells` lose their `quadrants` argument, `SECTOR_SORTS` loses `position` (so `GET /api/sectors?sort=position` is a 400 like any unknown sort), and `refresh_sector_system_counts` no longer works out a Quadrant for every sector. `spherical.sector_quadrant`, which only that column and filter used, is gone with `VIRTUAL_SCAN_ROWS` (the walk of the uncharted cells that filled a page the filter had thinned) and the Quadrant SQL. `star_systems.quadrant` and the octant `quadrant` columns are a different thing and stay. Reuse: the migration uses `helpers.drop_index_online` and `helpers.drop_column_instant`.

generation-output: unchanged
generation-affects: none
