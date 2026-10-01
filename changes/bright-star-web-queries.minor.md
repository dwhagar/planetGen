### Added
- Bright-star queries for the web pages: `adminStats.bright_star_counts` (placed, filled and unfilled pre-placed bright stars, also in `GET /api/admin/stats` as `database.bright_stars`), `queryDb.bright_stars_in_sector` (a cell's bright stars, unfilled only by default) and an `unfilled_only` option on `queryDb.galaxy_bright_stars_in_box`.
- `GET /api/galaxy/shape` now returns `bright_stars`: whether the bright-star scatter has run, its threshold and seed, and the default threshold (`queryDb.bright_star_scatter_status`).
