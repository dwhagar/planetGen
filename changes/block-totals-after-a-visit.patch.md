**Headline:** A sector the Galaxy Map generates when you visit it now counts in the block totals at once.

### Fixed
- **A sector generated on a visit, or generated again from its page, never reached the Galaxy Map's block totals (DB.83).** Only a `galaxy` run, a fill and a deletion kept `galaxy_block_stats` and `galaxy_block_lights` up to date, so cells above a sector showed stale or missing totals until `planetgen fill-block-totals` ran. `ensure_sector_generated` and `regenerate_sector_cell` now bring the blocks above the sector up to date after they save it; a failure only warns.
