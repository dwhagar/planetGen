**Headline:** A cell that expects more systems than the cell star ceiling is filled only with the stars at or above the fill luminosity cap, and the sector page and API report what is left unresolved.

### Added
- **The cell star ceiling and the fill luminosity cap (GEN.281).** Two settings of the galaxy, stored in `galaxy_shape` (migration 0098, schema v98) and set with `planetgen plan --cell-star-ceiling` (default 50,000 systems) and `--fill-luminosity-cap` (default 100 solar luminosities). A cell expecting more systems than the ceiling is filled only with the stars at or above the cap; ordinary space never reaches the ceiling and is filled as before.
- **Unresolved systems.** `planetgen/galaxy/resolution.py` is the one rule the fill, the API and the sector page read. A capped sector's page and `/api/sector` say how many systems and how much mass are left unresolved, worked out from the address and never stored.

### Changed
- **A capped cell's stats record what it resolves.** `sector_stats.expected_systems` is the resolved count, so the actual-to-expected ratio stays meaningful, and `level_before_fill_sol` holds the cap's level.

TODO items resolved: GEN.281. Added: none.
