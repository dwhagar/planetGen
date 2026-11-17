**Headline:** The sky draw of the bright-star backfill is about a fifth faster per cell: sectors that take one level or mass are written in one statement, and a cell's center is worked out once.

### Changed
- PERF.275: a sky share wrote one `UPDATE` of `sector_stats` for every cell it visited (21,000 statements for 22,000 cells in a profile), and worked out the center of the same cell again for every tier. The cells that take one level, or one mass and level, are now written 500 to a statement (`store.set_sector_bright_levels`, `set_sector_bright_masses`), and `spherical.sector_position_pc` is cached. The same rows are stored and the same stars are drawn. On a toy galaxy, two luminosity tiers over 22,402 cells took 10.1 s against 8.2 s (12.1 s against 15.0 s on a cold first pass), 19 percent less.

generation-output: unchanged
