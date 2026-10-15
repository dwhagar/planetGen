**Headline:** The Galaxy Map is drawn on the new shell, band and wedge grid, with the same drill-down from the galaxy to a sector.

### Changed
- The Galaxy Map is drawn on the shell, band and wedge grid: the block ladder (243, 81, 27, 9 and 3 shells deep, then sectors) nests wholly inside each block of the next size up, the whole galaxy is the level-243 blocks seen from above, and a window is the blocks inside a ball round its centre (MAP.205, MAP.206). The cylinder drawing (rings, layers, slabs and quarter stages) is gone.
- The map's data follows: the stored block totals (`galaxy_block_stats`, schema 92) are keyed by shell, band and wedge, the stage and window routes answer in them, and the tile cut covers the galaxy's ball. Stage links are `?win=m.shell.band.wedge`; the address bar reads `shell/band/wedge`, `shell N band N wedge N`, points, designations and names.
- `blockGenerateButtons` is removed with the cylinder drill-down.
- The Galaxy Map page and the sector page's map are drawn again, the "map is being rebuilt" notices are gone, and `planetgen.cli.warm_map` ends the sector-grid switch's maintenance mode (`maintenance.SECTOR_SWITCH`) once the opening view is built, so the site reopens when the new map works. An admin's own maintenance stays.

TODO items resolved: MAP.205, MAP.206 / Added: none
