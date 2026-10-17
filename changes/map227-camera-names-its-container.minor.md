**Headline:** The Galaxy Map's breadcrumb and a new "In Sector …" line follow the cell the camera is in, not the stage the map was opened at.

### Changed
- The breadcrumb is the trail of containers of the cell holding the view's middle (`stageAtPoint` in `galaxystages.js`, the cells are closed form), so moving the center into the next cell of the stage's size changes its last crumb; a crumb above it still flies out to that container. An opened sector or system keeps the stage it was opened from.
- A label over the map's top-left corner (`galaxymap3d-where`) names the sector the view's middle is in (`sectorAtPoint`), updated as the camera moves.

### Added
- Tests: `stageAtPoint` and `sectorAtPoint` through the stage view in `galaxymap.test.mjs`, and a browser test that moves the center into the next sector and checks the line, every crumb and the crumb out to the galaxy.

TODO items resolved: MAP.227 / Added: none
