**Headline:** With the cursor over empty space the wheel now zooms toward where the stars thicken along the cursor's ray, not toward a flat plane through the middle of the view.

### Changed
- Where the cursor ray meets nothing pickable, the point the wheel zooms toward (and a double click on empty space flies half way to) is the first place along the ray where the column of the galaxy's stars passes half the view's distance in local-average density (`columnAnchor` in `galaxystageview.js`, `columnDensityPoint` in `mapcontrol.js`). Before, it was the point of the focus plane under the cursor, however empty or dense the space along the ray. In a dense region (a bulge, an arm, the core) the point is nearer, in space of the average density it is half way to the focus, and where the ray finds too little it is the focus plane's point as before. A point is never nearer than a quarter of the way to the focus plane, so a wheel step still moves the eye.
- The density is the page's own galaxy model (`relativeDensity` in `galaxyprisms.js`), so no density field was added to the cell scene or the tiles. Not used in an opened system or before the galaxy has a shape.
- `docs/design/fly-through-view-distance.md` says what is built.

### Added
- A node test for `columnDensityPoint` (uniform space, a dense slab after empty space, a column that never passes, no density) and a browser test in `tests/test_web_browser_sector_reach.py`: with a sector open near the core, the anchor under two empty corners of the map lies between a fifth and a third of the view's distance from the eye, not on the focus plane.

TODO items resolved: MAP.247 / Added: none
