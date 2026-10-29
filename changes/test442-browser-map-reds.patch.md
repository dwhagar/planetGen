**Headline:** Four browser map tests that had gone red now test what they meant to, and the Galaxy Map exposes a hook that says when an opened sector can be clicked.

### Fixed
- **The Menu filters test matched two summaries (TEST.442).** The Charted regions list is a second `galaxy-sub` group in the Menu; the test now asks for the "Show on the map" group by its id.
- **The System Map drill test picked a system whose whole-system scene draws only stars (TEST.442).** A system of several stars draws its stars there and each star's planets in a scene of its own, so the test now picks a single-star system with a planet that has moons, and runs in a 1280x900 window because the planet's info card (MAP.332) covers the moons in 800x600.
- **A click on a star of an opened sector came before the sector's scene was up (TEST.442).** The panel names a sector from its address before its scene is fetched, and a click then fell on a cell of the galaxy. `canvas.galaxyOpenSector()` (a hook for the browser tests, like `galaxyAmbient`) now reports the open sector once its scene is whole, and the test waits on it.
