**Headline:** The wheel on the Galaxy and Sector Map zooms toward what is under the pointer, which stays under it.

### Changed
- The wheel zooms toward the star, remnant or cell under the pointer (else the point of the focus plane under it): the center moves toward it so it stays where it is on the screen. A step in moves the eye by at most half the way to it (`limitApproach`, `targetAboutAnchor` in `mapcontrol.js`), so repeated steps approach it without landing on it; a view that cannot be moved, or a wheel with no pointer over the map, zooms about the middle as before. The Shift-click, the zoom buttons and the slider are unchanged.

### Added
- Tests: `limitApproach` and `targetAboutAnchor` in `mapcontrol.test.mjs`, and a browser test that the star under the pointer stays under it while zooming in.

TODO items resolved: none (MAP.225 stays open: the wheel is built; the column-density ray fallback, the clearance from the nearest pickable thing and the half-way double click are not) / Added: none
