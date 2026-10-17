**Headline:** The wheel zoom toward a star in an opened sector now ends at the star instead of sliding past it, never passes through a star on the way, and a double click on empty space flies half way toward the point under the pointer.

### Changed
- In an opened sector the point the wheel zooms toward is the star under the pointer (within 14 px of it on the screen; one with a system, not left off the map). Before, the picker's cell behind the star (2.3 pc away when the star was 0.28 pc away) was taken for the point, so the camera slid sideways past the star; this was the real cause of TEST.213 (`test_a_systems_disc_fades_in_as_the_camera_comes_up_to_its_star`, "the system's scene never came in"). The sector's own cell is never the point inside the opened sector; the focus plane is, when no star is there.
- A wheel step in is also limited so the eye never comes within the clearance (2% of the view's distance) of any star of the sector on its way to the point (`limitPathClearance` in `mapcontrol.js`).
- A double click on empty space (nothing picked) flies half way toward the point under the pointer, which stays where it is on the screen; on a block it is still two clicks, on a star or other object the flight to it.
- `docs/design/fly-through-view-distance.md` says what is built.

### Added
- `canvasEl.galaxyAnchorUnder` (the point the wheel zooms toward at a screen point) for the browser tests; `sectorStage.hides(entry)`.
- A node test for `limitPathClearance`, and browser tests in `tests/test_web_browser_sector_reach.py`: the anchor is the star under and a few pixels from the pointer, twelve wheel ticks bring the eye nearer the star every time with the star staying under the pointer, and a double click on empty space halves the distance.
- Not built, as a new item for TODO: the column-density ray fallback. It would only act where the cursor ray meets no cell, and the cell scene has no density field to sample there.

TODO items resolved: MAP.245 (the clearance and the half-way double click; TEST.213 fixed) / Added: the column-density ray fallback as its own item
