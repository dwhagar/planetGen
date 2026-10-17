**Headline:** The Galaxy Map's address holds the camera, so a copied link, a bookmark and a reload open the view exactly where it was.

### Changed
- The address is `?cam=x,y,z,d,qx,qy,qz,qw` (the view's middle and distance in parsecs and its turn, ten significant digits). It is written in place once the view has been still for a moment, a hop to another stage adds a history entry, and the sector and display parts of the query (`sector=`, `open=`, `hide=`, `color=` and the rest) are left as they are. Opening one works the stage out from the camera (the finest window of 9 cells, centred on the cell holding the view's middle, that holds the view; the whole galaxy for a view of all of it) and puts the camera there. A history entry also keeps its stage, so Back, Forward and a reload return to it exactly.
- Bookmarks of a view keep the camera.
- The `?win=` and `&g=` tokens, and the `win=` that a sector link carried, are no longer read: those links open the galaxy.

### Added
- Tests: the camera query and its refusals (`test_galaxystages.py`), the address written in place and a link that opens the camera (`galaxymap.test.mjs`), and the browser tests of the stages now open a window through the page's own hook instead of `?win=`.

TODO items resolved: MAP.228 / Added: none
