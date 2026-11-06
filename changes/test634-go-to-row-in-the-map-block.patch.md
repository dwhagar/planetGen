**Headline:** On the Galaxy Map the Go to row and the map are one block that fills the window, so the page rests with the row in view under the header and the map is about 40 px shorter.

### Fixed

- The Galaxy Map scrolled itself to the map when it loaded (MAP.350), which left the Go to row (the address box and its button) above the window under the held header, over the site search link, Help and Menu; the layout test `test_controls_do_not_overlap[web.galaxy]` failed at 390, 600, 820 and 1280 px for it (TEST.634). `mapfit.js` now takes a block: the row and the map together are 95 percent of the window (less the held header), the map is the row's height shorter, and the page aligns to the row's top. A row that is not shown leaves the map to fill the window alone, as before, and the sector page's map has no row.

### Changed

- `tests/js/mapfit.test.mjs` and the MAP.350 browser tests (`test_web_browser_map_window.py`, the Galaxy Map window geometry test) measure the block; the layout test is unchanged.

generation-output: unchanged
