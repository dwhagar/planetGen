**Headline:** The Galaxy Map window-size test now expects the 95 percent of the browser window that the map is sized to, so its three sizes pass.

### Changed

- `test_galaxy_map_window_fills_the_page_width_and_the_height_the_window_leaves` still expected the map to be the browser window's height less the header and 40 px (the first MAP.350 rule). The map is sized to 95 percent of the window's height less the header (`mapfit.js` `fillHeight`, Boss's later rule), so the 1000x900, 1280x720 and 390x700 runs failed by 2 to 8 px. The test now expects `max(320, round((height - header) * 0.95))`. The three runs pass (TEST.575).
