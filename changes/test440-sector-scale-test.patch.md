**Headline:** The Sector Map arrow-key test waits for the opening flight to land before it reads the scale, so it passes again.

### Changed

- **`test_sector_map_arrow_keys_turn_the_view_and_keep_the_scale` failed on main (TEST.440).** Since the zoom redesign Z0 (MAP.333) a sector opens with a fly-in that lands after the page is ready, so the test read the scale (10 pc) before the landing and after the first key (1 pc) and called it a change. It now waits for the scale to hold still first. The map was right; the test was early.
