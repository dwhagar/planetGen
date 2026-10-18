**Headline:** The Galaxy Map's drill-down and forward-to-current tests no longer read the breadcrumb while the camera is still flying.

### Changed
- Tests: `_wait_settled` in `test_web_browser_maps.py` returned after 250 ms of an unchanged address and trail, but the trail follows the camera (MAP.227) and the software-drawn map's frames can be more than 250 ms apart, so a flight still in the air left both unchanged and then moved the trail to the cell it landed in. `test_galaxy_map_drill_down_by_clicks` recorded the trail of the first cell (for example `Block 124·3·238`) and Back returned to the cell the flight landed in (`Block 128·16·257`); `test_galaxy_map_forward_to_current_jumps_to_the_newest_view` failed the same way. It failed 2 runs in 3 on an idle machine. The wait now also needs the flight to have landed, through a `galaxyFlying()` hook on the map's canvas next to the other hooks for the browser tests; with it the three Galaxy Map tests that read the trail passed 12 runs in 12.

TODO items resolved: TEST.223 / Added: none
