**Headline:** The wheel-anchor, system-disc and breadcrumb browser tests wait for the camera and the address to stop changing before reading them, so they no longer fail when a loaded machine is slow.

### Changed
- `tests/test_web_browser_sector_reach.py`: `test_the_wheel_anchors_on_the_star_under_the_pointer_not_the_mesh_behind_it` and `test_a_systems_disc_fades_in_as_the_camera_comes_up_to_its_star` read the star's screen spot right after `galaxyCentreOn`, a camera flight. On a slow run the star had moved on by the time the anchor was asked for, so the anchor fell back to the focus plane's point (0.031 pc off the star on two axes) and the test failed alone in two of nine runs. They now read the spot with `_still_spot`: the same place on three reads 100 ms apart.

- `tests/test_web_browser_fixture_maps.py::test_galaxy_breadcrumb_always_matches_the_view_its_url_names` reads the page's URL only once it has stopped changing (`_url_settled`; TEST.126). The map writes the camera into the address 350 ms after it last moved, in its frame loop, so on a loaded run the scale bar held still before the address named where the camera ended and a fresh load of the stale address showed another sector (`Sector 1·0·2` against `Sector 1·0·1`, seen in four gates). Three runs at load 14 to 16 pass with it; the failure came in one of four runs without it.

- `test_a_systems_disc_fades_in_as_the_camera_comes_up_to_its_star` (TEST.213) stopped its wheel loop at a share of 0.999 and then asked for 1.0 within 1e-6, so it failed whenever the last step landed between the two (0.99975 in a loaded run). It now asks for the share it stopped at: at least 0.999, whole.

### Added
- Nothing.

TODO items resolved: TEST.231, TEST.213 (and the breadcrumb sighting on TEST.126; its other sightings and TEST.172 stay open) / Added: none
