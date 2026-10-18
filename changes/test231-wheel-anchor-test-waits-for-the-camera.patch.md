**Headline:** The wheel-anchor and system-disc browser tests wait for the camera to stop before reading where the star is on screen, so they no longer fail when a loaded machine is slow.

### Changed
- `tests/test_web_browser_sector_reach.py`: `test_the_wheel_anchors_on_the_star_under_the_pointer_not_the_mesh_behind_it` and `test_a_systems_disc_fades_in_as_the_camera_comes_up_to_its_star` read the star's screen spot right after `galaxyCentreOn`, a camera flight. On a slow run the star had moved on by the time the anchor was asked for, so the anchor fell back to the focus plane's point (0.031 pc off the star on two axes) and the test failed alone in two of nine runs. They now read the spot with `_still_spot`: the same place on three reads 100 ms apart.

### Added
- Nothing.

TODO items resolved: TEST.231 / Added: none
