**Headline:** the galaxy map camera preset test waits for the window's turn to end instead of reading the tilt early.

### Fixed
- **`test_galaxy_map_camera_presets_at_each_zoom_step` failed about four runs in five.** The test read the camera's tilt as soon as `_settle` returned, and `_settle` waits on the scale bar, which follows only the camera's distance. The last window's preset is a turn of the view that ends a third of a second later, so the test saw 59.7 degrees (the 5 degree manual turn of the step before added to the preset's 54.7) and the camera then held 54.7 as the preset says. The presets are right; the test now waits for the tilt to reach the preset (up to 8 s) before it asserts, and a turn that really carried over would still stay put and fail.

generation-affects: none
generation-output: unchanged
