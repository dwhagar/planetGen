**Headline:** The second system-disc browser test no longer fails when the disc's last share lands just under a whole one.

### Changed
- `tests/test_web_browser_sector_reach.py::test_a_systems_disc_fades_in_as_the_camera_comes_up_to_its_star` stops its wheel loop at a share of 0.999 and then asked for 1.0 within 1e-6, so it failed whenever the last step landed between the two (0.99988 in a gate, failing alone as well). It now asks for the share it stopped at: at least 0.999. The sibling test (`test_the_system_disc_...`, fixed under TEST.213) already did.

### Added
- Nothing.

TODO items resolved: none (the same defect as TEST.213, left in this second test by its first fix) / Added: none
