**Headline:** The neighbourhood test no longer fails when the random start lands in the dense core, where cells hold fewer systems than the plan's raw prediction by design.

### Changed
- Tests: `test_random_start_neighborhood_matches_the_real_skeleton_plan` counted a cell above the cell star ceiling at the plan's raw prediction (up to 250,000 systems), but such a cell is filled only down to its luminosity cap (GEN.281) and resolved about 100; the sum missed the 3x band in about one dense-core start in three (TEST.225). A capped cell now counts by `cell_resolution(...).resolved_systems`, and a second case starts in the dense core (`--min-start-density 1000`) so the case is covered on every run. The 3x bound is unchanged.
- Tests: the window-address pattern in `test_web_browser_maps.py` allowed a negative number in the wedge field, but the address is `win=size.shell.band.wedge` and it is the band that goes negative below the equator; it now matches `-?` in the band, as `test_web_browser_fixture_maps.py` already does (TEST.236).

TODO items resolved: TEST.225, TEST.236 / Added: none
