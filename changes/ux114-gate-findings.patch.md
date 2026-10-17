**Headline:** Two test failures from the UX.114 gate are filed as bugs.

### Changed

Core 4's full gate of UX.114 failed twelve tests; ten were its own and are fixed in its branch, and two are not: `test_systems_inside_a_nebula_point_at_it_and_the_innermost_wins` and `test_every_column_is_written` (the `moons.evolutionary_speed` column), both of which pass alone on main and failed only in the full parallel run. They are filed under the test-failures goal for the Bug Hunt.

TODO items resolved: none / Added: TEST.216, TEST.217
