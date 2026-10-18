**Headline:** The neighbourhood test's dense-core case no longer fails by chance: it uses a fixed galaxy and checks a capped cell against the ceiling, not against a count the run does not produce.

### Changed
- Tests: `test_random_start_neighborhood_matches_the_real_skeleton_plan[dense_core]` still failed in 3 runs of 10 after the first TEST.225 fix of #2227 (a start in a cell at 1,000 times the local density). One draw in three gave a galaxy too light to hold such a cell, so the run found no start and exited with `SystemExit: 1`; the other kind counted a capped cell by `cell_resolution(...).resolved_systems` (273 expected, 2 generated; 81 expected, 1 generated), which the run does not reach, since the scatter that holds the rest is not drawn in this test. The dense-core case now plans a fixed galaxy seed that holds such cells, a capped cell is left out of the sum and checked to hold no more than the cell star ceiling (GEN.281), and the other cells keep the 3x band against the plan's prediction.

TODO items resolved: none (completes TEST.225) / Added: none
