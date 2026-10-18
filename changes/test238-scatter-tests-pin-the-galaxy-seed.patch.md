**Headline:** The phenomenon scatter tests no longer fail now and then on a galaxy seed that places no hypervelocity star.

### Changed
- Tests: `test_phenomenon_scatter.py` seeds its toy galaxy from a fixed galaxy seed. The toy galaxy reaches about 36 pc, so the hypervelocity count is 0 for about one seed in six and `test_a_plan_run_stores_the_scatter_and_its_seed_and_a_rerun_replaces_it` failed on that draw (TEST.238); the file's 80 tests pass on the fixed seed.

TODO items resolved: TEST.238 / Added: none
