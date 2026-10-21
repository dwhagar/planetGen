**Headline:** If the core-ball batch-fill test ever counts a fourth system again, its failure now names the systems and the galaxy seed.

### Changed
- Tests: `test_core_ball_batch_fill_makes_planets` failed once with 4 systems where `--num-systems 3` asked for 3 (Core 5, on an older main). It does not come back on current main: 69 random-seed galaxies in one process and 40 separate test processes (the reproducer on the item) all gave 3. Nothing in the generator is changed on a guess; the assertion now prints the systems' ids, sectors and names and the galaxy's seed, so the next sighting carries the cause. TEST.254 stays open as "not reproduced".

TODO items resolved: none / Added: none
