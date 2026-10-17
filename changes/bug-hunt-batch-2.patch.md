**Headline:** Four load-sensitive tests no longer fail for reasons that are not in the code they check.

### Changed
- Tests: `test_whole_run_is_fast` measures the process's CPU time, not the clock, so a loaded machine cannot fail it (TEST.170).
- Tests: the shared bright-star scatter galaxy is seeded with a fixed galaxy seed, so the forced-scatter and backfill tests cannot draw a galaxy whose first wedges hold no star (TEST.171, TEST.173).
- Tests: the Galaxy Map click waits up to 90 s for the software-drawn map to answer instead of 15 s, which a loaded machine passed and then reported "did nothing" (TEST.172, TEST.202).

TODO items resolved: TEST.170, TEST.171, TEST.173, TEST.172, TEST.202 / Added: none
