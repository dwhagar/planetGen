**Headline:** Four load-sensitive tests no longer fail for reasons that are not in the code they check.

### Changed
- Tests: `test_whole_run_is_fast` measures the process's CPU time, not the clock, so a loaded machine cannot fail it (TEST.170).
- Tests: the shared bright-star scatter galaxy is seeded with a fixed galaxy seed, so the forced-scatter and backfill tests cannot draw a galaxy whose first wedges hold no star (TEST.171, TEST.173).
- Tests: the Galaxy Map click waits up to 90 s for the software-drawn map to answer instead of 15 s, which a loaded machine passed and then reported "did nothing" (TEST.172, TEST.202).
- Tests: the process-race test waits until every starter has imported before it lets them go, so a starter that came more than 15 s late on a loaded machine no longer finds the first job's lock stale and starts a second (TEST.134).
- Tests: the clean-up of the test Redis leaves a queue whole when it lists a job made in the last 15 minutes, so a second run starting beside a first can no longer delete the first's queued jobs (TEST.203).
- Tests: the page accessibility test no longer counts a request the page cancelled itself (ERR_ABORTED, the map dropping its tile fetch) as a failed request, as the settings browser test already did (TEST.174, TEST.197).

TODO items resolved: TEST.170, TEST.171, TEST.173, TEST.172, TEST.202, TEST.134, TEST.203, TEST.174, TEST.197, TEST.129, TEST.130 (the last two were already fixed on main) / Added: none
