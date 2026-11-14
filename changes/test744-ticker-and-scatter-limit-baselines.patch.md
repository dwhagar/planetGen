**Headline:** Two tests that were red on main follow the changes they missed: the one progress cadence of the layered cleanup, and the 18,000 solar luminosity floor of GEN.1681.

### Fixed
- `test_ticker::test_no_place_keeps_its_own_progress_clock` found `block_lights.PROGRESS_SAY_SECONDS` still there. OPS.185 (#6440) added it for the block-light loops; the layered cleanup (#6454), which made every place say where it is at the one `tuning.PROGRESS_LOG_SECONDS`, took away the code that read it but not the constant. It is deleted.
- `test_db_scatter_limits::test_the_core_follows_the_baseline_times_a_tenth_until_set_and_the_halo_has_its_own_floor` expected a core floor of 9,900 solar luminosities (9,000 times 1.1). GEN.1681 (Train 26, Boss: twice as high) raised the default baseline floor to 18,000 on purpose, so the core floor is 19,800; the test missed the change when the train's tests followed it. The expected number now follows `BRIGHT_STAR_MIN_LUMINOSITY_SOL` times `CORE_TIER_FACTOR`, so the next change of the default does not break it again. No product code reads the wrong number: the code was right and the test was stale. Part of TEST.744.

generation-output: unchanged
generation-affects: none
