**Headline:** the heartbeat rules test no longer expects a duration in raw seconds that the Queue page now says in minutes.

### Fixed

- TEST.760: `test_api_web_and_tile_prefetch_jobs_are_judged_by_the_same_rules_on_the_same_clock` looked for "600" in the quiet reason, but `beatjudge` writes durations through `format_duration_seconds` ("10 minutes"). The limit (`QUIET_FLOOR_SECONDS`, 5 minutes) and the test's clock agree; only the wording moved. The test now builds the expected text from the same formatter.
