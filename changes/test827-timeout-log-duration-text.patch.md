**Headline:** The test of the statement-timeout line in the activity log reads the duration as the clock text the log now writes.

### Fixed
- `test_a_statement_that_hits_the_time_limit_is_named_in_the_activity_log` turned the logged `seconds` into a number, but the log writes it as clock text ("201 ms", since the duration format of UX.319), so it failed on main and on Train 60. The test now reads that text back to seconds (milliseconds, seconds, minutes, hours, days) and keeps its limit of under 4 s.
- The docstring of `_put_irregular_moons` lost a line break in PERF.288; restored.
