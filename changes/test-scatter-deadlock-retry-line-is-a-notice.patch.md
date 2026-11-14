**Headline:** The deadlock retry test listens where the retry line is logged now.

### Fixed
- `test_scatter_write_deadlock::test_a_deadlock_is_retried_with_a_logged_wait` collected the retry line through `log.normal`. OPS.167 (#6279, warnings at the right level) moved that line, "MySQL error 1213 on try 1 of 8; trying again in ...", to the `generate` component logger's `notice`, on purpose: a retried deadlock is an expected, notable event, not ordinary progress. The test now collects `Logger.notice` as well as `log.normal` (the "asking for it" and "granted exclusive use" lines are still normal lines). The code was right; the test was stale.

generation-output: unchanged
generation-affects: none
