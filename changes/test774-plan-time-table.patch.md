**Headline:** the plan timing test reads the durations the table prints, so it is green again.

### Fixed

- TEST.774: `test_the_plan_says_where_its_time_went_step_by_step` parsed the plan's timing table as plain seconds, and the table prints durations in the clock form (`1h 23m 45s`, `850 ms`) since UX.319. The test now reads each row and the named total through a small helper that is the inverse of `format_duration_seconds`, with a tolerance of half a second per whole-second row. All the checks stay (sorted, at least eight rows, the sum equals the named total, the steps named).

generation-output: unchanged
