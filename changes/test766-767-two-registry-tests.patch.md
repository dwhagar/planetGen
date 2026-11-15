**Headline:** Two registry tests that were red on main are green: the deep database check's step kind and the two job pages of a run.

### Fixed
- `test_step_registry` knows `db-check` is used (the deep check names it by the constant `DEEP_STATS_KIND`), TEST.766.
- `admin_queue_jobs` and `admin_queue_job` are named in the admin registry as pages reached from a run, TEST.767.

generation-output: unchanged
