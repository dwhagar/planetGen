**Headline:** The galaxy-run timing test makes the control schema its run records its plan in, so it passes alone and in any order.

### Fixed
- OPS.297: `test_object_timing.py::test_a_galaxy_run_saves_each_sectors_timings_under_its_run` ran `plan` and `sectors` on the plain `mysql_config` database, which has no `engine_runs` table (the control schema holds it), so the run's plan rows failed with "Table ... engine_runs doesn't exist" unless an earlier test of the same worker had left the control schema in the database. It takes `control_config` (OPS.278), the fixture for a test that runs a stage. Every assertion is kept.

generation-output: unchanged
