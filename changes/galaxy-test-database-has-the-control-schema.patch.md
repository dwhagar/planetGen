**Headline:** A test that runs a stage gets the job-engine tables it reads, so the plan and scatter tests that failed on main pass.

### Fixed

- OPS.278: one `control_config` fixture in `tests/conftest.py` (the leased test database as the control database, with the control schema in it and `PLANETGEN_CONTROL_DATABASE` pointed at it) replaces the two copies in `test_stages` and `test_job_runner_step_times`, and `test_galaxy_seed` (its three `plan` tests) and `test_progress_eta` use it instead of pointing the control database at a database with no `engine_runs`, `engine_jobs` or `engine_attempts` (OPS.215, OPS.216), which failed them with error 1146 on main. No product code changes.

generation-output: unchanged
generation-affects: none
