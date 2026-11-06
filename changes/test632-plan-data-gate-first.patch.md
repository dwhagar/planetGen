**Headline:** `planetgen plan` without the downloaded star data is now refused by the data gate with its usual message, instead of failing with a missing-isochrone trace.

### Fixed
- Checking the arguments of the `plan` command draws the galaxy's shape, which reads the MIST isochrones (the mean mass of a system), and that check ran before the data gate in `main`. With the data missing the command died with `DataMissing` and a traceback; the gate's message ("isochrones (SC1) ... nothing was generated") only showed when an earlier call in the same process had cached the mean masses, which is why `test_datasets_gate` and `test_disregard_missing_data` failed with `-n 3` and passed in a plain run.
- `generate.refuse_if_data_missing` now holds the gate (the `--disregard-missing-data` setting, the refusal of a bulk run, the notice lines of a run that goes on), `main` calls it as before, and the `plan` validation calls it before it draws the shape. Under `--disregard-missing-data` the shape guard is skipped when the isochrones are missing, since it cannot be made; the run itself says what it lacks. A bulk run without the option still refuses and writes nothing.

### Added
- `test_a_plan_in_a_new_process_is_refused_by_the_gate_and_not_by_a_missing_isochrone_trace` (a fresh process, no cache; it failed before the fix).

generation-output: unchanged
TODO items resolved: TEST.632
