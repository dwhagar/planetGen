**Headline:** The old-schema shape test starts from the baseline and the newest three stored schemas instead of all 137, about 1,350 s less on a full run.

### Changed

- `tests/test_db_old_schemas.py`: `test_old_schema_migrates_to_the_fresh_shape` runs for the baseline (v130) and the `NEWEST` (3) versions before the current one (`checked_versions()`), not for every stored schema. The run from the baseline passes through every revision and the fixtures hold schema only, so the start points between add nothing; the newest revisions are the ones Boss's server has not run, so they are also started from the schema just before each. The fixture-presence test now asks for those four.
- `tests/fixtures/old_schemas/`: the 84 fixtures no test loads by number are removed with their README lines (the 53 that remain are the baseline, the newest three, and the ones `test_downgrade_to` and the `test_migration_*` files still load; they go with those tests under TEST.731 and TEST.733).

generation-output: unchanged
generation-affects: none
