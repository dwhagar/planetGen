**Headline:** The three test files that read the small "rich" galaxy now share one build per test run instead of building it each, cutting about 12 minutes to under 5 when they run together; no generated output changes.

### Changed
- **One rich galaxy per test process (TEST.417).** `test_db_column_round_trip.py`, `test_db_check_constraints.py` and `test_cellkey_motion.py` each built the same seeded galaxy (about four minutes: three filled sectors, system shapes, phenomena, two orbit updates). The `rich_galaxy` fixture now lives in `tests/conftest.py` at session scope (the builder is `scratch_rich_galaxy` in `tests/rich_galaxy_support.py`), so every module in a process shares one build. Nothing is dropped from the galaxy: a cheaper build (no deepest-floor draw, one sector) was tried and lost the rogue planets and ocean columns the column checks need, so the galaxy is unchanged and the three modules only read it.

generation-output: unchanged

TODO items resolved: none (TEST.417 continues).
