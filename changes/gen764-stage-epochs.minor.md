**Headline:** A galaxy now records which epoch each scatter, each sector's contents and each species was made under, so a later run can tell which parts are behind the code.

### Changed
- **Stage epochs are recorded (GEN.764, step 2).** Migration 0164 adds the `scatter_epochs` table (one row per kind of scatter: mass, luminosity, band, phenomena, nebula-field, with the epoch it was made under) and `species.generator_epoch`. The scatter passes write their row as they finish and the clears forget it; the population pass stamps each species with the running epoch. The plan (`galaxy_shape`) and the sectors already recorded theirs.
- **An existing galaxy gets the epochs it can.** The migration gives each kind of scatter that has rows the plan's epoch, and each species the epoch of the sector its homeworld's system is in; it prints nothing it cannot know. Both new records are left out of the fingerprint.
- **`planetgen.galaxy.epoch_record`** reads these records and, for an epoch, says which stages are behind (`behind(conn)`), using the stages each epoch affected (step 1). Nothing redoes a stage yet; that is the `bring-to-epoch` command, next.
- The golden galaxy is read as written under schema 160 (`tests/golden_seed.py`, `GOLDEN_SCHEMA_VERSION`): a system row records the schema version that wrote it and the fingerprint reads it, so every migration would have moved the golden digests; now only a change of what the generator makes does.
- `tests/test_scatter_epochs.py` and `tests/test_migration_stage_epochs.py` cover the records, the migration from v155 and the stale-stage report.

TODO items resolved: none (GEN.764 continues). Added: none.
