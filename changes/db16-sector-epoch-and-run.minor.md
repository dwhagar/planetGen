**Headline:** Each sector now records the generator epoch and the run that saved it, in place of four version text columns.

### Changed
- **`sectors.generator_epoch` and `sectors.run_id` replace `version_key`, `planetgen_version`, `python_version` and `platform` (DB.16, schema v121).** A sector keeps a small epoch number and a reference to its `generation_runs` row, which already holds the version key, release, Python and platform, so the same text is no longer repeated on every sector. A sector saved outside a run points at one shared row for the running code. The migration gives each distinct set of the old values one run row (command `migrated`), points the sectors at it and drops the old columns and their index, rebuilding the sector table twice; a sector with no recorded version keeps an empty run. The mixed-version warning and the version-key check read the key through the run.

TODO items resolved / Added: DB.16 / Added: none
