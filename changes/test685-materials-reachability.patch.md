**Headline:** A test now fails when a row of the materials table is read by no draw and has no reason for it.

### Added
- **Materials reachability test (TEST.685).** `tests/test_materials_reachability.py` walks the materials table against the modules that use it (the ground recipes, the small-body components, the compaction and bulk-structure readers and `tuning.py`) and fails on a row that no draw names and that is not tagged vocabulary-only, on a tag that a draw has since made stale, and on a new double (the same formula as a powder or grain-size variant). The 36 rows nothing reads today (the high-pressure silicate phases, three ice phases, the envelope and element rows and a few minerals) are tagged with why in `VOCABULARY_ONLY` (`scripts/materials_data.py`, carried in `materials.json` and `materials.VOCABULARY_ONLY`); the loader refuses a tag that names no row or gives no reason.

`generation-output: unchanged`
`TODO items resolved: TEST.685`
