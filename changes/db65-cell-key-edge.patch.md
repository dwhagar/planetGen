**Headline:** The cell key no longer wraps for a position within 3e-11 parsecs of the 262 kpc box edge: it takes the last cube instead of the opposite edge's.

### Changed
- `cellkey.cell_index`, `cell_keys` and the SQL expression (`sql_cell_key`) hold the cube number at 2^21 - 1 (`cellkey.LAST_CELL`). Migration 0106 recomputes the generated `cell_key` of all 14 tables that have one (a table rebuild per table); `planetgen.cli.migrate` shows its progress.
- `store.detect_schema_version` reads generated-column expressions to recognise v106.

No generated object lies that close to the reach, so no stored value changes and the generator epoch stays.

TODO items resolved: DB.65. Added: none.
