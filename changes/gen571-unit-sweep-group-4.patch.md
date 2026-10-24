**Headline:** The database layer and the edit pages convert AU, km, years and parsecs through the one unit registry, with the same numbers.

### Changed
- `db/store.py` (the 67 AU-to-km and km-to-AU conversions of star, pair, comet and binary rows, the orbit clock and the SQL year constants), `db/query.py`, `db/edits.py`, `db/precision.py`, `db/orbit_trees.py`, `db/point_columns.py`, `db/sector_paths.py`, `db/column_units.py` and `admin/edits.py` call `unit_registry.converter` instead of a typed factor (GEN.571, fourth file group). Every replaced expression gives the same float bit for bit, so the stored values and the SQL text are unchanged. `db/objectindex.KM_PER_LY` is left for the epoch-changing pass: the registry's light-year-to-km factor differs from it in the last digit.
