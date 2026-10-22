**Headline:** The precision table lists the ships and stations table, so its floor body is checked like the other galaxy-frame tables.

### Fixed
- **`vessels` (DB.63) had no row in the precision test's floor list (DB.77).** Both landed within the hour, so main failed `test_the_floor_of_every_table_with_a_position_is_the_documented_one`. The table is a galaxy-frame point in parsecs like the other free bodies, with the same 9,000 km floor, and the design table in unified-unit-quantity.md names it.
