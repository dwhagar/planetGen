**Headline:** The habitable-systems table is on the lists the schema tests keep, so three tests that failed since MAP.220 pass.

### Fixed
- **`test_every_table_holding_an_object_is_on_exactly_one_list` and `test_no_table_inside_a_system_stores_an_absolute_position` failed on main (TEST.343).** `habitable_systems` (MAP.220) holds a copy of a system's galactic point for the Galaxy Map tiles and is not an object, so it is now on `point_columns.NOT_BODIES` with its reason, next to `object_index`.
- **`test_sql_portability` failed on the `VALUES(habitable_worlds)` in migration 0148 (TEST.343).** The migration is shipped and stays as it is; its one use is on the lint's list of allowed `VALUES()` uses, with the others.
TODO items resolved: TEST.343. Added: none.
