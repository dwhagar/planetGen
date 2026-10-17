**Headline:** UX.114 is retired after its merge, and the round-trip column failure TEST.217 is closed by the same pull request.

### Changed

- UX.114 (members and orbits of an uncharted scatter row, #1989) is retired.
- TEST.217 is retired: #1989 put `moons.evolutionary_speed` on the list of columns that are NULL unless the moon bears life, and the failing order (`test_workers`, `test_update_guard_work`, `test_api_jobs`, `test_db_column_round_trip`) passes on main with it (67 passed). TEST.216 records that it did not recur in the second UX.114 gate.

TODO items resolved: UX.114, TEST.217 / Added: none
