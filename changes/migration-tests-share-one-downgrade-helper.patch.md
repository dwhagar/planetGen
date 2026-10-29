**Headline:** The migration tests no longer carry their own hand-written undo of every later migration; one helper takes a database back to any released version.

### Changed
- **One `downgrade_to(config, version)` replaces the per-test undo blocks.** It loads the stored `schema_v<N>` fixture into a scratch database, compares it with the test database table by table (columns, keys, constraints, views), and alters the test database to match, keeping the rows and stamping the old version. The seventeen migration tests that each listed the columns and tables later migrations had added (and had to be edited every time a migration landed) now call it; the old per-feature `drop_*` helpers and `undo_parent_link` are gone. A new test checks the result against every fixture it is asked for, and that the real migrations bring a downgraded database back (TEST.328).

TODO items resolved: TEST.328.
