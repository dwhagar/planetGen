**Headline:** The registry test of Object ID tables knows the two scatter tables that carry the ID of the system or phenomenon built from them, and passes again.

### Changed

- `bright_stars` and `phenomenon_scatter` gained a `uid` with GEN.1254 (the Object ID the built system or phenomenon takes over) but were on neither the `uid` group of `db/kind_tables.py` nor `point_columns.NOT_BODIES`, so `test_every_table_with_an_object_id_is_in_the_uid_group` failed on main (GEN.1301). `kind_tables.SCATTER_UID_TABLES` names them: they are not objects of their own (the `uid` group is what `api.ids` resolves) and not `NOT_BODIES` (their point is a source on `point_columns.LAYOUTS`). The test accepts exactly those two and fails if one loses its `uid`. No generated value changes.

generation-output: unchanged
