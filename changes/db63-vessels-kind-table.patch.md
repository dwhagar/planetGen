**Headline:** The ships and stations table is in the registry of kind tables, so an object ID can find a ship.

### Fixed
- **`vessels` (DB.63) has an Object ID column but was not in the kind-table registry (DB.74).** Main failed `test_every_table_with_an_object_id_is_in_the_uid_group`. The table is now a `vessel` kind in the `uid` group and in the `ids` order, so the public ID lookup also probes it.
