**Headline:** An object that is deleted keeps its ID and name on record, so a link to it answers "gone" (410) and its ID is never issued again.

### Added

- GEN.1141: `designation_rows.tombstone` and `tombstone_tree` (migration 0265, schema v265: `designations.reason` takes `deleted`). When `DELETE /api/systems/<id>`, `db/edits.delete_phenomenon` or a system edit deletes an object, its current ID (and name) is written as a retired designation and every row of the object is stamped `gone_epoch`; a deleted system takes its bodies' IDs with it and a deleted planet its moons'. The 410 answers of GEN.1139 already read the stamp; nothing wrote it until now. A sector emptied for regeneration is not tombstoned (the cell is made again). Reuse: `body_rows`, `retire_name` and the `designations` table and rules of GEN.1138. The three captured-body functions of `object_uid` (`captured_member`, `outside_reservoir_member`, `parse_captured`) are marked by design in the wiring register: they wait for the captured-asteroid-moon columns; irregular moons print through `irregular_designation`.

generation-output: unchanged
generation-affects: none
