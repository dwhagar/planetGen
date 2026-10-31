**Headline:** An object's page now lists the names and Object IDs it used to have, its JSON carries them, and a search finds an object by a name it gave up.

### Added
- **`former_designations` on the JSON of an object (GEN.1140).** `GET /api/objects/<ref>` (every kind), `GET /api/systems/<id>`, `GET /api/sectors/<id>` and `GET /api/phenomena/<type>/<id>` list what the object used to be called, newest first, each `{designation, kind, reason, retired_epoch, retired_on}` (`kind` is `name` or `uid`), `[]` for an object that never changed. The 301, 409 and 410 bodies of GEN.1139 use the same list, now with `retired_on`. Additive, so the API version stays.
- **A "Formerly" line on the item, system, sector and phenomenon pages**, under the badges: "Formerly: Old Name (name, renamed, epoch 130, 2026-10-31)". An object that never changed shows none.
- **Search by a past name.** The sector, system, star, planet and moon name searches also match a name an object gave up (only for objects that still exist), and a row found that way carries `formerly` and shows "(formerly Old Name)" after its name on the Search page.
- **`GET /api/lookup?text=` (new route).** An Object ID or a name answers `current`, or `redirect` with the object that has it now (`retired` is what was asked for, `also` lists other objects that once had the same name), `410 gone` for a deleted object's designation, `404` for one nobody had. The Search page asks it when the box holds an Object ID or a name no object has now and writes "X is now Name" above the results.

### Changed
- `asteroid_belts` leaves `designation_rows.NAMED_TABLES`: a belt has an ID and no name.
- `designations.former_lines` becomes `designations.former_json` (dicts with the date, not tuples); `designations.resolve_name` is read now, so its register line is gone.

No generator output changes.

generation-output: unchanged

TODO items resolved: GEN.1140.
