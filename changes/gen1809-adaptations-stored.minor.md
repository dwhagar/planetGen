**Headline:** The adaptations the life of a world needed are now kept as rows in new tables `planet_adaptations` and `moon_adaptations` (nothing reads them back yet)

### Added

- GEN.1809: the record of GEN.1808 (`planet.life_adaptations`, a list of dicts: life_class, axis, habitat, demand, grade, rows; a body with life of several classes has the entries of each, one class after the other) is kept with the body in two new tables, one row for each need (`entry_number`, `life_class`, `axis`, `habitat`, `demand`, `grade`, `row_ids`), written with the body, loaded back with it (`store.body_adaptations`) and kept by an edit. A need no grade of the table reaches is a row with no grade and no row ids (the life there is past the table); a habitat row has no axis and no demand. Indexed on the axis and grade, so "the worlds whose life needed the extreme grade of cold tolerance" is a short read. A rogue planet has none. Migration (a placeholder number), schema, models, column units, fingerprint tables and sections, command touches and docs. Nothing reads the rows back yet.

generation-output: changed
generation-affects: fill
