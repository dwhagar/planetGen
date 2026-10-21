**Headline:** The database gets the object index table, where every list and the name search will read from.

### Changed
- **`object_index` (DB.79, schema v113).** One row per object of every kind the site lists (systems, stars, planets, moons, belts, comets, facilities and each phenomenon kind): kind and 96-bit object ID, name, sector, system, parent, class, mass and radius, with a keyset key per sort and a FULLTEXT index on the name. The migration only makes the table; it is empty until the index builder (PERF.140) fills it, and no other table is touched.

TODO items resolved / Added: DB.79 / Added: none
