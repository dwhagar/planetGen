**Headline:** What is nearby is read by the cell key, so it finds every object (also those with no sector) and its cost follows the rows near the place, not the size of the galaxy.

### Changed
- `objects_within` (the nearby page, `GET /api/near`, `planetgen query near`) reads each positioned table by its cell key: star systems, stand-alone facilities and the eight phenomenon tables, through the new `planetgen/db/region.py` (a list of at most 256 key runs joined to the table's covering index, the exact distance test in the same statement, a `UNION ALL` branch for the rows outside the key box). The `sectors` table is no longer searched; a sector is read by its id only to name the sector of a row on the page. A system or facility with no sector is found like any other and has no sector shown.
- The answer no longer reports `sectors_in_range` and `sectors_generated`. It reports `not_generated` (pre-placed bright stars and scatter objects in the sphere that are not built yet, counted by kind) and `not_generated_total`; the nearby page and the command line say how many are not generated yet.

### Added
- Tests: the same rows as brute force on the seeded galaxy, the NULL-key branch, an object with no sector, the not-generated counts, and the `EXPLAIN FORMAT=JSON` plan shape.

TODO items resolved: NAV.60 / Added: none
