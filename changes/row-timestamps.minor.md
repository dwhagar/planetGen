### Added
- **Created and modified times on the main tables.** `sectors`,
  `star_systems` and the seven exotic-phenomenon tables now have
  `created_at` and `modified_at` columns, with an index on `modified_at`
  (schema v27). MySQL keeps `modified_at` current on every edit. Changing
  a planet or moon bumps its system's `modified_at` instead of the child
  row getting its own timestamp, and deleting a system bumps its sector's.
  Orbit updates from `updateOrbits.py` don't count as a modification. The
  migration adds the columns without rebuilding the tables where MySQL
  supports it, and builds the indexes online, so it's safe to run on a
  large database. Existing systems and sectors are backfilled from the
  systems' original creation times, in small batches; existing
  phenomena, which have no record of when they were made, get the
  migration's own time.
