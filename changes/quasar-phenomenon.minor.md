### Added
- **Quasars.** A galaxy's nucleus can now be active: 10% of the time
  (`QUASAR_ACTIVE_NUCLEUS_CHANCE`) the first core sector gets a quasar at
  the exact galactic center, so a galaxy has at most one. Each has a
  supermassive black hole (1e8-1e10 solar masses), a luminosity set by
  its Eddington ratio, the matching accretion rate and broad-line-region
  size, and ~10% are radio-loud with jets. It shows on the Sector Map,
  in the sector's Contents table, in the phenomena list and on its own
  detail page, and `generate.py phenomenon --type quasar` makes one on
  demand (`--sector-id` must be a shell-0 sector). New `quasars` table,
  schema v31; run `migrateDb.py` (update.sh does).

### Fixed
- **A black hole's accretion-disk temperature could render one kelvin
  low.** Loading an anchored black hole back from the database truncated
  its fractional disk temperature, so the rendered page could read e.g.
  3,676,064 K instead of 3,676,065 K.
