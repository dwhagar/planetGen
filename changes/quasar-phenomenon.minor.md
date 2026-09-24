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
  schema v30; run `migrateDb.py` (update.sh does).
