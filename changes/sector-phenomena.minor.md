### Added
- **Every kind of stellar phenomenon now appears on the Sector Map.**
  Supernova remnants, rogue planets and interstellar comets had no galaxy
  position, so they were missing from the Sector Map, the sector's
  listing and NAV. They now have one (schema v28), drawn as a glowing
  shell, a dim world and an icy coma, and each is clickable like the
  others. A supernova remnant's leftover black hole or neutron star sits
  at the remnant's center. `migrateDb.py` gives existing ones a random
  spot inside the sector they were generated in.

### Changed
- **The sector page lists its phenomena alongside its systems.** One
  "Contents" table replaces the separate Systems and Nearby Exotic
  Phenomena tables, nearest the sector's center first, with a distance
  column, 50 rows a page. A phenomenon generated as part of a sector is always listed
  there, even if an older placement put it outside the cube.
