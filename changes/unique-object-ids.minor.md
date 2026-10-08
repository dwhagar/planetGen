### Added

- Every sector, star system, star, planet, moon, belt, comet and phenomenon
  has a unique ID (GEN.69, schema v58): a `uid` column on each table, written
  when a sector, system or phenomenon is saved. A sector's is its designation,
  so a sector nobody has generated already has one; a system's or
  phenomenon's is 96 bits and a star's, planet's, moon's, belt's or comet's 64
  bits (unique under its system), hashed from the galaxy seed, the parent's ID
  and the object's slot in it, so saving the same sector again gives the same
  IDs. An interstellar object or bright-sweep system keeps its position ID.
  Rows saved before this have none until their sector or system is saved
  again. The design is in `docs/design/object-ids.md`.
