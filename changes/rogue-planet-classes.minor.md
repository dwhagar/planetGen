### Added

- Rogue planets have a planet class (GEN.8). Each planet class now says
  whether a rogue planet can have it (a new `"r"` zone flag in
  `PLANET_CLASSES`): dead worlds (C), icy bodies (D), ice giants (I), gas
  giants (J) and gas dwarfs (T). A rogue draws its class from those that fit
  its type, radius and mass; a brown dwarf has none. The class shows on the
  rogue's page (linked to the class page), in the sector's Contents and on
  the Sector Map, and the class reference lists "Interstellar space" as a
  zone. Schema v47 adds `rogue_planets.planet_class` and gives stored rogue
  planets their class, so run `update.sh` (or `migrateDb.py`).
