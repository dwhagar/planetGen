**Headline:** A script now draws the same few hundred seeded systems on two builds and counts how many bodies a change moved in ocean tide, tidal heating, locking, spin, tilt and eccentricity; no generated world changes.

### Changed

- GEN.1448: new `scripts/census_interactions.py` (`draw` then `compare`). It draws seeded systems of the six census star types, records for each planet and moon the spring ocean tide range, tidal heating flux, locking state, spin period, tilt and eccentricity (and any Kozai loss a build records), and a second run on another build is compared body by body; a system whose body count changed is counted apart. It is run before and after GEN.1446, GEN.1449, GEN.1454 and GEN.1450 land so each says what it did to the world.

generation-output: unchanged
generation-affects: none
