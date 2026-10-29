**Headline:** The Galaxy Map can mark the star systems that hold a world of a chosen planet class, one optional marker for each class, off until you turn it on.

### Changed

- **Markers for chosen planet classes (MAP.220).** The Menu has a World classes group with one marker for each of the twenty planet classes, all off until chosen (`show=marker-class-K` in the address for class K). A system with a planet or moon of the class gets a disc with the class letter, hover text counting the worlds ("Fixture Eden, 2 class K worlds") and the class in its card. The marks come from the blocks of 3 shells or fewer, so they appear once the map is close enough to have a system's sector; a system a world reason that is on already marks (inhabited, habitable, gas giant) is still marked once. Migration 0172 (schema v172) fills the new `class-` rows of `system_interests` from the planets and moons already stored; nothing a generator draws changes.
