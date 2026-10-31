**Headline:** The life paces a star supports now come from its own lifespan, not from a table by spectral letter.

### Changed
- **GEN.942: a planet is no longer asked to wait past its star's death.** The table that gave each spectral letter a lifespan and the paces of evolution it supports (`tuning.STAR_EVOLUTION`, from the old `10 Gy x M^-2.5` law) is removed. The paces a star supports are the ones whose technological civilization comes within the star's own modelled lifespan (`stellar_object.supported_evolutionary_scales`; a white dwarf uses the time it has had). A G star that lives 10 Gy supports the fast and normal paces, not the slow pace that ends at 35 Gy; A stars (model lifespan 0.65 to 2.2 Gy) support only the fast pace; M dwarfs support all three. A star's minimum age for the worlds it hosts follows, so some stars' ages and so their planets change (generation-output: changed). The star's page says what its lifespan allows in one sentence of its own in place of the per-letter note. Stand-in registry row LF3 is retired with the table.

generation-output: changed
