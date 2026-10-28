**Headline:** The Galaxy Map can mark the systems that hold a gas giant or a Yellow-or-better world, each with a reason you can read in its infobox.

### Added
- **System interests (MAP.305, MAP.303).** A new `system_interests` table holds one row per star system and reason (`habitable`: a planet or moon scoring Yellow or better, `gas-giant`: a gas giant planet), with the count of bodies and the best habitability score. Each reason is one entry of the registry in `db/interests.py` and one in the map's `INTEREST_REASONS`, so adding a reason needs no change to the table, the writers, the query or the node. The Galaxy Map has a marker kind for each reason (off until asked for in the layers menu); a system with several reasons is drawn once and its infobox lists them all. The infobox line carries the tier of the best world ("2 habitable worlds, best Green").

### Changed
- The habitable-world marks now include Yellow worlds, not only Blue and Green ones (the Blue/Green rule for a habitable world, GEN.643, is unchanged).
- The map's node payload carries `interests` in place of `worlds` (tile format 7).
- Schema v165: `system_interests` replaces `habitable_systems` (migration 0165 fills it from the stored planets and moons; no epoch, no regeneration).
