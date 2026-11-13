**Headline:** The ten heaviest round bodies of every asteroid belt are now planets, made by the same planet generator as every other planet.

### Added
- **Belt bodies become planets (GEN.1510).** The ten heaviest round members of each belt (`tuning.BELT_PLANETS_PER_BELT`) are run through the planet generator with the member's orbit, mass and composition (rock, metal, water and ice shares) in place of the planet's own draw, so a dwarf world has an air, an interior, a surface, a day, a field, a dose, a class and habitability scores like any other planet. They take their place in distance order among the star's planets (the numerals of the planets after them shift), the belt keeps the rest, and the promoted mass comes out of the belt's mass. `planets.belt_order` and `planets.belt_rank` (schema v244) keep the belt and the rank so the belt designation still resolves; a promoted body does not count against the orbit spacing of the slot planets. A body the generator gives no class stays in the belt.
generation-output: changed
generation-affects: fill

### Changed
- Planet names after a promoted body shift by one numeral. A system with a fixed number of orbits (`SLOTS`, `NUM_ORBITS`) promotes none.
