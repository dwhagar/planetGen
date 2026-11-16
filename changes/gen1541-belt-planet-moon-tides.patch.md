**Headline:** The moons of a belt member that became a planet (GEN.1510) are checked to feel the tide of their own dwarf planet and the star the way any planet's moons do; no generated value changes.

### Changed
- GEN.1541 (GEN.1464 part C): `tests/test_belt_planet_moon_tides.py` (new, 4 tests, 366 promoted planets in 60 seeds have moons): a promoted planet is raised on by the star and its sibling planets, its moon by the planet and the star, and the promoted body is not counted again in its belt. Promoted planets are in the star's list of planets, so `physics/tide_sources` already finds their sources and the sources of their moons: no code change was needed. What is left open is only the small unpromoted belt members (asteroid binaries), which have no satellite system yet.

generation-affects: none
