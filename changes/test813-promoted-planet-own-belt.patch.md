**Headline:** The promoted-planet belt test looks in the planet's own belt only, as a member's rank belongs to its belt.

### Fixed
- TEST.813: `test_the_promoted_planet_is_not_counted_again_in_its_belt` looked for the promoted planet's `belt_rank` among the members of every belt of the system. A rank is a member's place in its own belt, so a second belt (a system with an inner and an outer belt) holds a member of the same rank and the test failed on seeds 0 and 6 with nothing wrong in the generator: the planet is in no sibling list of its own tide sources and has left its own belt's members. It now looks in the belt the planet came from (`irregular_moons.belts_of(...)[planet.belt_order]`).

generation-output: unchanged
