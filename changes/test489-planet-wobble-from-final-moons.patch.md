**Headline:** A planet's wobble from its moons is made from where the moons end up, not where they were first placed.

### Fixed
- **The reflex offset of a planet with moons (TEST.489).** It was made when the moons were, and their positions moved afterwards (the orbits settle and the bodies are made again at the star's final age), so on about one planet in five the stored `reflex_offset_*` disagreed with the final positions of its moons by up to a percent. `StarSystem` now makes it again from the final moons, in the pass that makes the stars' own wobble, so `test_moon_having_planet_reflex_offset_matches_manual_formula` passes (0 of 104 planets with moons off, from 24 of 117). The moons themselves are as they were.
- generation-output: changed

TODO items resolved: TEST.489 (the reflex offset test).
