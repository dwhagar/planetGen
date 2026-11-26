**Headline:** The Kozai pump no longer leaves a planet on an orbit more eccentric than every other orbit of a planet is allowed to be.

### Fixed
- A planet a companion's or a sibling's Kozai cycle pumped could be left at the eccentricity of its place in the cycle with no limit (a class C planet at 1.06 AU of a 0.48 solar mass star at 0.96), past the 0.9 every drawn orbit is held to and the periapsis limit `orbit_elements.max_planet_eccentricity` gives. The eccentricity it takes, and the birth eccentricity the tide wears down from, are now held to that bound. A planet whose peak puts its periapsis inside the star is still lost as before.
- The two-star test of the pump read seed 26, which no longer gives a pumped planet once the draws moved; it takes the first seed that does.

generation-output: changed
generation-affects: fill
