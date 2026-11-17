**Headline:** A companion star now pulls on the orbits of a star's planets (their periapses, forced eccentricities, and orbit planes); a lone star's planets are unchanged.

### Changed
- **A companion star forces its star's planets (GEN.1499 slice 1).** `secular_chain.evolve` takes the companion's terms (`secular_chain.StellarTerms`, made by `belt_forcing.stellar_terms`): its `A_c` on the planets' apsidal and nodal diagonals and a push along its periapsis and its orbit plane from one extra held body at rest, so each planet circles the fixed point of the planets' own exchange plus the companion's (the forced eccentricity of GEN.1492 slice 0, now with the planets' own pulls, relativity and the belt rings in it). `planets.secular_exchange` builds them from `radiation_field.companions_of` and now also runs for a star with one planet and a companion. A companion with no known plane turns no node. Reused, not rebuilt: slice 0's `companion_rates`, `forced_planet_vectors` (the closed-form check of the tests), `secular_chain._drive`. A lone star gives today's planets to the last digit (tested). Census of 200 generated pairs (279 stars with a companion, 2,626 planets): the terms move 1,457 planets in 275 stars, the mean eccentricity of a moved star's planets goes from 0.135 to 0.138 and the mean inclination from 10.6 to 14.4 degrees; none of the exchanges turns "outside the theory" (`tests/audit/planet_forced_census.py`). No migration.

generation-affects: fill

generation-output: changed
