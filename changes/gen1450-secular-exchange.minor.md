**Headline:** The planets of a star now pull on each other's orbits: their eccentricities and inclinations take the Laplace-Lagrange secular exchange to the star's age.

### Changed

- A star with two or more planets runs the first-order secular theory (`physics/secular_chain.py`) over all of them at once, with general relativity on the inner ones, in the pass that settles the orbits at the star's final age. Each planet's eccentricity, longitude of periapsis, inclination and node move from what the generator drew to where its siblings leave them; a heavy giant pumps the eccentricity of the planets beside it and lowers its own a little, and a chain shares its angular momentum deficit. A planet alone, and a star with a single planet, are unchanged.
- The spacing the generator keeps between neighbours is not broken: an eccentricity the exchange would raise is raised only as far as the spacing allows (a chain packed to the limit keeps its eccentricities and takes the inclinations and turned periapses), and a planet beyond the theory's limits (eccentricity over 0.3 or inclination over 30 degrees) is held where it is and pulls on the others as it is. The moons' stable range is checked against the new eccentricity.
- On 120 seeded systems the 56 stars with several planets change 644 planets' orbits; the mean eccentricity of the rocky planets in them falls from 0.048 to 0.039 (the exchange moves the angular momentum deficit of the hot draws to the neighbours and the spacing holds back the rise). Stand-ins SC1 to SC5 are registered (section 6.244).
- `secular_chain.precession_rates` gives each planet's own apsidal and nodal rates for the Kozai screen of GEN.1447 and the spin-orbit drift of GEN.1454.

generation-output: changed
