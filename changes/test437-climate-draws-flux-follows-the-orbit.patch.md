**Headline:** The climate draws test measures the starlight the way the generator does, so it passes again.

### Fixed
- **The test that checks a planet's stored temperature against its stored draws read the light of a circle (TEST.437).** The generator has used the light of an eccentric orbit, `a (1 - e^2)^(1/4)`, since GEN.727, so the test saw a temperature 0.03 K off for every class. The test now uses `orbit_elements.flux_distance_au`.
