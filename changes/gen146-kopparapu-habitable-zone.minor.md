**Headline:** A star's habitable zone now follows its temperature, from the Kopparapu fits, so a cool star's zone sits where the physics puts it.

### Changed
- **The habitable zone is the Kopparapu zone (GEN.146).** `orbits.calculate_habitable_zone(luminosity, temperature)` replaces the fixed fluxes (1.1 inside, 0.53 outside) with the runaway greenhouse limit (inner edge) and the maximum greenhouse limit (outer edge) of Kopparapu et al. (2014, ApJL 787, L29, Table 1, 1 Earth mass), each a quartic in `Teff - 5780 K`. The Sun's zone is 0.95 to 1.68 AU (it was 0.95 to 1.37). A 3,200 K star's outer edge was a third too near; a hotter star's zone now lies farther out for the same luminosity. A temperature outside the fitted 2,600 to 7,200 K takes the nearest end's edges (a remnant's zone stays as small as its luminosity makes it). A close binary's zone uses its combined luminosity and the pair's temperature weighted by each star's light. The mass-dependent runaway sets (0.1 and 5 Earth masses) are in the paper and not used yet.
- A seed generates different systems from here on: the zone decides where planets are placed and which classes they take.

generation-output: changed
