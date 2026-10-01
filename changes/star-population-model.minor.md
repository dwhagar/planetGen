### Changed
- **Random stars are now drawn as physics: a mass, an age, then how the
  star has evolved.** Masses come from the Kroupa (2001) initial mass
  function and ages from a 0-10 Gy star-formation history; each star is
  main sequence, subgiant, giant, supergiant or white dwarf according to
  how far its age is into its own lifetime, and its spectral letter and
  Yerkes class follow from its temperature and luminosity. Every O, B and
  A star used to be a supergiant or subgiant and no white dwarfs or red
  giants were ever made; now about 74% of stars are M dwarfs, 14% K,
  5% white dwarfs, a quarter of a percent giants, and supergiants are
  about one in a million, as in the real galaxy. `+large_star` uses the
  same model above 1.4 Msun. A specified `--star-type` draws its
  luminosity log-uniformly within its class instead of linearly.
- **A binary's secondary is born with its primary.** It takes a mass
  ratio of 0.1-1.0 of the primary's initial mass and the primary's age,
  and evolves by the same model, so it's never brighter or hotter than
  its mass allows.
- **Planets follow their star's age and history.** A star younger than
  10 Myr keeps belts only; no habitable-class world or moon forms around
  a star younger than 0.1 Gy; a giant has engulfed everything inside
  twice its radius, and a white dwarf's progenitor everything inside
  1.5 AU. A required habitable world gets a star old enough to host one.

### Added
- **Stellar populations by position.** `galaxyDensity.population_densities`
  splits a point's star density into young (0-0.1 Gy, a thin disk that
  crowds the spiral arms), intermediate (0.1-3 Gy), old (3-10 Gy, the
  thick disk) and bulge (8-12 Gy) stars, summing to the same total as
  before. A system config's new `POPULATION` draws its star's age from
  one, and `MAX_STAR_LUMINOSITY_SOL` keeps it dimmer than a threshold.
- **Bright-star sampling for pre-placement.** `stellarPopulation` gives
  the share of a population's stars at or above a luminosity
  (`bright_star_fraction`), draws stars conditional on being that bright
  (`sample_bright_stars`, exact and about 10 microseconds a star) or
  dimmer (`sample_dim_star`), and `Star.from_params` /
  `StarSystem(primary_star_params=...)` rebuild a stored star without
  re-rolling it. `SpaceSector.add_preplaced_system` places one at its
  stored position before the rest of the sector fills around it.
