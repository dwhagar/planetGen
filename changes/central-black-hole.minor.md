### Added

- **Schema v36: a supermassive black hole at every galaxy's center.** When
  the nucleus roll finds no quasar (90% of galaxies), a quiescent
  supermassive black hole (1e6-1e8 solar masses, like Sagittarius A*) is
  placed at the galactic center instead. It has no galactic orbit, a faint
  accretion flow far below its Eddington limit, and a sphere of influence
  of `G*M/sigma^2`.
- `black_holes.mass_class` (`stellar`, `intermediate`, `supermassive`);
  the migration fills it from each row's mass.

### Changed

- Intermediate-mass black holes now span 1e2-1e5 solar masses
  (log-uniform) instead of 100-1,000, still 2% of black holes
  (`BLACK_HOLE_INTERMEDIATE_MASS_CHANCE`), and their text says they grew
  in a dense star cluster rather than in one supernova.
