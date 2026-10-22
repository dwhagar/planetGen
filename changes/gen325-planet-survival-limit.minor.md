**Headline:** A planet or belt now keeps only the orbits its place's passing systems have not stripped, so a nuclear disc holds planets to a few AU, a dense core to a fraction of one, and the Sun's neighbourhood to a few hundred.

### Changed
- **The planet survival limit (GEN.325).** `physics/planet_survival.py` takes any density and dispersion and answers the largest semi-major axis whose expected number of ionising passes over the star's age (at most 12 Gy, a 0.5 solar mass host) is one: `N(a) = n pi a^2 v t (1 + 2 G M / (a v^2))` with `v = sqrt(6) sigma`, solved in closed form. It reproduces the table of the core note (0.03 AU at 0.5 pc from Sgr A*, 0.4 AU at 4 pc, 2 AU at 10 pc, 3 to 11 AU in the nuclear disc) and the expected passes at the four rows of the cluster note (188, 11, 2.6 and 2.1).
- **`PlaceTerms`** carries the mean system mass of its place, so `systems_per_pc3` follows from the mass density; `place_terms.planet_limit_au(terms, age)` is the limit of a place, and `StarSystem._orbit_ceiling_au` takes the tighter of it and the stable orbit. Comets keep their own ceiling.
- The metallicity factor (GEN.236) is unchanged and still thins a metal-poor star's giants and planets; the cull is that factor and this limit together.
generation-output: changed
