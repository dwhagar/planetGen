**Headline:** The galactic centre holds the observed nuclear star cluster and nuclear stellar disc, in place of ADM.49's Gaussian core.

### Added
- **The nuclear star cluster and the nuclear stellar disc (GEN.280).** Two more terms of the density model, scaled to the galaxy's bulge: a Hernquist cluster of 2.5e7 solar masses (half-light radius 3.8 pc) with its cusp held flat inside the radius where its density equals the core ball's mean, so the cell at the centre holds the cluster's mass inside 2 pc (about 1.5e7 systems), and an exponential disc 90 pc long and 35 pc high of 1.5e9 solar masses. They are the `nuclear_cluster` and `nuclear_disc` populations of a sector's profile (and `bright_stars.population`), count in the outline's bound, the bright-star scatter's bound and the Galaxy Map's prisms, and give a profile the region `nuclear`. A cell this dense is the case the cell star ceiling (GEN.281) is for.
- **`planetgen plan --nuclear-amplitude`** (default 1; the Generate page's "Nuclear structures") scales both terms, 0 being none.
- **Five math checks** compare the default galaxy's cluster mass, half-light radius, density at 1 parsec, the disc's mass inside 200 parsecs and its height with the published figures.

### Changed
- **`galaxy_shape.core_amplitude` and `--core-density` are gone** (migration 0099, schema v99): the column is renamed `nuclear_amplitude`, set to 0 for a galaxy planned before it, so its density is unchanged. `bright_stars.population` becomes `VARCHAR(16)`.

### Fixed
- **Filling a dense cell no longer takes quadratic time.** Each pre-placed star was checked for overlap against every earlier system of its cell one at a time in Python, which is hours for a cell of 40,000; it is one NumPy pass now (the sector's existing neighbour arrays).

TODO items resolved: GEN.280. Added: none.
