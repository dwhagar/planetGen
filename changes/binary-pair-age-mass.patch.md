### Fixed

- Both stars of a binary now share one age (GEN.53). A `--star-type`
  primary's companion is born at the primary's age, and the age adjustment
  for planets now settles one age for the pair: old enough for the planets
  around either star, and no older than either star's generated state
  allows. Before, a wide pair's stars were aged separately (an M2V of 13.8
  Gy next to an M8V of 0.45 Gy) and a `--star-type` companion rolled its
  own age.
- A `--star-type` primary's companion is now a real star for its mass
  (GEN.54). It comes from the population model at a fraction of the
  primary's mass, so its type, temperature and luminosity follow from that
  mass. Before, it kept the requested type's temperature and luminosity
  with an unrelated mass (a "G2V" of 0.17 Msun).
