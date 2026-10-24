**Headline:** A drawn comet is now tested against the chance it would have a published coma composition and the heliocentric distance it would be measured at.

### Added
- **The comet admission selection model (GEN.507).** `planetgen/physics/comet_admission.py` holds the share of known comets with a measured composition by group and perihelion distance, and a distance model (parent molecules and radicals separately, truncated at the perihelion), fitted to the Robinson et al. 2024 database joined to the JPL Small-Body Database orbits (`scripts/fit_comet_admission.py`). It is registered in the tuning harness as `comet-admission`. It conditions on a known comet and does not use nucleus size or activity.

### Changed
- The tuning harness imports the module that builds the comet selection model on first use. No draw changes.
