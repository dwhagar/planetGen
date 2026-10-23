**Headline:** A drawn planet is now tested against the chance it would have been measured to mass 25 percent and radius 8 percent, the cut of the published density samples.

### Added
- **The mass-radius sample selection model (GEN.505).** `planetgen/physics/mass_radius_sample.py` holds an error model of the radial-velocity amplitude and the transit depth fitted to the NASA Exoplanet Archive's transiting planets with a measured mass (`scripts/fit_mass_radius_sample.py`), the chance a planet passes the cut, and one planet measured with its errors. It is registered in the tuning harness as `mass-radius-sample`. The cut is a weak filter on mass, radius and period because observers follow a planet up until it is precise; it conditions on a published mass and radius, and detection is GEN.504.

### Changed
- The tuning harness imports the module that builds the mass-radius selection model on first use. No draw changes.
