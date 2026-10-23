**Headline:** A drawn planet is now tested against the chance the Kepler pipeline would have found it: the orbit's geometry, its signal against the host's noise, and a recovery curve fitted to the DR25 injections.

### Added
- **The Kepler transit-plus-recovery selection model (GEN.504).** `planetgen/physics/kepler_recovery.py` holds the geometry and duration of a circular transit, the multiple event statistic against a host's noise (a noise level from the host's Kepler magnitude and radius, fitted to the DR25 stars), and a recovery curve fitted by maximum likelihood to the 146,257 DR25 on-target injections (`scripts/fit_kepler_recovery.py`). It is registered in the tuning harness as `kepler-recovery`. It is the detection step only (the vetting removes more) and holds for periods of 0.5 to 500 days.

### Changed
- The tuning harness imports the module that builds the Kepler selection model on first use. No draw changes.
