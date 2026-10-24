**Headline:** A tuned atmosphere share is now compared with the planets JWST has observed, with a model of which planets those are.

### Added
- **The JWST target-selection model (GEN.506).** `planetgen/physics/jwst_selection.py` holds a logistic regression of "has a published JWST spectrum" on the transmission and emission spectroscopy metrics of Kempton et al. 2018 and a radius-class offset, fitted to the NASA Exoplanet Archive (`scripts/fit_jwst_selection.py`), host magnitudes for a generated planet, and the chance a planet is on the list. It is registered in the tuning harness as `jwst-target-selection`. It reproduces the number observed in each radius class; it conditions on a known transiting planet with a mass.

### Changed
- The tuning harness imports the module that builds the JWST selection model on first use. No draw changes.
