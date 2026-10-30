**Headline:** The planet architecture is built as a tested model, with the observed planet occurrence the generator will be held to.

### Added
- **The planet architecture** (GEN.728, GEN.729, GEN.730, GEN.828, GEN.829, GEN.920): a Kepler-like inner system drawn at the observed share, oligarch chains thinned by Pu and Wu's instability law, gas giants with metallicity, and hot-Jupiter and inward migration (`physics/planet_architecture.py`). Not wired into generation yet; it rides a later merge train.
- **The observed bands** of planet occurrence with their sources (`physics/data/planet_occurrence.json`, `physics/planet_occurrence.py`) and a Kepler multiplicity harness (`physics/kepler_multiplicity.py`).
- **TEST.395**: the occurrence census now checks Kepler-like share, multiplicity, giants, hot Jupiters, habitable-zone Earths and spacing against those bands for M, K, G and F stars.
- Registry rows PA1 to PA12 for the new stand-ins (GEN.725).

### Changed
- The planet mass law's solids and gas steps are separate functions (`solids_surface_density_gcm2`, `gas_outcome`) with the same results.
