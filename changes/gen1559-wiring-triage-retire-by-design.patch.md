**Headline:** The wiring register marks seven reference helpers by design and retires two duplicate helpers.

### Changed
- **Wiring triage, retire and by-design rows (GEN.1559).** `photometry.bolometric_correction_v` and `blackbody_colour` (the reference fits the sky and star colour are held to) and five `clathrate_brine` functions (the reference curves the tests hold the wired ones to) are marked by design in the wiring register, each with its reason. `comet_species.low_over_h2o` and `high_over_h2o`, one-line copies of `bracket_over_h2o`, are deleted; the register now names `bracket_over_h2o` for GEN.1021. Left for their owners: `planet_zones.flux_at` (asked of Lane 11) and `clearance_ledger.log_uniform_quantile` (belt birth, GEN.1068).
generation-output: unchanged
