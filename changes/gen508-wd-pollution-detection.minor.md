**Headline:** A tuned white-dwarf redox width is now compared with what the polluted white dwarfs would show, not with the bodies, and the model says the shipped width is too wide.

### Added
- **The white-dwarf pollution detection model (GEN.508).** `planetgen/physics/wd_pollution.py` holds the 24 polluted white dwarfs with a measured oxygen fugacity (Doyle 2020 and 2023), the papers' rule for what is reported (both bounds, an upper bound only, or an upper limit only), an absolute abundance error fitted to them, and `expected_sample`, the sample a population would show. It is registered in the tuning harness as the selection model `wd-pollution-detection`; the harness imports it when a target names it. The spread of the 22 reported medians is explained by the measurement, so the intrinsic width of the redox rows a3 and a4 is not identified (0 to about 0.15) and their shipped 0.45 is too wide; the retune is left to GEN.513.

### Changed
- The tuning harness imports the module that builds a selection model on first use (`tuning_harness.selection_model`). No draw changes.
