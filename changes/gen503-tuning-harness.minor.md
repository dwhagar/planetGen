**Headline:** A tunable composition draw is now tested against the published statistic it is tuned to, and a tunable row with no target is refused.

### Added
- **The shared tuning harness (GEN.503).** `planetgen/generation/tuning_harness.py` and the target table `tuning_targets.json` hold the published statistics (value, sigma, source, label, validity range, selection model) that tunable draws are tuned to; `simulate_and_compare` draws a row at many seeds and passes it within 2 sigma. The loader refuses a tunable row with no target, a Solar System value used as a target, and a snippet or recalled value used to tune. The white-dwarf redox rows a3 and a4 are the worked example.

### Changed
- A possibility-table row can be marked `tunable`; it changes no draw.
