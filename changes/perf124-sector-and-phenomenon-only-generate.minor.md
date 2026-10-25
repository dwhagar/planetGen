**Headline:** `planetgen sector` and `planetgen phenomenon` now only generate; the population pass and the sector paths are commands of their own.

### Changed
- **`sector` no longer takes `--population` (PERF.124).** A `sector` run saves its sectors and stops; the population pass is `planetgen population`, or a stage of a `galaxy` run.
- **`phenomenon` no longer saves sector paths, and `--no-settle` is gone (PERF.124).** `planetgen phenomenon --sector-id N` saves the phenomenon and nothing else.

### Added
- **`planetgen settle --sector-id ID...`.** Saves the paths of the named sectors and the sectors around them (the step that follows `phenomenon --sector-id` for a black hole or neutron star); it takes the place of `--since` and the two cannot be combined.
