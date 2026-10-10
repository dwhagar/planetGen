### Added
- Generation rates are recorded per kind of work and worker count (PERF.32): sector fills, bright-star layers and now phenomenon layers. The Stats page shows the worker count of each rate and has a Reset stats button; the time estimate reads the rate for the run's own worker count, blending the neighbouring counts when it has none.

### Changed
- Recorded rates are deleted by the first run of a new version, since they describe the release, Python and machine that measured them. Control schema v12 replaces the old `generation_stats` table.
