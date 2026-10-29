**Headline:** `planetgen benchmark` times the real runs instead of a sequence of its own.

### Changed
- **The benchmark is `plan` and then `galaxy`.** The scatters, the backfill, the population and the settle are stages of the galaxy run (PERF.172), so the benchmark no longer runs the scatters before the fill as commands of its own.
- **The stage timings come from the runs' records.** The report reads each stage's wall time, units and worker seconds from the run records in Redis (the figures the Queue page shows) instead of the stage rows in the control database.
