### Fixed
- **The bright-star backfill after a `galaxy` run has its own progress
  bar with a working ETA** (PERF.28). A `--slot` run (the Sector Map's
  "Generate neighborhood" among them) used to backfill inside the
  sector's bar, which sat at 0 of 1 with no ETA for minutes; the backfill
  now runs at the end of the run, counting the sectors it visits.
