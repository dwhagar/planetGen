### Added
- `planetgen galaxy` and `planetgen plan` draw one "Whole job (stage N of M)" bar on the command line, with the time elapsed and the time left across all their stages, from the stored stage times (PERF.55). A finished run also stores its whole-job layers per second, counting every layer the scatter stages visited, even those that placed nothing.
