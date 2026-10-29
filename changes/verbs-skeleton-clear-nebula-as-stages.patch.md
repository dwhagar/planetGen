**Headline:** The skeleton, clear-scatter and clean-nebula-field commands now run as stages of the shared run, so they can share one run with other commands.

### Changed
- **Three one-step commands start their run through the shared stage list (PERF.202).** `planetgen plan` (the skeleton), `planetgen clear-scatter` and `planetgen clean-nebula-field` each used to start a one-step run of their own; they now follow a stage of the shared run like the other commands, so two of them can run in one run and keep their refusals, in order. The output is unchanged; a stopped run of these now also records the galaxy seed and build so `planetgen finish-run` can carry it on.
