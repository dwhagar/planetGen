**Headline:** The galaxy skeleton is built by a command on the work queue, like the scatters.

### Changed
- **`plan-skeleton` is a command (PERF.106).** `planetgen plan` queues one unit that builds the skeleton (the shape, the clusters and satellites, the outline and the naming key) through the wrapper, as a run of its own, and logs the lines the unit returns; a unit runs alone from its JSON arguments (`planetgen command-run`). A reach or seed the build refuses still ends the plan with the same message and exit code 1. The stage's row in the stage history is stored by the unit that completes it.
- **`run_plan.build_skeleton` is gone** with no replacement; the command takes the plan's shape options as plain parameters.
- **`docs/design/redis-and-the-work-queue.md` describes it.**

TODO items resolved: none (slice four of PERF.106). Added: none.
