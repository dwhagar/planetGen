**Headline:** A setting runs the units of a generation run as forked jobs handed out by the supervisor, off by default.

### Added
- `queue.forked_jobs` (off): when on, the cut of a stage of a generation run makes its unit rows `ready` and queues nothing in RQ; the supervisor's new `Launches` duty hands them out with `dispatch.dispatch` and runs each as a job forked from the lane's warm template (OPS.228, slice B). Runs that carry on from an earlier one, and runs with no ledger, stay with the workers.
- The `open` row of a stage names the Redis run id (`run_id`) beside the stage's index, so the launcher can build the unit's job from the ledger alone; `planrows.release_units` makes a stage's unit rows ready.

### Changed
- The launcher leaves a unit's attempt to the unit (`units.run_unit_on`, as in a worker), marks a row `running` when it forks, and settles a process that died with its row still `running`: handed out once more, then failed.
