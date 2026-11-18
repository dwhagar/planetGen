**Headline:** Every run is a run of the job ledger: the units that went to RQ workers are gone, and the update migrates the database before the math check.
### Changed
- A run's units are rows of the control database's ledger that the supervisor's launcher hands out in forked jobs. The RQ unit path (`units.enqueue_units`, `run_unit`, `open_stage`, the unit hooks, the warm worker's unit preparation) is deleted, with no fallback for a run with no ledger (OPS.271).
- A failed or waiting unit row is made ready again by `held.retry_held`; the opening of a stage that does not exist is the end of the run.
- `update.sh` migrates the database before it runs the math check, because the math check is a run and has no ledger before the migration.
