**Headline:** A generation run's record is kept in the control database, so a restart of the machine or of Redis no longer loses the run (OPS.284).

### Fixed

- OPS.284 (bug): a run's record (its steps, the galaxy seed and build it began on, each stage's clocks and conclusion) lived only in Redis, so a restart of the machine or of Redis left `planetgen finish-run` with "nothing to finish" and the run's ledger rows `running` for ever (found by the OPS.247 trial, which lost two runs of 11 and 3 hours to container restarts). The record is now written to the control database (`engine_run_records`) at each change that matters, and a new supervisor duty (`LostRuns`) and `finish-run` put back the runs Redis lost and end the ones that were live as runs whose process is gone (failed, killed, ledger rows ended); a run the Queue page deleted is not brought back. Redis stays the cache the pages read; there is no fallback.

generation-output: unchanged
generation-affects: none
