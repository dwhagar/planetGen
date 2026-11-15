**Headline:** An update no longer cancels a generation run: it holds the run "stopped for update", and the end of a good update carries it on.

### Changed

- OPS.222 part A (job engine stage 3, L8): when the operator agrees to stop the work, `update.sh` (`updateguard`) drains a generation run in the run record instead of canceling it: `held.drain_run` holds it with the reason "stopped for update" (the machine is free, no worker takes a unit of it), each running unit stops at its next piece (the engine finds the run no longer live) and goes back to `waiting` in the ledger, not failed. The run is drained before the process that follows it is stopped, so that process cannot end it. At the end of a good update `held.resume_after_update` carries each drained run on (`retry_held`): every unit the ledger does not say is done is queued again unless the queue still has it, the opening of a stage the hold skipped is queued, and the output says what was queued. When the generator epoch changed with the update (`generator_epoch`) the run stays held and the output says to run it again. A run started by systemd is still never stopped (the update refuses while it generates), and a Generate job stays protected as before.
- `retry_held` queues every unfinished unit of the run, not only the held ones, so it also carries a drained run. Reused: `retry_held`, `hold_run`, `units.requeue_unit`, the engine's stop between two pieces, `generator_epoch`, the guard's source and cancel path.

Part B (finish-run plans from the rows, `REWIND_WITH` deleted) waits for the checkpoints of OPS.219.

generation-output: unchanged
generation-affects: none
