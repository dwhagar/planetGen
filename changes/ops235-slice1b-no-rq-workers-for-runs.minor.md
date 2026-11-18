**Headline:** Nothing starts RQ workers for a run or the maintenance lane any more: the supervisor's launcher forks every unit, and the code that asked for workers is gone.

### Changed
- **OPS.235 slice 1b.** Deleted the supervisor's `WorkerStarts` duty and `nudge` (the asks for workers), `workers.ensure_workers`, `start_workers`, `registered`, `maintenance_cap` and the starting counters, and `wrapper.stage_ensure` (it was already a no-op for any run with an id). The wrapper's follower calls a plain optional `wrapper.pump` while it waits for a stage; production leaves it unset and the test harness (`tests/inline_queue.py`) sets it to hand out the ready rows in the test's process. `held.retry_held` wakes the supervisor with `ensure_started`. The Generate page's own worker (`cli/worker.py`, `web/jobs.py`) stays until OPS.244 moves it to the ledger; slice 2 removes the `rq` dependency then.
- Reuse: the launcher's hand-out and `supervisor.ensure_started` do what the asks and starts did; `inline_queue.run_ready_rows` is the pump.

### Removed
- `WorkerStarts`, `supervisor.nudge` and `nudges_key`, `workers.ensure_workers`, `start_workers`, `registered`, `lane_of`, `maintenance_cap`, `wrapper.stage_ensure` and their tests (`test_supervisor` ask tests, `test_workers` start tests, `test_maintenance_lane` and `test_queue_lanes` start tests, the `stage_ensure` cases of `test_stage_worker_cap`).
