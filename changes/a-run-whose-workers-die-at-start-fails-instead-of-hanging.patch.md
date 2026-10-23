**Headline:** A run whose workers keep dying at start now fails with the exit code instead of starting workers for ever.

### Fixed
- **The queue executor started a new worker every few milliseconds when its workers could not start (TEST.295).** `RQExecutor` starts workers for as long as tasks wait in the queue. A worker that exits at once (an import error, a process killed for lack of memory, a Redis it cannot reach) is replaced at once, so a run on a machine under pressure could spin for ever, with its workers' output discarded: 752 workers in 21 s in a test, and `test_one_worker_draws_what_two_and_four_workers_draw` reached the 900 s timeout in a loaded gate while it waited with tasks in the queue and workers being started. The executor now counts the workers that exit with a failure in a row (reset when any task finishes) and, at `WORKER_CRASH_LIMIT` (10), fails the run with a `WorkerDied` naming the exit code and the command that shows why. A new test with workers that exit at once fails without the limit.

TODO items resolved: TEST.295. Added: none.
