**Headline:** A run or job whose process was killed no longer blocks every new run or holds up an update, and `update.sh --reset-jobs` resets every job and the queue when something else is stuck.

### Fixed
- A run that owns the machine is the Redis key `planetgen:run:active`, which had no expiry and was freed only by the run's own process. A run killed with `kill -9`, by the out-of-memory killer or by a restart left it set, and every new run was refused as busy for ever. A run now records which process follows it and that process says it is alive at every poll (`runrecord.touch`); when a new run finds the machine held, it ends a run whose process is gone (same host: the process id is not running; another host, or a run recorded before this: no sign of life for two minutes) as failed, frees the machine and starts. The waiting units of the dead run do nothing.
- A queued job (a Generate page job, an API job) is queued with no time limit, so RQ never gave up the entry of one whose worker was killed: it stayed "running" in Redis, `update.sh` listed it, could not stop it, waited ten minutes and refused to update. The update guard and every worker now clear such jobs when they look (`workers.reap_dead_jobs`): a job whose worker is gone fails with a message saying so (or is queued again when it has a retry left), and a registry entry whose job has expired is dropped.

### Added
- `python3 -m planetgen.cli.reset_jobs` (`--yes`, `--dry-run`) resets every job and the queue: it stops what the update guard lists, kills the queue workers on the machine, removes the RQ keys and the run records in Redis (the admin's worker count and the caches stay) and the rows of the jobs table, frees and un-pauses the work lease, and prints what it found and removed.
- `sudo ./update.sh --reset-jobs` runs it first, without asking, before the guard is taken.

TODO items resolved: PERF.142, PERF.143 / Added: none
