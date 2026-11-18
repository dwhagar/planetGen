**Headline:** The `rq` library is gone from the install: nothing queues a job on Redis for a worker any more, so the update installs one package less.

### Removed
- **OPS.235 slice 2.** The `rq` dependency (`setup.py`, `scripts/install-python-deps.sh`, `requirements.lock`, `requirements-server.lock`, with `croniter`, `python-dateutil` and `six`, which only it brought) and every module that used it: `planetgen.cli.worker` (the queue worker process), `redisqueue.TimedQueue`, `queue`, `requeue_job`, `fetch_job`, `worker_argv` and `serializer_for`, `allowlist.GuardedJob`, `workers.live_workers`, `cancel_now` and `reap_dead_jobs`, and `claims.claim_for`, `default_claim` and `lane_of`.
- The shared queue list of the update guard: `QueuedJobSource` now lists only the site's light jobs (a run's units, a page's steps and the light jobs are rows of the ledger or records the supervisor forks); stop-all no longer empties RQ queues or clears dead RQ workers.
- `docs/reference/cli.md`'s `worker` command.

### Changed
- `util/memwatch` counts the supervisor's template processes (`planetgen.queue.template` and the jobs forked from it) as the "workers" of the memory peak line, and `cli/benchmark_suite` ends those between steps, where it ended queue workers.
- `reset-jobs` still removes the `rq:*` keys an older release left in Redis, read as plain keys, and ends the worker processes their records name; it no longer needs the `rq` library for either.
