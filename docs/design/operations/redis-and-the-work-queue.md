# Redis and the work queue: how planetGen uses it, how it should be used, and what to change

Boss (2026-10-14 20:01Z) asked for three things in order: look at how planetGen uses Redis now, research how Redis is best used and what the industry standards are, and recommend how to change the task process to use Redis fully and correctly. This note is that, in one piece. It is the research PERF.86 (redo the work queue, GitHub issue [#1229](https://github.com/dwhagar/planetGen/issues/1229)) asked for before the rebuild starts, so section 4 gives Boss options first and section 5 recommends one.

Informs: PERF.86 (and the items it folds in or drops: PERF.29, PERF.30, PERF.40, PERF.78, ADM.15), PERF.20, PERF.33, ADM.22, UX.3, DOC.20, DOC.44, OPS.21, OPS.46.

Status (2026-10-29): research of 2026-10-14, then built as PERF.86 (PERF.95 to PERF.100 and the command steps PERF.106 to PERF.108, PERF.118, PERF.119, PERF.121 to PERF.128; the MySQL job tree and its tables are gone, PERF.98 and PERF.100). The findings and the options are as written then; the "As built" paragraphs in section 5.12 say what was built. Decisions marked "Boss" are his; everything else is a recommendation.

It builds on [work-queue-audit.md](work-queue-audit.md), [performance-eta-queue-and-caching.md](performance-eta-queue-and-caching.md) and [sampling-backfill-and-resume.md](../generation/sampling-backfill-and-resume.md) and does not repeat them. Where this note and those disagree, this note is newer and was checked against `main` at 8.0.1089.

Evidence tags: [C] computed or measured for this note (scripts in `/mnt/project-files/research/redis-review/`); [S] read in the repository, in a downloaded package, or in documentation or a search result (URLs under Sources); [R] recalled and unconfirmed. Measurements come from one shared sandbox (Redis 7.0.15 from the distribution on loopback, Python 3.13, `rq` 2.12.0, `redis` 8.1.0, no persistence unless stated). They show orders of magnitude and ratios, not the production server's speed.

## Findings in short

1. **Redis does not hold the plan.** A run's tasks wait in the supervising process's memory; only `2 * workers` of them are ever in Redis at once ([S] `queue/work.py` `submit`, `_dispatch`). So Redis gives no durability, no visibility of what is left, and no way for anyone but the supervisor to cancel or resume. That, not Redis, is where the "task and subtask" complexity comes from.
2. **There are three stores of job state, not one.** RQ in Redis, the MySQL job tree (`work_jobs`, `work_tasks`, `work_lease`, `work_job_args`) and files in the jobs directory (`active`, `job.json`, `state.json`, `progress.json`, `cancel`, `output.log`). Each has its own cancel path and its own "one at a time" lock (a file and a MySQL row).
3. **Redis is used well in four small places** (rate limits, login lockouts, the count cache, tile single-flight) and as a process launcher in the big one (the queue). Those four are fine; the big one is what to redo.
4. **The server's own configuration is unmanaged.** The installer apt-installs `redis-server` and starts it; nothing sets persistence, memory limit or a password. With the distribution defaults a crash loses everything written in the last minutes: 0 of 100 keys survived `kill -9` shortly after writing, against 100 of 100 with `appendonly yes` [C]. The login lockouts live only in Redis, so a crash also forgets who was locked out.
5. **Two real defects in the current use** (section 1.4): `api_jobs.ensure_workers` counts dead workers as alive forever after a worker is killed (measured), and workers that serve more than one queue lose RQ's reliable-dequeue path (source).
6. **RQ 2.12.0 already has what a plain-Redis design needs**: enqueue many jobs at once, stop and cancel, callbacks, results history, a killed-work-horse handler, `max_idle_time`, `unique` job ids. Section 2.2 lists them with versions. A hand-rolled queue on lists or streams would rebuild them.
7. **Two traps in "run the stages as queue items", both measured** [C]. Burst workers exit when the queue is momentarily empty, which a stage barrier makes happen: 24.0 s against 10.6 s with idle-timeout workers (9.6 s ideal) in a 3-stage test, so workers must outlive barriers (`max_idle_time`, not `--burst`). And RQ's `depends_on` is quadratic when one job waits on many: a gate that depends on N units took 2.3 s for 200, 10.1 s for 500 and 45.4 s for 1,000 units, against a flat 7 ms per unit when the last unit of a stage advances the run itself with an atomic Redis set.
8. **Recommendation:** keep RQ and Redis, make Redis the only place a run lives (a run record that lists its stages, each stage's units as RQ jobs, the last unit of a stage queuing the next, progress as counters plus an event stream), drop the MySQL job tree, the jobs-directory files and the supervisor's queue, and configure the Redis server properly. Section 5.

## 1. How planetGen uses Redis today

### 1.1 Every use

Update 9:25 AM EDT 2026-10-24 (OPS.66, as the lane reported): `WorkQueue`, `RQExecutor`, `RQFuture` and `Channel` are gone. A run's units are RQ jobs queued by `queue/units.py` and `queue/wrapper.py` on the shared queues; the rows and steps below that name the removed classes describe the design before that change.

| Use | Code | Redis data | Lifetime | If Redis is down |
|---|---|---|---|---|
| Per-run task queue | `queue/redisqueue.py` `RQExecutor`, `queue/work.py` | RQ queue `planetgen-run-<16 hex>`, one job per task, pickled arguments | deleted when the run ends | `WorkQueue` logs it and runs every task in the supervising process, one at a time |
| Run progress back to the supervisor | `redisqueue.Channel` (scatter layers, sector fills) | list `<queue>:<name>`, pickled items, `RPUSH`/`BLPOP` | deleted at the end, 1 day expiry as a backstop | the run fails with `OSError` |
| Generate page jobs | `web/jobs.py`, `web/job_runner.py`, `cli/worker.py` | one RQ queue per job, `planetgen-web-<job id>`, one job, one detached burst worker | `result_ttl=0`, `failure_ttl` 1 day | `start_job` raises `OSError`: no job starts |
| API jobs (neighborhood, regenerate, wiki, estimates, one-off system) | `queue/api_jobs.py` | shared queue `planetgen-web`, jobs `planetgen-api-<16 hex>`, counter `planetgen-api-workers-starting` (15 s) | `result_ttl` and `failure_ttl` 1 day | `NoQueue`; `command_and_wait` runs the command in the request instead |
| Request rate limits (Flask-Limiter) | `api/config.py` | `planetgen:limiter:*` (the `limits` library's keys) | per window | falls back to counting in memory |
| Login lockouts (SEC.30) | `admin/throttle.py` `RedisStore`, `api/loginguard.py`, `cli/lockouts.py`, `cli/lift_lockout.py` | one JSON string per address or username, `planetgen:limiter:login:<scope>:<subject>`, written under `WATCH`, 7 day expiry | until idle | counts in memory, warns once |
| Stored counts (PERF.64) | `db/countcache.py` | `planetgen:counts:<db>:<digest>` (6 h) and `...:lock` (10 min), `SET NX EX` | expiry | in-process cache only |
| Tile single-flight (PERF.38) | `web/lib/singleflight.py` | `planetgen:build:<name>`, `SET NX EX 60` | expiry | every caller builds |

All of it is in database 0 of one server at `redis.url` (default `redis://127.0.0.1:6379/0`, no password) [S]. The tests need it (`PLANETGEN_TEST_REDIS_URL`) and CI runs `redis:7` [S].

Not built, although earlier notes describe them: the Redis progress record (`planetgen:progress:<run>`), the epoch counters and the in-flight dedupe keys of [performance-eta-queue-and-caching.md](performance-eta-queue-and-caching.md) sections 2.5 and 4. `GENERATION_QUEUE` in `redisqueue.py` is declared and listed as a default queue but nothing is ever queued on it [S, grep]. `cli/worker.py`'s docstring names `redisqueue.ensure_workers`, which does not exist.

### 1.2 How a generation run flows

1. **Admin presses Generate.** `web/jobs.py` `start_job` makes a directory under the jobs directory, takes the job lock (the hash `planetgen:job:active` in Redis, PERF.154), writes `job.json`, queues `web.job_runner.run(<dir>)` on a queue of its own and starts one detached burst worker for it.
2. **The job runner** (inside that RQ job) runs each step as a subprocess, `planetgen plan` or `planetgen galaxy ...`, writes `state.json`, appends to `output.log`, watches for a `cancel` file, and hangs a "web-job" node and one "step" node in the MySQL job tree.
3. **The step's `planetgen` process** is the supervisor. Its `WorkQueue` opens a MySQL node, takes the `work_lease` row (one run's pool on the machine at a time; heartbeat every 5 s, stale after 30 s), starts a per-run RQ queue and up to `worker_count()` burst workers, and then feeds tasks to them: at most `2 * workers` are in Redis at once, the rest wait in a Python list ([S] `submit`, lines 1498 to 1535).
4. **Each finished task** comes back to the supervisor by polling RQ every 50 ms (`fetch_many`), is recorded in `work_tasks`, and runs its `on_done` callback in the supervisor, which updates the progress bars and ETA. The scatter and sector fills also push per-layer progress through the Redis `Channel`.
5. **The Queue page** reads the MySQL tree; the Generate page reads `state.json` and `progress.json`; the API's `/api/jobs/<id>` reads RQ.

So Redis carries tasks one hop to a worker and results one hop back. Everything else about the run (what is left, who is paused, what to show) is the supervisor's memory, MySQL rows and files.

### 1.3 Where each piece of state lives

| State | Redis / RQ | MySQL | Files |
|---|---|---|---|
| What is left to do | at most `2 * workers` tasks | task counts (`set_total`) | |
| What is running, who runs it | RQ started registry, worker registry | `work_tasks.state`, `work_lease.holder` | `runner.pid`, `active` |
| Pause, cancel | RQ `cancel()` on shutdown only | `work_jobs` control rows, `work_lease.paused` | `cancel` file |
| Progress, ETA | `Channel` lists (live) | node totals, ETA roll-up | `progress.json` |
| Timings, stats | | `work_tasks`, `generation_stats`, `generation_stage_runs` | |
| Log | | | `output.log` |
| Final result | RQ result (1 day) | node state | `state.json` |

### 1.4 Findings

**F1. The queue cannot survive its supervisor.** The plan is in memory. A killed `planetgen` process leaves RQ jobs for the workers to finish and nobody to collect them, a stale lease for 30 s, and a half-run with no record of what remained. PERF.29 and PERF.30 (record what a sector still needs, resume) exist to rebuild by other means what a queue that held the whole plan would have. Boss has since decided there is no automatic resume (PERF.86), which removes the need but not the problem: after a crash the admin cannot see what was left.

**F2. Three stores, three cancels, two locks.** Cancel is a file for web jobs, a MySQL control row for work queues, and `job.cancel()` at shutdown for RQ. One-at-a-time is a file (`active`) for web jobs and a MySQL row (`work_lease`) for pools. Each store needs its own staleness rule (runner pid and a 300 s grace, 30 s heartbeat, 15 s starting grace).

**F3. Bookkeeping costs more than the work it records.** The `work.py` header counts about 60 database statements per task [S, comment at `CHUNKS_PER_WORKER`]; that is why `submit_each` sizes chunks by time (PERF.79, PERF.83). The Redis side of a task is cheap: an RQ job costs about 1.8 KB queued and about 7.2 KB with its result kept an hour, so all 2,041 layers of a default-scale scatter are 3.7 MB queued and 15 MB finished if results are kept an hour [C]; the memory follows the number of jobs, so a design must keep the number of jobs bounded (section 5.1). Queuing 2,041 jobs takes 0.90 s one at a time and 0.36 s with `enqueue_many` [C].

**F4. Dead workers are counted as alive (defect).** `api_jobs.ensure_workers` starts `min(waiting, worker_count()) - Worker.count(queue) - starting` workers. `Worker.count` is the size of a Redis set that only a worker's start-up maintenance cleans [S, `rq/worker_registration.py`, `rq/worker/base.py`]. Measured: after `kill -9` of a worker the count stayed at 1 for the whole test, and still stayed at 1 after the worker's own key had expired (worker TTL 20 s, key gone at 80 s) [C]. When the server's OOM killer takes an API worker (it has taken Apache before, see [server-checklist.md](../../deployment/server-checklist.md)), queued API jobs wait on a phantom, and with a cap of one nothing replaces it. Not reproduced end to end through planetGen's own code; the cause and the measurement are enough to file it.

**F5. A killed job is noticed late.** A job running under a worker that was `kill -9`'d stayed `started` for about 90 s (heartbeat TTL = monitoring interval 30 s + 60 s), then went back to the queue and ran again because of `Retry(max=1)`; that cleanup only happened because the test polled the registry [C]. `RQExecutor` covers this inside a run by watching its own child processes. A design that loses the supervisor needs the same covering elsewhere (a worker-killed handler, section 5.6).

**F6. Multi-queue workers skip RQ's reliable dequeue.** `Queue.dequeue_any` uses `BLMOVE` into an intermediate queue only when a worker serves exactly one queue; with several it uses `BLPOP`, and a worker that dies between the pop and "started" loses the job id [S, `rq/queue.py`, version 2.12.0]. Today's task and web-job workers serve one queue each (good). `cli/worker.py` with no arguments, and the `[interactive, bulk]` worker of PERF.40, would serve two.

**F7. The server is configured by default.** `update.sh` installs and starts `redis-server` and only checks that `PING` answers (`scripts/deploy-common.sh` `ensure_redis`). The distribution's `redis.conf` has `appendonly no` and, with the save lines commented out, Redis's built-in save points `3600 1 300 100 60 10000` [C, `CONFIG GET save`]. There is no memory limit (so no eviction, which is right for a queue, but also no guard), no password (it binds loopback and is in protected mode, which is adequate for a single-user host), and nothing documents any of this: DOC.44 notes that no install guide lists Redis.

**F8. Pickle on an unauthenticated store.** RQ jobs and the `Channel` use `pickle`. RQ's own documentation warns that its pickle default is not secure and that RQ should only run against a Redis you trust [S]. On loopback with one user it is acceptable; on a shared host any local user can queue code that a worker runs.

**F9. Client settings are ad hoc.** `redisqueue.connect()` builds a new client (a new pool) on every call and sets no timeouts, so a hung Redis hangs the calling web thread; `countcache` and `singleflight` each build their own client with 0.5 s timeouts; `api_jobs.wait` polls `status` every 50 ms with a fresh fetch and, for a queued job, a worker check ([S]). RQ itself raises a worker connection's `socket_timeout` above its blocking time when needed [S, `rq/worker/base.py` lines 430 to 440], so timeouts belong on the web-side clients, not the workers'.

**F10. Locks are efficiency locks, which is fine.** The single-flight and count claims are `SET NX EX` on one server, and a lost claim only means a duplicate build. That is the use Kleppmann's analysis allows a single-instance Redis lock for; the run-exclusivity lease is the only one that approaches correctness, and it guards CPU share, not data (section 2.6).

## 2. How Redis and RQ are best used

### 2.1 Redis as a queue

- **Lists** are the basic queue: `LPUSH`/`BRPOP`. The Redis documentation says lists suit "queue management for background worker systems" and points to `LMOVE`/`BLMOVE` for "safer queues": the job moves atomically to a processing list, so a crashed consumer does not lose it [S].
- **Streams with consumer groups** add explicit acknowledgement, a list of pending messages, claiming of messages from dead consumers (`XAUTOCLAIM`), and history [S]. They are the right primitive for an event log that several readers follow, and an alternative for the queue itself.
- **RQ** is a job library over lists: a queue is a list of job ids, each job a hash, with registries (sorted sets) for started, finished, failed, deferred, scheduled and canceled jobs, worker heartbeats, and an intermediate queue for reliable dequeue [S]. It is at-least-once: a job can run twice after a worker death, so jobs must be safe to repeat.
- **Delivery guarantees are a property of persistence too.** A queue in an unpersisted Redis loses queued work on a crash (2.3).

### 2.2 RQ features that matter here

Versions from the RQ documentation on GitHub (master) and checked against the installed 2.12.0 source where noted. PyPI's latest is 2.12.0 (2026-08-30), which is what `requirements.lock` pins [C, PyPI JSON].

| Need | RQ feature | Since | Notes |
|---|---|---|---|
| Queue many jobs at once | `Queue.enqueue_many` with `prepare_data` (one pipeline) | 1.9 | 2,041 jobs in 0.36 s [C]; `prepare_data` takes `timeout=`, not `job_timeout=` |
| Run B after A | `depends_on` (one job, a list, or `Dependency(allow_failure=...)`) | 1.x | dependents stay `deferred` until all succeed; a failed dependency leaves them deferred forever; **cost grows with the square of the dependency count** (200 deps: 2.3 s total, 500: 10.1 s, 1,000: 45.4 s [C]), fine for a chain, wrong for "after all N units" |
| Stop a running job now | `send_stop_job_command(connection, job_id)` | 1.7 | kills the work horse in under a second; status `stopped`; not retried |
| Remove a waiting job | `job.cancel()` | 1.10 | does not stop a running job (measured earlier, [performance-eta-queue-and-caching.md](performance-eta-queue-and-caching.md) 3.3) |
| Run code when a job ends | `on_success`, `on_failure`, `on_stopped`, `Callback(timeout=)` | 1.9, 1.14 | 60 s limit; do not run if the horse is hard-killed |
| Handle a killed horse | `work_horse_killed_handler=` | 1.13 | the hook for OOM kills |
| Retry | `Retry(max=, interval=)` | 1.5 | an `interval` needs a worker running the scheduler; interval-less retries are immediate and work on any worker |
| Do not queue twice | `unique=True` with an explicit `job_id` | 2.8 | raises `DuplicateJobError`; not allowed with `depends_on` |
| Results history | `job.latest_result(timeout=)`, `job.results()` | 1.12, 1.16 | the timeout form blocks instead of polling |
| Track a set of jobs | `Group` | 2.0 | removed when all its jobs expire |
| Workers that wait for work | `max_idle_time`, `--max-idle-time` | 1.14 | exits after N seconds with nothing to do |
| Several workers in one command | `rq worker-pool -n N` | 1.14 | also `--burst` |
| Completion pings | `Webhook` | 2.10 | not needed here |
| Cap how many jobs of one kind run at once | `RateLimit(key, concurrency)` | documented as 3.0 | **not released**: absent from the 2.12.0 source [C, grep]; do not plan on it |

### 2.3 Server configuration

- **Persistence.** The Redis documentation: RDB snapshots can lose "the latest minutes of data"; AOF with `appendfsync everysec` loses at most about a second; "use both persistence methods" for safety comparable to a database; `always` is "very slow" [S]. Measured on this note's sandbox: defaults lost all 100 keys on `kill -9`, `appendonly yes` + `everysec` kept all 100 [C]. For planetGen the queue and run records want AOF; the lockouts want it too.
- **Memory and eviction.** `maxmemory` unset means no limit; `noeviction` makes writes fail with an error at the limit instead of silently dropping keys; `allkeys-lru` and the other policies evict any key, including queue entries and lockouts. The documentation says a single instance used for both a cache and persistent keys should be split into two instances "if possible", and the `volatile-*` policies are the compromise [S]. Celery's Redis guide tells users to set `noeviction` or `allkeys-lru` to avoid its own broker keys vanishing [S].
- **Security.** "Designed to be accessed by trusted clients inside trusted environments"; bind to loopback or firewall the port; protected mode is the default safety net; ACLs (Redis 6+) are the recommended authentication, `requirepass` the legacy one; neither encrypts, so use TLS across a network [S]. Loopback only is a sound posture for planetGen's single-host install.
- **Monitoring.** `INFO` gives `used_memory`, `evicted_keys`, `rejected_connections`, `rdb_last_bgsave_status`, `aof_last_write_status`, keyspace hit and miss counts [S]. RQ adds `rq info` and a worker registry.

### 2.4 Client configuration

The redis-py production guide: set `socket_timeout` and `socket_connect_timeout`, set `health_check_interval` for connections that sit idle, retry `ConnectionError` and `TimeoutError` with backoff, and size one pool per process (`BlockingConnectionPool` waits instead of raising when the pool is full) [S, search result]. The numbers in community write-ups (5 s timeouts, 30 s health check, 50 connections) are reasonable starting points but are third-party advice.

### 2.5 Job design

Sidekiq's best-practices page, the best-known write-up of this for a Redis queue: pass small arguments (ids, not objects), because they are serialised into Redis; make jobs idempotent and transactional, because retries re-run half-finished jobs [S, search result]. planetGen already has the hard part: `task_seed(run_seed, key)` makes a task's output independent of worker count and order, and a sector save or backfill chunk is one transaction, so a repeated or killed task leaves nothing half-written ([sampling-backfill-and-resume.md](../generation/sampling-backfill-and-resume.md) 7.1). The change is to carry those arguments as plain JSON (numbers, strings, lists), which also removes pickle (F8) to the extent that matters.

### 2.6 Locks

Kleppmann's analysis of Redlock: a lock used for efficiency can live on a single Redis; a lock used for correctness needs fencing tokens or a store that enforces order, because a paused client can outlive its lease [S, search result]. Applied: single-flight and count claims are efficiency locks (a duplicate build is harmless); "one run at a time" is a resource guard, so a TTL lease on one Redis is enough; anything that must be exactly-once is guarded by a MySQL transaction or unique key, as the sector save is.

### 2.7 Libraries and the server itself

| Option | Fit |
|---|---|
| **RQ** (kept) | Minimal Python/Redis job library; fork per job; already pinned and used; has every feature above except concurrency limits. Slower than Celery per job (one benchmark: 20,000 small jobs, 10 workers, 51 s against Celery's 12 s [S, search result]), which does not matter at thousands of jobs per hour. |
| Celery | Full workflow engine (canvas, ETA, routing); heavy; on Redis it redelivers any task not acknowledged within the visibility timeout (default 1 hour), and tasks longer than that run again and again [S]. Wrong tool for hour-long stages. |
| Dramatiq, Huey, arq | Lighter alternatives with retries, middleware and (Dramatiq) rate limiting; none gives something RQ cannot, and moving costs a rewrite of the worker, pre-warm and test harness. No reason to switch. |
| Hand-written lists or streams | Exactly the registry, heartbeat and retry code RQ already has, to be owned and tested. |
| MySQL `SELECT ... FOR UPDATE SKIP LOCKED` queue | Supported in MySQL 8 for "queue-like tables" [S] (MariaDB from 10.6 [R]); durable and queryable, but each claim is a transaction on the database that is already the bottleneck, and Boss's direction is Redis. |

Redis 8 is tri-licensed with AGPLv3 added (2025-05) [S, search result]; Valkey is the BSD fork of Redis 7.2.4, and RQ supports Redis 5+ and Valkey 7.2+ [S]. planetGen uses only core commands, so either should work; this was not tested. The `redis-server` package version differs by distribution (7.0.15 here) [C]; `BLMOVE` needs 6.2 or later [S].

## 3. Measurements

All [C], scripts in `/mnt/project-files/research/redis-review/`.

| Question | Result |
|---|---|
| Queue 2,041 jobs (a default-scale scatter), one by one / `enqueue_many` | 0.90 s / 0.36 s |
| Redis memory per job, queued / finished with `result_ttl` 3600 | 1,831 B / 7,229 B |
| RQ fork worker, no-op job | 14.2 ms per job (70 per second); `result_ttl=0` 15.4 ms |
| `SimpleWorker` (no fork, no heartbeat, no isolation), no-op job | 3.2 ms per job (309 per second) |
| A chain of 200 jobs, each depending on the last | 16.2 ms per link, all finished |
| 75 jobs (3 stages of 24 units 400 ms each, a gate job per stage that depends on its units, the next stage depending on the gate) queued up front | 0.06 to 0.09 s; 51 deferred |
| Same run, 3 burst workers | 23.96 s: after the first barrier all later units ran on one worker (8, 8, 8 units per worker in the first stage, then 24 and 24 on one) |
| Same run, 3 workers with `max_idle_time=3` | 10.59 s (ideal 9.60 s, serial 28.8 s) |
| Cancel that run during the second stage: cancel waiting and deferred jobs, `send_stop_job_command` the running ones | 44 canceled, 4 stopped, 27 already finished; quiet 0.48 s after the request |
| A gate that depends on N no-op units, 3 workers, N = 200 / 500 / 1,000 | 2.3 s / 10.1 s / 45.4 s to finish (quadratic); queuing the gate itself is cheap (0.4 s for 4,000 dependencies) |
| Two stages of N units, no gate: each unit adds its index to a Redis set, the one that completes the set queues the next stage with `enqueue_many`; 3 workers; N = 1,000 / 4,000 | 14.3 s / 55.0 s, 7.1 and 6.9 ms per unit (linear; includes a new connection per unit and the fork) |
| `kill -9` of worker and horse under a running job | job stays `started` about 90 s, then returns to the queue and retries (when something polls the registry) |
| Dead worker in `Worker.count` after `kill -9` (worker TTL 20 s) | still counted after 60 s, and after its key expired at 80 s |
| `kill -9` of Redis shortly after 100 writes: distribution defaults / `appendonly yes`, `everysec` | 0 of 100 keys / 100 of 100 |
| Default save points, `CONFIG GET save` | `3600 1 300 100 60 10000` |

## 4. Options for the task process

What Boss decided for PERF.86 (2026-10-14 13:24Z) and what the design must therefore do:

1. Each scatter, and each phase of anything, is its own queue item, run one after another.
2. A "cancel now" that stops a running job.
3. Progress bars and ETAs per stage; per-stage timing records; the job log with its download link.
4. An admin starts or cancels a run; no automatic resume after an interruption.
5. Several workers on one job are fine.
6. Plain Redis queue items instead of the task and subtask logic.

Standing wishes it must not break: a worker count the admin can change (ADM.15), an interactive lane so edits do not wait behind bulk work (PERF.40, PERF.78), a banner for every visitor while a job changes the galaxy (UX.3), and the same galaxy for the same seed at any worker count.

| | A. The run is a list of stages in Redis, each stage's units RQ jobs | B. One RQ job per stage, a process pool inside | C. Hand-written Redis lists or streams | D. MySQL queue | E. Keep today's design, fix defects |
|---|---|---|---|---|---|
| Each scatter is a queue item | yes, and so is each chunk of layers | yes (the whole stage) | yes | yes | no (one run, many tasks in the supervisor) |
| Plan visible and durable in Redis | yes (run record; the current stage's units) | stages only | yes | in MySQL | no |
| Cancel now | stop command, 0.5 s measured | stop command plus killing the pool's process group | yours to write | row flag plus a kill | three paths, as today |
| Several workers on one stage | native | needs a pool inside the job (what we have now) | native | native | native |
| Replaces MySQL job tree and jobs-directory files | yes | mostly | yes | partly | no |
| New code | smallest of the full designs: library calls plus about a few hundred lines of glue [R] | smallest overall | largest | large | none |
| Risk | workers must outlive barriers (solved, section 5.4); a failed unit stalls its stage until `on_failure` ends the run; `depends_on` must not be used for "after all N" | the pool inside is the old complexity again | reimplements registries, heartbeats and retries | puts queue traffic on the busiest resource | the supervisor stays the single point of failure |

**A** matches Boss's decision as written: queue items in Redis, run in order, several workers on one stage. **B** is the smallest step if Boss would rather change as little as possible: a stage is one RQ job, a single worker per bulk queue runs stages in queue order (which also serialises runs without a lease), and the stage keeps its process pool. It fails requirement 5 only in the sense that parallelism is not a Redis feature there. **C** and **D** are not recommended: C rebuilds RQ, D moves the traffic onto MySQL. **E** keeps the defects of F1 and F2.

## 5. Recommendation: option A

### 5.1 Shape

A run is declared once, from the stage list that already exists (`generation/stages.py`: `skeleton`, `phenomena`, `mass`, `luminosity`, `sectors`, `link`, `backfill`, `population`, `settle`, ...), and stored as a run record in Redis that lists every stage with its parameters. The first stage's units are queued at once:

```
run record  planetgen:run:<id>   (hash: title, who, database, argv, seed, state, stage list, current stage, times)
stage k     N_k unit jobs   run:<id>:s<k>:u<i>   (JSON arguments: stage, index, layer or sector numbers)
            planetgen:run:<id>:s<k>:done   a set of the indexes that have finished
```

A unit that finishes adds its index to the stage's `done` set. The unit that makes the set complete (`SCARD == N_k`, and `SET ... NX` on an `advanced` flag so only one unit does it, even if a retried unit runs twice) records the stage's end, checks the run is not being cancelled, and queues the next stage's units with `enqueue_many`. The last stage's last unit marks the run done. This is the classic fan-in with an atomic counter; it is linear in the unit count, and it avoids RQ's `depends_on` for "after all N", which is quadratic (section 3).

A unit is the smallest schedulable piece: a chunk of scatter layers, a chunk of sectors, a backfill block. A one-shot stage (`link`, `population`, `settle`) is one unit. Chunks are static, decided when the stage starts (as built, PERF.97: about 512 units by weight, at most 64 shells and 4,000 units, so it never depends on the worker count), which replaces the adaptive, time-sized chunking of `submit_each`. The number of units per stage is capped (about 4,000; a larger sector fill uses bigger chunks), because Redis memory and queuing cost follow the unit count, not the object count: 4,000 queued units are about 7 MB [C], and at 14 ms of fork cost per job the overhead is under 1% for a unit that runs a few seconds. The run's seed is in the run record, and `task_seed(run_seed, (stage, index))` seeds each unit, so the result does not depend on which worker takes it or in what order.

Stages run one after another because a stage's units are only queued when the previous stage is complete. Within a stage any number of workers take units in parallel (requirement 5). The plan is durable in the run record, not in the supervisor's memory; the units of a later stage appear as jobs when their turn comes.

### 5.2 What Redis holds

| Key | Type | Purpose | Lifetime |
|---|---|---|---|
| `planetgen:run:<id>` | hash | run record and state | no expiry while live; 30 days after it ends |
| `planetgen:run:<id>:s<k>:done` | set | indexes of the finished units of stage k; also the barrier | with the run |
| `planetgen:run:<id>:s<k>` | hash | per-stage counters: units total, units done, cost done, seconds summed, first start, last end | with the run |
| `planetgen:run:<id>:events` | stream | the log and timeline: stage start and done, unit done, log lines, errors; `XADD` with `MAXLEN ~ 100000` | with the run |
| `planetgen:run:active` | string | the run that owns the machine; `SET NX`, so a second Generate gets "busy" as today | until the run ends; cleared by the last unit, `on_failure` or the stale check (5.6) |
| `planetgen:workers:wanted` | string | the admin's worker count (ADM.15) | none |
| RQ keys (`rq:*`) | RQ | queues, jobs, registries, workers | RQ's ttls |

Progress for a bar is one `HGETALL` of a stage hash (the units add to it as they finish), so a page does not read the stream. The event stream is for the log view, its download, and a live tail by `XREAD` with the last seen id. The ETA is the cost-weighted estimator already designed in [performance-eta-queue-and-caching.md](performance-eta-queue-and-caching.md) 2.3, as a pure function of those counters.

### 5.3 What MySQL keeps

Only history and statistics, written once per stage by the unit that completes it: the existing `generation_stage_runs` and `generation_stats` rows (per-stage timing, requirement 3). The per-task rows, the node tree and the lease go: `work_jobs`, `work_job_args`, `work_tasks`, `work_lease` and their migration (one Alembic revision, Foundations lane 1 only). Per-object timing (PERF.92) stays its own item.

### 5.4 Workers

- **Workers must survive barriers.** Run them with `max_idle_time` (60 s is plenty; 3 s matched the ideal in the test), not `--burst`. A burst worker exits when the queue is empty, and the queue is briefly empty at every barrier, before the last unit has queued the next stage [C].
- **Still no service to install.** `ensure_workers(queue)`, as `api_jobs` has it today, starts workers when work is queued, up to the wanted count, and is called by whoever queues (the CLI, the web) and by every status poll, so a worker that dies is replaced. Gates call it too, which keeps the pool topped up across stages. Workers exit after the idle time; nothing runs when there is no work. If Boss ever wants always-on workers, `rq worker-pool` under systemd or launchd is the upgrade; the queue design does not change.
- **Count liveness correctly** (F4): live means the process id on this host exists and the worker's key exists; also run RQ's `clean_worker_registry` before counting.
- **Keep** the lowered priority, the pre-import warm-up (PERF.42) and the stale-release restart in `cli/worker.py`.
- **One queue per worker** (F6): `planetgen-bulk` workers serve only that queue, `planetgen-interactive` workers only theirs. A reserved interactive worker (PERF.40, PERF.78) is then a separate, mostly idle worker, not a priority order, and every worker keeps the reliable `BLMOVE` path.
- **Unit length** (PERF.155): a unit is meant to run no longer than `queue.unit_ceiling_seconds` (300 s) unless its stage is one unit by nature (`CommandSpec.whole`), which bounds how long a worker is told to finish its unit before it leaves, and how long a cancel waits. A stage whose longest unit ran longer is written to the run's log; see `docs/design/generation/generation-performance.md`.
- **Worker count** (ADM.15): the Queue page writes `planetgen:workers:wanted`; `ensure_workers` starts up to it, and lowering it sends `send_shutdown_command` to the extras, which finish their unit first. "Ludicrous Speed" is a number, not a mode.

### 5.5 Cancel, pause, failure

- **Cancel now:** set the run's state to `cancelling` (the unit that completes a stage will not queue the next), `cancel()` every waiting job of the current stage, `send_stop_job_command` every started one, clear `planetgen:run:active`, write the final event. The job ids follow from the run id, stage and indexes, so no list needs keeping. Measured at 0.5 s for a 75-job run [C]. Units are idempotent and transactional, so a stopped unit leaves nothing half-written (2.5).
- **Pause** is not in Boss's list; if wanted later it is a state the units check at start (they requeue themselves while paused) and costs no new machinery. Left out.
- **A unit fails:** its index never joins the `done` set, so the stage never completes and the run cannot move on; an `on_failure` callback marks the run failed, cancels the stage's other waiting units and records the error. `Retry(max=1)`, interval-less, stays for worker deaths (RQ retries a job whose horse was killed; a stopped job is not retried); the set makes a repeated unit harmless.
- **A hard-killed horse (OOM):** `work_horse_killed_handler` marks the unit and the run, so the failure is seen at once instead of after the 90 s heartbeat expiry (F5).
- **A dead supervisor** is no longer possible: there is none. The CLI is a client that queues the run and follows its events; closing it does not stop the run. `planetgen galaxy` without `--detach` follows and prints; Ctrl+C detaches, a second command (or the page) cancels.
- **Timeouts:** give every job a `job_timeout` (a multiple of its stage's expected time) instead of `-1`, so a hung unit ends and shows as failed.

### 5.6 Stale runs

`planetgen:run:active` is cleared by normal endings. For the abnormal one (Redis restarted without persistence, everything killed), the Queue page shows the run as stale when no event arrived for `max(120 s, 3 * the longest unit so far)` and no worker is busy, and offers "Clear" (an admin action, as Boss wants for resuming and cancelling).

### 5.7 What goes

`queue/work.py`'s `WorkQueue`, `_ControlStore`, the job tree and the lease; `RQExecutor`, `RQFuture` and `Channel` in `redisqueue.py`; `web/jobs.py`'s files and `web/job_runner.py`; `queue/progress_file.py` (the page reads Redis); `submit_each`'s adaptive chunking; the in-process "no Redis" fallback (Redis becomes required for every run, including one worker, which `OPS.21` and `DOC.44` already treat as true); `GENERATION_QUEUE`. Those files are about 3,300 lines today (`work.py` 1,711, `jobs.py` 728, `job_runner.py` 420, `redisqueue.py` 322, `progress_file.py` 143; `progress_rate.py`'s 203 is reused) [C, `wc -l`]. A replacement of a few hundred lines is plausible but is an estimate [R]. PERF.29 and PERF.30 are dropped as queue items (their sector-level "what is still missing" check is a database check, not a queue feature, and can stay if Boss wants it); PERF.40 and PERF.78 become the two queues of 5.4; ADM.15 becomes the wanted-count key; PERF.33's estimator is reused.

### 5.8 The Redis server

Add `examples/redis/planetgen.conf` and have `update.sh` check and warn, never edit (the OPS.46 rule for Apache): 

```
bind 127.0.0.1 -::1
protected-mode yes
appendonly yes
appendfsync everysec
save 3600 1 300 100 60 10000
maxmemory 512mb
maxmemory-policy noeviction
```

`maxmemory` is a guard: with at most about 4,000 units queued per stage (section 5.1) the queue is under 10 MB, so 512 MB leaves a wide margin and a runaway shows up as write errors, not as the host swapping. `noeviction` keeps lockouts and queue entries safe. Keep one instance: queue, lockouts, counts and locks are all small; a second instance would split persistence needs that are identical here. A password (`requirepass`, or an ACL user in `redis.url`) is optional while Redis binds loopback and Boss is the only user; document how to set one. Document all of it in the install guides (DOC.44).

### 5.9 Clients and security

- One module builds clients: `redisqueue.client(role)` with a cached pool per process and role. Web and API roles: `socket_connect_timeout=2`, `socket_timeout=5`, `health_check_interval=30`, retry on `ConnectionError` and `TimeoutError` with backoff. Cache roles keep their 0.5 s. Workers use RQ's own connection handling.
- RQ jobs and `Queue` use `serializer='json'`; function names are strings. This removes pickle but not the underlying exposure, since a job still names the function a worker calls; the protection that matters is that only planetGen can reach Redis (loopback or a unix socket, optional password).
- Wait for a result with `job.latest_result(timeout=)` or the event stream, not a 50 ms poll.

### 5.10 The other uses stay

Rate limits, lockouts, the count cache and single-flight need no redesign. Keep their fail-open behavior, and with AOF on the lockouts now survive a restart. The unbuilt epoch counter of PERF.20 and PERF.25 (`INCR planetgen:epoch:<db>`) is still a good idea and independent of this work; the unbuilt progress record is replaced by 5.2.

### 5.11 Suggested items

For the TODO thread to file if Boss agrees. Each is one PR; the numbers in brackets are what they fold in.

1. **Bug (no dependency on the rebuild):** `api_jobs.ensure_workers` counts dead workers; count live ones and clean the registry first (F4). Bugfixes lane.
2. Redis server configuration: `examples/redis/planetgen.conf`, an `update.sh` check that warns when `appendonly` is off or `maxmemory-policy` is not `noeviction`, and the install guides (F7, DOC.44, OPS.21).
3. One Redis client module with timeouts, health checks and a cached pool; replace `connect()` call sites (F9).
4. PERF.86 rebuild, in steps that each leave `main` working, Foundations lane 1:
   1. Queue spine: `planetgen-bulk` and `planetgen-interactive`, `ensure_workers` with idle-timeout workers, one queue per worker, the wanted-count key, API jobs moved to the interactive queue [PERF.40, PERF.78, ADM.15].
   2. Run record, stage counters, event stream and the ETA function, as a library with a fake clock; no UI yet [PERF.33].
   3. One stage on the new spine end to end (the mass scatter): static chunks, the done-set barrier, `planetgen galaxy` follows the events, `--workers 1` uses the same path. **Built (PERF.97):** `queue/units.py` and `run_plan.scatter_shells_unit`. A unit's arguments are the database's name, host, port and user (the worker reads the password from its own settings, so none is queued), the shells, their expected stars, the sectors to leave out, the floors and the seed; the result and the unit's progress are hashes in the run's keys (`runrecord.results_key`, `progress_key`) read by the starter, which replaces the `Channel` for this stage. The cut is `units.static_chunks` (about 512 units by weight, 64 shells at most, 4,000 units at most), so it never depends on the worker count; `--workers` is the most workers the run starts (`ensure_workers(count=)`). The starter ends the run `cancelled` on Ctrl+C and `failed` on any error, which frees the machine; the Queue page's "Clear" (PERF.99) is for a starter that was killed outright.
   4. The remaining stages and the sector and backfill units.
   5. Generate page, Queue page, banner and job log read Redis; cancel through the stop command; log written to the jobs directory at run end for the download link [UX.3, ADM.22].
   6. Delete the old machinery, with the schema migration for the `work_*` tables [PERF.29, PERF.30].
   7. Rewrite the notes and docstrings that describe the old model [DOC.20, done].
5. A health check on the admin page: Redis answering, AOF on, memory policy, `used_memory`, last save status, worker count (monitoring, 2.3). Optional.

### 5.12 Every stage a command, one wrapper (Boss, 2026-10-17, 4:12 PM EDT)

Boss: every part of the generate process, all the scatter phases included, is an individual command, and a general process wrapper makes the list of the commands that need to run and puts them into the Redis queue individually. The wrapper builds the command list, each command is one queue item.

This is the shape of 5.1 made explicit. A command is a kind in the command registry (`queue/commands.py`, `REGISTRY`) with JSON arguments and a JSON result, and can be run alone. The wrapper is a pure function from the arguments of `planetgen galaxy` or `planetgen plan` to the ordered stages, each a list of commands with predicted cost; it writes the run record and queues the first stage, and the barrier of 5.1 queues each next stage. A chunk of a chunked stage is its own command. Built by PERF.118 (the wrapper), PERF.106 to PERF.108 (the stages as commands) and PERF.119 (every start path goes through the wrapper).

Decided (Boss, 2026-10-17, 4:15 PM EDT): from the command line Boss runs `plan`, then `scatter`, and each scatter kind is its own command even where the root command is the same: a scatter command for only mass, only luminosity, only phenomena, only the nebula field, only the backfill, never one command that groups them. Each command has its own parameters, set by the UX when it runs the commands through the API and by Boss on the command line. When a new generate job is called, the wrapper only selects the commands that will run and never touches the others; they run in sequence; the progress bars are dynamic to the current context and show only what is going to run.

Governing principle (Boss, 2026-10-17, 4:18 PM EDT): modular commands, each command does one and only one thing and does it exceptionally well. A scatter for mass does that and only that and never worries about anything but what it does. No command knows about another stage, about what ran before it or what runs after it, or about the order. The list of steps is built from what the user requested and handed to the wrapper, which handles the sequence, so the order lives only in the wrapper's step list, kept as data, and the sequence can be changed without touching any command.

Scope (Boss, 2026-10-17, 4:19 PM EDT: "redesign ALL the commands like that so we can build whatever sequences we want"): the principle holds for every command and entry point, not only the generate stages. The inventory, with what each command does today and the item that reshapes it, is [command-inventory.md](../architecture/command-inventory.md); the umbrella is PERF.120.

As built (PERF.118, first slice). `queue/commands.py` holds the registry: every stage of a `galaxy` or `plan` run has a kind (`scatter-mass`, `scatter-luminosity`, `scatter-phenomena`, `scatter-backfill`, `scatter-band`, `plan-skeleton`, `generate-sectors` and so on; `planetgen commands` lists them). A kind is either queued (`unit`, the function a worker runs from the unit's JSON arguments, and `expand`, which cuts the stage into units from its parameters and the database alone) or inline (still done by the process that started the run, until PERF.107, PERF.108 and PERF.121 to PERF.128 give it units). `commands.build_steps` is the pure function from parsed arguments to the ordered `Step`s (key, kind, label, flat JSON parameters): exactly the stages the request selects, never a skipped one, the same list for the same arguments. `queue/wrapper.py` starts the run: it stores that list in the run record (`steps`), declares a stage for each queued step, expands the first stage, records the cut and a JSON `plan` of it (`<PREFIX>:run:<id>:s<k>:plan`) before queueing its units, and follows the units as the starter. The unit that completes a stage calls `wrapper.advance`, which reads the next step from the run record, expands and queues it (a stage of no units is passed over), or finishes the run after the last; no command knows the order. The skeleton is on it as one unit (`plan-skeleton`: `run_plan.plan_skeleton_unit`, whose parameters are the database, the seed asked for as 32 hex digits and the plan's shape options, `run_plan.skeleton_params`; the unit returns the lines the run says about it for the process that started the run to log, since a worker's output reaches nobody, and a parameter it refuses ends the plan with its message). The three shell scatters are on it (PERF.97 and PERF.106), and so are the shells of the phenomenon scatter (`run_plan.expand_scatter_phenomena`, `scatter_phenomena_unit`: the same rows as before at any worker count; its clearing, the special phenomena, the nebula field and the indexes are still done around the stage by the starter until PERF.121 and PERF.122, so the stage's row stays the starter's): the massive-star pass, the luminosity pass and the band (`run_plan.expand_scatter_mass`, `expand_scatter_luminosity`, `expand_scatter_band`, all run by `_scatter_on_queue`, with `scatter_shells_unit` as the unit of each). Each pass is cut into the same units at any worker count and draws the same stars as before the move; the luminosity pass reads the sectors the massive pass marked when its stage starts, between the passes. The clearing and recording around the passes stay in the starter (`scatter_bright_stars`, `add_bright_star_band`) until PERF.119 and PERF.121 move them. Three things every stage reuses came with them (PERF.106). A unit's `job_timeout` is `redisqueue.job_timeout` of its stage's expected seconds (the recorded time of the stage with its settings, `wrapper.expected_seconds`, also the stage's ETA before any unit reports): `redisqueue.TIMEOUT_MULTIPLE` (10) times them, never less than 15 minutes or more than 6 hours, and 6 hours for a stage with no recorded time; no job queued to RQ has the limit -1 (PERF.144: `redisqueue.TimedQueue` refuses one, and the Generate page's jobs, the API's jobs, the work queue's tasks and the tile prefetch get a limit from their own estimates). `units.horse_killed`, the workers' `work_horse_killed_handler`, fails the run at once when a unit's worker process is killed (out of memory or `kill -9`), naming the unit, instead of when the worker's heartbeat runs out; RQ's retry then finds the run over and does nothing. The unit that completes a stage calls `wrapper.finish_stage` before the next stage is queued: the command's `finish` (`CommandSpec.finish`, the shell scatters' `run_plan.finish_scatter`) records the units' speeds in `generation_stats` and returns the stage's metrics, and the stage's `generation_stage_runs` row (the `row` the starter put in the run record with `stages.hand_over_open_stage`) is written with the seconds from the stage's start to its last unit's end. The starter no longer writes that row, but still adds the stage to the whole job's totals (the `job` and `scatter` rows). A stage's seconds count from before it is cut into units, so reading the galaxy is part of them. A stage whose units are only part of what the stage does keeps its row in the starter, so its recorded time still covers all of it: the band, which tops up backfilled sectors after its units, until PERF.108 moves that work too. `planetgen command-run KIND ARGS` runs one of a command's units alone. The map luminosity fill and the habitability re-score are on it too (`fill-map-luminosity` and `rescore-habitability`: a unit per range of ids, `run_plan.expand_fill_map_luminosity` and `expand_rescore_habitability`, whose units write only their own range and so give the same rows at any worker count), and the nebula field is two commands: `clean-nebula-field` (one unit: the dry-run check of the stored field clouds and, with `apply`, the delete and the containment refresh; it returns the lines the run says) and `draw-nebula-field` (a unit per slab of field cells, `run_plan.expand_draw_nebula_field`; the starter clears `nebula_field` first and links the stored clouds back afterwards, as it also does inside the phenomena stage, where the draw is a run of its own after the phenomena run). PERF.202 moved the one-step starters onto the shared run, family by family: the skeleton, `clear-scatter` and `clean-nebula-field` (first family), and the map luminosity fill, which is a `MapLuminosityFollower` stage followed by `BlockLightsFollower`, and the habitability re-score, a `HabitabilityFollower` stage whose block totals are made by the command's `conclude` (`run_plan.conclude_rescore_habitability`). The other stages are inline until their items move them over.

As built (PERF.107, the sectors). The three stages of a `galaxy` run that generate sectors, the starting sector (`generate-start-sector`), the neighborhood (`generate-neighborhood`) and the rest (`generate-sectors`, also the shell, band, wedge, span and course modes), share one cut and one unit (`run_galaxy.expand_generate_sectors`, `generate_sectors_unit`): the starter still enumerates the sectors, works out which are empty, estimates their cost and orders them as before, then submits them as flat `shell, band, wedge` triples in the order to generate them (`run_galaxy._submit_batch`), and the cut is about `SECTOR_UNITS_MAX` (4,000) static chunks of consecutive sectors, the same at any worker count. A unit resolves each sector's density from its sector profile (`_BatchDensity.resolve`), generates and saves it in its own transaction (`generate_and_save_sector_at`) and returns, per sector, what the starter logs: the starter follows the units' results, logs each saved sector, counts it and adds its time to the speed stats, so the console reads as before. Every draw for a sector comes from its own seed (GEN.39), so the stars, planets and moons are the same rows as the old process pool wrote, at any worker count (checked on the toy skeleton against `main`). The link stage (`link-neighbours`, chunks of `LINK_UNIT_SECTORS` (500) sectors since PERF.113, which take the neighbour lock in turn) and the settle stage (`settle-sector-paths`) is a unit per `SETTLE_UNIT_SECTORS` (100) sectors since PERF.115 (the expansion finds the sectors to settle once; each unit saves its chunk's paths): they work on the sectors saved since the run started, found from the database's clock at the start (`since`), so a unit needs nothing the starter held in memory. A failing unit stops the run, which raises `units.RunStopped`, where the old pool re-raised the sector's own exception; the link stage still only warns, and the sectors already saved stay saved, so a rerun resumes. The galaxy job tree no longer has a queue node of sector tasks: the sectors are units of the run record on Redis.

As built (PERF.108, the backfill). The backfill stage of a `galaxy` run is the `scatter-backfill` command: the starter passes the database, the database's clock at the start of the run (`since`) and the worker count, and the stage's cut (`run_galaxy.expand_backfill_bright_stars`) finds the sectors to draw once, around the centers of the sectors the run created (`backfill_targets`: the mass steps outward, the band mass limits, the filled sectors and the sectors already taken that deep), and cuts them in address order into static chunks of at least `BACKFILL_UNIT_CELLS` (10,000) sectors and at most `BACKFILL_UNITS_MAX` (2,000) chunks, the same at any worker count; a stage with nothing to draw is one empty unit, so its row is still stored. A unit (`backfill_bright_stars_unit`) draws its chunk with the same function as before (`run_plan._draw_sector_masses`, a few hundred sectors a transaction, each rechecked under its lock), every sector from its own seed, so the stars are the same rows whichever sector reached a sector first and at any worker count; it returns the stars placed per shell and kind for the starter's lines. The starter keeps the bar of the sectors visited and the per-shell lines. The backfill around one requested sector (`ensure_sector_generated`, the neighborhood estimate) is a single in-process call of the same draw, not a run.

As built (PERF.121 and PERF.122, one root command per kind of work). `planetgen plan` is the skeleton and nothing else (`run_plan.run_plan`: `plan_stages` is the one `skeleton` stage; the new outline also empties the bright stars, the phenomena and the nebula field drawn on the old one). `planetgen scatter` takes exactly one of `--mass`, `--luminosity`, `--phenomena`, `--nebula-field` or `--band --down-to N`, and only the options of that kind (`generate.validate_scatter_args`); `run_plan.run_scatter` runs that kind's stage (`stages.scatter_stages`) and the wrapper's step list for it is that one step (`stages.SCATTER_STAGE`, the nebula field drawn by `draw-nebula-field`). The three that fill or clean (`fill-map-luminosity`, `rescore-habitability`, `clean-nebula-field`, `stages.SINGLE_STAGE_COMMANDS`) are roots of their own and take no setting of a scatter. `commands.STAGED_COMMANDS` names the six roots whose stages the wrapper turns into steps (and `galaxy`'s). A scatter draws only where none of its kind is stored (`run_plan._refuse_stored_scatter`, `store.scatter_stored`), and there is no redo flag: Boss chose (F7) a separate `planetgen clear-scatter` that empties one kind (`--mass`, `--luminosity`, `--phenomena`, `--nebula-field`; the `clear-scatter` command, one unit, `run_plan.clear_scatter_unit`), run before the scatter. The luminosity pass works at the mass limit the mass pass stored and refuses to run with none; the mass clear refuses while the bright stars, which are split from the massive stars at that limit, are stored. The Generate page builds a list of argv steps, one per scatter (`generate_page.scatter_steps`), where it once passed one `plan` run a mode flag; `galaxy --then-scatter` still lists the scatters as steps of the galaxy run. The Generate page's Redo scatters form builds the clears then the scatters as steps. The same galaxy comes out of `plan` followed by `scatter --mass`, `--luminosity` and `--phenomena` as out of the old combined `plan`.

As built (PERF.119, first slice: `sector` and `system`). `planetgen sector` and `planetgen system` are runs of one command each: the entry point (`run_sector.run_sector`, `run_system.run_system`) draws the run's seed, builds the step list with `wrapper.steps_for` (`commands.ONE_COMMAND_RUNS`) from what was asked and calls `wrapper.run`; it does no generating itself. `generate-unplaced-sector` cuts the sectors into static chunks of consecutive numbers (`run_sector.expand_unplaced_sectors`, at most `UNPLACED_UNITS_MAX` units) and each sector draws from a stream of its own, seeded from the run's seed and its number, so the sectors are the same at any worker count; the unit (`unplaced_sectors_unit`) generates and saves them and returns what the starter logs. `generate-system` is one unit (`run_system.system_unit`) that renders (`--output`) or saves the system and returns the page or the id and, for a forced body the system could not hold, the warning the starter refuses or logs by `--strict`. A `system --output` run saves nothing, so it is a maintenance run and does not take the generation slot. The old process pool path for these two commands (`run_common._work_queue`) is deleted. The Generate page, the API's jobs and the one-off reset and check-math steps follow in the next slices.

As built (PERF.119, second slice: the API's jobs). The API jobs that change the galaxy are the starters of runs, not stages of their own: `api_jobs.generate_neighborhood` is `run_galaxy.generate_sector_neighborhood`, which queues the neighborhood's sectors as a run of `generate-sectors` (`_submit_batch`, the helper the `galaxy` command uses), then the bright-star backfill around the requested sector as a run of `scatter-backfill` (`center_pc` in its parameters in place of the sectors made since the start) and the paths of the new sectors as a run of `settle-sector-paths` (`sector_ids` in its parameters in place of those made since the start; `expand` false settles exactly the named ones); `api_jobs.regenerate_sector` deletes the sector and runs `run_galaxy.regenerate_sector_cell` (the sector and its backfill) and the settle of the new sector and the old neighbours; the settle run alone has no job of its own (TEST.531 removed `api_jobs.settle_sectors`, which only a test helper reached). The job is on the interactive queue and follows the runs from its worker, which answers `GET /api/jobs/<id>` as before. A job takes its turn the way every generation run does: a run that finds the machine taken waits in line for it (`runwait.start_when_free`, PERF.153, up to `queue.wait_minutes`). The sector the map generates on a visit (`ensure_sector_generated`, inside the web request) is not a job and still generates in the request; whether it becomes a run is an open question. The API's other jobs (edits, uploads, page builds, the one-off page's commands) start no stage and stay as they are.

As built (PERF.119, the one-unit commands). `reset-galaxy` (`cli.reset.reset_galaxy_unit`, a generation run), `check-math` (`physics.mathcheck.check_math_unit`, a check run) and `check-datasets` (`cli.datasets.check_datasets_unit`, a check run) are registered commands of one unit each, with the database or the options as plain JSON in the step's parameters. Their verbs, `python3 -m planetgen.cli.reset` (after the typed confirmation or `--yes`; `--dry-run` starts no run), `planetgen check-math` and `python3 -m planetgen.cli.datasets status`, build the step and run it through the wrapper, and print what the unit returns (a worker's output reaches nobody). `datasets fetch` is a setup download and stays direct. The Generate page's reset and data-check steps still run these verbs as steps of the job until PERF.203 puts the page on one wrapper run.

As built (PERF.172, `planetgen galaxy` as one run). The galaxy is one run in the queue over one step list. The process that starts it plans the sectors in the mode asked for (the estimate and the questions of `_plan_sectors`: random start, single address, neighborhood, shell, band, span, course or block), and then `generation/followed.run_followed` puts the stages in one `wrapper.run`: the sector stages (`start` and `neighborhood`, or `sectors`, or `sectors` and `neighborhood` for `--wedge` with `--radius-pc`), `link`, `mass`, `luminosity`, `phenomena`, `phenomena-index`, `nebula-draw`, `block-lights`, `backfill`, `population`, `territories`, `settle` and `run-blocks` (`stages.galaxy_stages`). Each stage's `Follower` says what its process shows and counts as the plan arrives and the units come back (bars, lines of the log, the summary); it does none of the stage's work. The work that is not a unit moved into the commands: a command may name a `prepare` (once, in a worker, when its stage opens, before the cut: empty the table the units refill) and a `conclude` (once, in the unit that completes the stage, before the next stage opens: store the record the later stages read, build the indexes, link the nebula field); the starter waits for the conclusion before it closes the stage's display. A hook or the opening of a stage that fails fails the run at once, naming the stage ("stage X could not start", "the end of stage X failed"), where it used to be lost as a unit done when Redis failed it. A stage with nothing to do (no sector saved since the run started, nothing to backfill) stores its row as skipped, so it adds nothing to the estimates. The same seed gives the same galaxy at one worker and at four: every draw still comes from the seed of the thing it draws, and the cuts do not depend on the worker count. The estimate (`planetgen galaxy --estimate-only`, the run's listing) shows the whole list. A link that fails now fails the run, where it only warned; a failed galaxy prints the `planetgen link --since` line to run it again.

As built (PERF.182, finish an incomplete run). A run can be carried on after it stopped. The run record already holds the steps (a command and its complete parameters per stage) and each stage's `started`, `ended` and `concluded` clocks; it now also holds the galaxy's seed, the build (`galaxy.version_key`) and, for a run that carries another on, `finishes`. `generation/finish_run.py` reads the record of the newest generation run that stopped (a live record whose holder is dead is marked failed first), refuses when the seed, the build or the database is not the one the run recorded, when the run was already carried on, or when a stage it would run is not in `FINISHABLE`, and gives the steps of a new run: a stage the old run finished carries `finished_by` and is *carried over* by the wrapper (`wrapper._carry_finished_stage`: it opens and ends at once with no unit and stores the earlier conclusion), the stage that stopped carries `again` and runs from its start with its `prepare` (the scatters empty what they rewrite; the sector stages skip the sectors already saved, `run_galaxy.generate_sectors_unit`), and the stages after it run as planned. `REWIND_WITH` runs the scatter before a stopped `finish-phenomena`, and the mass scatter before a stopped luminosity scatter that draws both passes, again with it: their rows are written into what only the earlier stage's `prepare` empties. `planetgen finish-run` follows the new run as the other commands do (`RecordedStage` stands in for the followers, which hold the stage state only the starting process knew); the Generate page's `finish_run` job is that command, then the later steps of the job the stopped run was part of (`jobs.remaining_steps`).

As built (PERF.213, a parameter worked out when its stage opens). A step's parameters are complete before the run starts, but a value that belongs to another stage cannot be: the stages that work on the sectors a run saved (`link`, `backfill`, `settle`, `run-blocks`) need the database's clock when the first sectors stage opened (`since`), not the moment the followers were built (a run that waits for the machine, or one whose job begins with other stages, would otherwise take in sectors that are not its own). Such a parameter is a *deferred* one, `wrapper.opened_by(key)`, a dict with a `deferred` name that stands in the step until its stage opens; the run record keeps the step with the marker. When a stage opens, `wrapper._expand` records the database's clock for it (`runrecord.set_opened_clock`, only for a stage some step refers to) and gives the command's `prepare` and `expand` the parameters with every marker replaced (`wrapper.resolve_params`, through `wrapper.RESOLVERS`: marker name -> `resolver(conn, run_id, marker)`; a command that needs another value worked out when its stage opens adds a resolver there). The conclude and finish hooks see the resolved parameters too. `finish_run.plan` replaces the markers whose stage opened in the run being carried on with the clock that stage recorded (a stage run again keeps the clock of its first opening, so the sectors the stopped run saved are still the run's) and leaves the ones whose stage never opened for the new run to resolve. `planetgen galaxy` and `planetgen sectors` take the clock of the first sectors stage this way; the single-verb runs (`planetgen link --since`, a settle with named sectors) pass a clock or a list of sectors as before.

As built (PERF.214, the sectors planned when their stage opens). `run_galaxy._plan_sectors(args, edge_pc)` counts, estimates and asks about the sectors in the process that starts the run, so a run could only be built once the galaxy it reads was there. With `deferred=True` it plans nothing: `_deferred_sectors` gives the same stages (`stages.sectors_stages`: `start` and `neighborhood`, `sectors` and `neighborhood`, or `sectors`), each told only the arguments that name its sectors, and the stage works them out when it opens in a worker (`planned_items`). The shell, band, span, course and block modes and the neighborhood of a centre sector were already cut that way; a single sector (`planned = "single"`) and the random starts (`"random-start"`, from the run seed the plan draws, the one draw the starter's plan makes too) now are. The neighborhood after a random start takes its centres from the plan of the `start` stage (`wrapper.plan_detail("start", "centers")`, a deferred parameter resolved when the neighborhood opens); a single sector's neighborhood takes the address. Nothing is estimated or confirmed (`--estimate-only` is refused): the caller has its own estimate, which is why this is for the one run of a whole Generate job. `finish_run.plan` turns the starts a stopped run picked into the explicit sectors of the stage it runs again, so picking again cannot skip the start it saved and choose another. The stored sectors are the ones the starter's plan makes.

As built (PERF.99, the pages). The Queue page lists the runs in Redis (`queue/runview.py` turns a run record into the words and numbers every page shows; `queue/runcontrol.py` acts on one). `/admin/queue/<run id>` shows the run's facts, its stages and the tail of its log; cancel sends the stop command to the started units and takes the waiting ones off the queue (pipelined, 0.1 s for 4,000 units), a run that shows no sign of life or whose follower is gone can be cleared (ended `failed`, with who cleared it), and finished runs can be deleted. The wrapper tees the `planetgen` logger into the run's event stream and writes `<jobs dir>/runs/<run id>.log` when the run ends, the last 20 kept, for the download link. A run started under a Generate page job (`PLANETGEN_JOB_ID`, set by the job runner) is listed under that job (`runrecord.job_runs`), the job page shows its runs, and cancelling the job cancels them at once. The API is version 6: `/api/admin/work` lists runs, `/api/admin/work/<run id>` and `/log` read one, `/control` takes `cancel` or `clear`, and the job tree routes and the lease clear are gone. Left for PERF.100 with the `work_*` tables: `progress.json` and `progress_file.py`, the `work.py` writers, the job runner's step estimates from `work_jobs`, and the update guard's reads of them.

As built (DB.29, the update guard). `update.sh` and `planetgen.cli.migrate` (the only caller of `migrate_database`) take the guard (`admin/updateguard.py`) before they change anything. It lists the live runs of the run record (`RunRecordSource`: a run in state running or cancelling, cancelled through `runcontrol.cancel_run` with who `update`), the Generate page's job, the jobs on the shared queues (the units of a run are left out, `UNIT_JOB`: the run lists them), the orbital update and the warm-up, and goes on only with an agreement: `--yes` (it replaces `--cancel-running`, no shim) or a y typed at the question. Nobody to ask is a no: `started_unattended` (systemd by `JOURNAL_STREAM`, or `INVOCATION_ID` with no terminal on stdin; cron by an ancestor process named cron, crond, anacron or atd) and a run with no terminal never ask and refuse with exit status 1, and the question waits `PROMPT_SECONDS` (30) for a line (`answer_within`, `select` on stdin) and takes silence as no. The old work table source (`WorkRunSource`) stays until PERF.100 deletes the tables. `install.sh` does not take the guard.

What the old command structure becomes (Boss, 2026-10-17, 8:06 PM EDT: the rework removes or adapts every old entry point, with no compatibility shims). The registry has one path: a command is either queued (it has a `unit` and an `expand`) or inline (it has neither, and the process that started the run still does it). A follow-up moves a command over by giving its registry entry a `unit` and an `expand` and deleting its inline function and its call in the same change, so no stage ever has an old and a new way side by side. What each follow-up replaces:

| Old entry point | Replaced by | Item |
|---|---|---|
| `units.HANDLERS`, `units.run_stage`, `run_plan._scatter_shells_on_run_record` and the `on_run_record` switch (the first draft of the new queue) | `commands.REGISTRY`, `wrapper.run`, `run_plan._scatter_mass`; already gone in PERF.118 | PERF.118 |
| `run_plan._scatter_shells` with `work.py`'s task queue (the bright-star pass, the band) | the `scatter-luminosity` and `scatter-band` units; already gone in PERF.106 | PERF.106 |
| `plan`'s flags and in-process stage calls (`scatter_phenomena`, `redo_scatters`; the fills and the nebula field have their commands, PERF.106) | one `scatter` command per kind, the skeleton, the fill and the re-score as their own commands; `plan` becomes `build_steps` plus `wrapper.run` | PERF.121, PERF.122, PERF.106 |
| `galaxy`'s in-process start, neighborhood, sectors, link, backfill, population and settle | the same stages as registered commands; `galaxy` becomes a step list | PERF.123, PERF.107, PERF.108 |
| `sector` with `--population`, `phenomenon` with its settle | `sector` and `phenomenon` only generate; the settle after a phenomenon is `settle --sector-id` (done) | PERF.124 |
| `population` with `--territories-only` and `--rescan` | two commands, `run-population` and `run-territories` (done) | PERF.125 |
| `check-db` behind flags, `warm_map --only`, `benchmark`'s own sequence | a command per check, per cache family, and a step list; `check-db` is the nine `check-*` commands run as a check run (PERF.126, done); `warm_map` is `planetgen warm` and the six `warm-*` commands, run as a maintenance run (PERF.127, done); `benchmark` is `plan` and then `galaxy` handed to the wrapper in scratch databases, its stage timings read from the run records (PERF.128, done) | PERF.126, PERF.127, PERF.128 |
| The Generate page, API jobs and one-off pages starting their own runs | build a step list and call the wrapper | PERF.119 |
| The Generate page, Queue page and job log reading the old tables; the `work_*` tables | read the run record in Redis; deleted with their migration | PERF.99, PERF.100 |
| `update.sh` and migrations starting while jobs run | refuse while the run record shows a live run (done) | DB.29 |

## 6. Decisions for Boss

Each has a default; the default is what a "yes" means.

1. **Option A** (default), B if the smallest step is preferred.
2. **Redis required for every run**, including one worker; the in-process fallback goes (default yes).
3. **Workers: on demand with an idle timeout, no service** (default), against an always-on systemd or launchd worker service.
4. **Turn on AOF `everysec`** and set `noeviction` and a memory guard, through a documented example config and a warning from `update.sh` (default yes; `update.sh` never edits his Redis config).
5. **One Redis instance** (default) against two (queue and lockouts apart from caches).
6. **No Redis password while it binds loopback** (default); an ACL user documented for later.
7. **Stage timing into the existing MySQL stage tables, from the completing unit** (default yes); live task rows in MySQL end.
8. **The log:** the event stream is the live log; at run end it is written to the jobs directory for the download link (default), with the last 20 kept as today.

## 7. What could not be measured

- The production server: all numbers come from a 4-core sandbox with Redis 7.0.15 and no MySQL load.
- A full planetGen run on the new design. The prototypes used `time.sleep` and no-op units, not generation, so they show queueing, barriers, worker lifetime, cancel and per-job overhead, not generation speed. Retry of a unit interacting with the done-set was reasoned, not run.
- The roughly 60 database statements per task are the `work.py` header's own count; they were not re-measured.
- Redis 8 and Valkey were not tested; neither were macOS (RQ's fork worker is not safe there; `SpawnWorker`, RQ 2.2 and later, is the documented replacement) or the idle memory of a pre-warmed worker (195 MB earlier, before the numba work).
- The size of the replacement code is an estimate.
- The end-to-end effect of F4 through planetGen's own `api_jobs` was inferred from the library's behaviour and the measurement, not reproduced.

## Sources

Documentation read for this note (Redis and RQ pages were read from their documentation repositories on GitHub, `redis/docs` and `rq/rq`, `main` and `master` on 2026-10-14; the RQ pages describe master, which is ahead of the released 2.12.0 in places, so features were checked against the installed 2.12.0 package where noted):

- Redis persistence: https://redis.io/docs/latest/operate/oss_and_stack/management/persistence/
- Redis eviction: https://redis.io/docs/latest/develop/reference/eviction/
- Redis security: https://redis.io/docs/latest/operate/oss_and_stack/management/security/
- Redis lists: https://redis.io/docs/latest/develop/data-types/lists/
- Redis streams: https://redis.io/docs/latest/develop/data-types/streams/
- RQ documentation: https://python-rq.org/docs/ (workers, jobs, results, exceptions, scheduling); `rq` 2.12.0 source (`queue.py`, `worker/base.py`, `worker_registration.py`, `intermediate_queue.py`)
- PyPI JSON for `rq`, `redis`, `celery`, `dramatiq`, `huey`, `arq`, `valkey` (versions and dates): https://pypi.org/pypi/rq/json and its siblings
- redis-py production usage (health checks, retry, timeouts): https://redis.io/docs/latest/develop/clients/redis-py/produsage/ (read through a search result)
- Kleppmann, "How to do distributed locking" and Antirez's reply: https://martin.kleppmann.com/2016/02/08/how-to-do-distributed-locking.html and https://www.antirez.com/news/101 (read through search-result summaries, not the originals)
- Sidekiq best practices: https://sidekiq.org/wiki/Best-Practices (search result)
- Celery Redis broker notes (visibility timeout, eviction): https://docs.celeryq.dev/en/main/getting-started/backends-and-brokers/redis.html (search result)
- Task queue comparison: https://judoscale.com/blog/choose-python-task-queue (search result; one benchmark, treat as indicative)
- MySQL locking reads, `SKIP LOCKED`: https://dev.mysql.com/doc/refman/8.4/en/innodb-locking-reads.html (search result)
- Redis 8 and AGPLv3: https://www.phoronix.com/news/Redis-8.0-Goes-AGPLv3 and https://www.infoq.com/news/2025/05/redis-agpl-license (search results); distribution packaging was not confirmed from an official package page

Repository: `src/planetgen/queue/redisqueue.py`, `work.py`, `api_jobs.py`, `cli/worker.py`, `web/jobs.py`, `web/job_runner.py`, `admin/throttle.py`, `api/loginguard.py`, `api/config.py`, `db/countcache.py`, `web/lib/singleflight.py`, `generation/stages.py`, `db/control_schema.sql`, `scripts/deploy-common.sh`, `docs/TODO.md` (PERF.86 and the items named above).

Scripts and raw output: `/mnt/project-files/research/redis-review/` (`b1.py`, `killtest.py`, `phantom.py`, `persist.sh`, `proto.py`, `proto_tasks.py`, `scale5.py`, `fan.py`, `fanrun.py`, `tasks.py`, and the `.out` files).
