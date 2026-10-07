# What the site starts, and what moves to the work queue (PERF.19)

The audit PERF.19 asks for. Boss (2026-10-01 22:13Z): "I want
_everything_ that communicates with the API to use the work queue
whenever possible." This lists every path where the API or the website
generates something or writes to the database. For each one it says
where that work runs today and where it goes once PERF.24 puts the queue
on Redis with RQ. The moves themselves are done in PERF.24.

Audited on 2026-10-07 against `main` after the package move (OPS.24).

## Today: two kinds of background work, and the rest in the request

- **Web jobs** (`planetgen.web.jobs`, run by `planetgen.cli.job`). A job
  is a list of command lines (`python -m planetgen.cli.generate ...`,
  `planetgen.cli.reset`) that a detached runner process works through,
  with its state and log in the jobs folder. Only one runs at a time
  (`JobBusy`).
- **The generation work queue** (`planetgen.queue.work`). Inside a
  bulk command, it spreads sectors, layer scatters and backfill blocks
  over worker processes, using leases, heartbeats and reclaim in MySQL
  (`work_lease`). The web never calls it directly. It reaches it only
  through a web job that runs the CLI.
- **Everything else runs inside the request**, in the web server's own
  process, while the browser or API client waits.

## Every path

| Path | What it does | Where it runs today | Where it goes |
|---|---|---|---|
| Generate page (`/admin/generate`): plan, galaxy (ring, column, shell, center, block, random), sector, backfill, population, reset | Bulk generation and the galaxy reset | Web job, then the CLI, then the work queue | **RQ.** Each unit (sector, layer scatter, backfill block) is an RQ job under one parent per run. |
| Galaxy Map and Sector Map Generate buttons (`generatebuttons.js`) | The same forms, posted to the Generate page | Web job | **RQ**, as above |
| Sector page "Generate the neighborhood" (ADM.11) | `planetgen galaxy --center-sector` | Web job | **RQ**, as above |
| Queue page actions (`/admin/queue/action`): retry, requeue, pause, resume, delete | Starts a web job again, or changes queue state | Web job (retries); in request (state changes) | **RQ** for retries; the state changes become RQ calls, still made in the request |
| Generate page size and time estimate | `planetgen ... --estimate-only` as a subprocess, up to 120 s | **In request** | **RQ**, a short job the page waits on with a spinner. It reads the database hard on a big ring. |
| One-off system page (`/admin/generate/system`) | `planetgen system` as a subprocess, up to 120 s, never saved | **In request** | **RQ**, a short job. The page shows the result when the job finishes. |
| `POST /api/sectors/<id>/generate-neighborhood` | Fills every sector within the radius, up to 652 ly (minutes to hours) | **In request**: holds a web worker, and the docstring warns proxies time out | **RQ**: returns `202` with the job's id at once; `estimate_only` stays in the request |
| `POST /api/sectors/<id>/regenerate` | Refills one sector with its bright-star backfill | **In request** | **RQ**: a sector fill is what the queue already runs |
| `POST /api/systems`, `PATCH /api/systems/<id>` with `regenerate` | Generates one star system and stores it | **In request** | **RQ**, with the request waiting a short time for the result. One system is quick, but a forced large system isn't. |
| `POST /api/{planets,moons,...}/<id>/regenerate`, `POST /api/phenomena/<type>/<id>/regenerate`, `POST /api/planets/<id>/class`, `POST /api/moons/<id>/class`, `POST /api/systems/<id>/star` | Rolls one body, phenomenon or star again | **In request** | **RQ**, waited on in the request, as above |
| `POST /api/systems/<id>/wiki`, `POST /api/sectors/<id>/wiki` | Uploads a page to Wiki.js or MediaWiki | **In request**, with an outside HTTP call | **RQ**: an outside server can be slow or down, and RQ retries it |
| Deletes: `DELETE /api/sectors/<id>`, `/sectors/<id>/contents`, `/systems/<id>`, `/{planets,moons,...}/<id>`, `/phenomena/<type>/<id>`, `/facilities/<id>` | Removes rows | In request | **Stays in the request.** A few row deletes in one transaction; queueing them adds a round trip and nothing else. |
| Renames and small edits: `PATCH /api/sectors/<id>`, `/systems/<id>` (name only), `/stars/<id>`, `/planets/<id>`, `/moons/<id>`; `POST /api/sectors`, `POST /api/facilities` | Updates or inserts a row | In request | **Stays in the request**, for the same reason |
| Logins, two-factor, credentials, API keys, lockouts (`/api/auth/*`, `/api/admin/lockouts/lift`, the admin and account pages) | Session and account rows | In request | **Stays in the request.** They must answer the person at once, and the rate limits sit on them. |

No GET route writes anything: map visits and page views generate
nothing. The Galaxy Map's tiles come from `tilecache.py` and the
database only.

## The plan for PERF.24

1. **One queue for everything.** RQ on Redis replaces both web jobs
   (`planetgen.cli.job`, the jobs folder) and the generation work queue
   (`planetgen.queue.work`, `work_lease`). A Generate page run is a
   parent job whose children are the sector, scatter and backfill
   units. The admin pages read status, progress and logs from RQ, and
   ADM.22 streams them.
2. **Long work leaves the request.** The neighborhood route, sector
   regenerate and wiki uploads return `202 Accepted` with a job id and a
   status URL (`GET /api/jobs/<id>`). The website shows the job the way
   it shows a Generate page run.
3. **Short generation runs on the queue too, waited on.** Single-object
   generation (one system, body, phenomenon or star), the one-off system
   page and the estimates are queued, and the request waits up to a few
   seconds for the result. A quick answer still comes back in the same
   response. A slow one answers `202` with the job id, so no web worker
   is held. This puts all generation in one place, with one worker
   count, one log and one set of limits. It is also what PERF.20's
   short-term caching builds on.
4. **Plain row writes stay in the request.** Deletes, renames, facility
   rows and anything about logins and accounts. They are one quick
   transaction each. Queueing them would only add a Redis round trip
   and a second place to fail.
5. **Same results.** A queued run gives the same galaxy as today for
   the same seed at any worker count (PERF.24's test). The estimate and
   disk-space refusal (PERF.3, ADM.33) and the math gate (TEST.68) run
   before anything is queued, so a refused run never reaches a worker.

## The first choice: short generation

Point 3 is the one choice here. Single-object generation is quick
today, so it could stay in the request. Queuing it follows Boss's
"everything ... whenever possible" and keeps all generation behind one
queue, at the cost of a Redis round trip on each edit. Boss approved
the plan on 2026-10-07 with the recommendation: queue it, waited on, as
above.

## How PERF.24 gets built

RQ needs worker processes that wait on Redis. Today nothing runs in the
background between runs: a `planetgen` run supervises its own pool and
exits, and the web starts a job runner per job. That raises the second
choice below. The rest of the build doesn't depend on it.

**Workers (the second choice).** The recommendation is **burst
workers, started by whoever queues work**. When a run, the web app or
the API queues jobs and fewer than `worker_count()` workers are alive,
it starts the missing ones as `python -m planetgen.cli.worker --burst`
at lowered priority. RQ's burst mode exits once the queue is empty.
That keeps today's model: no service to install, nothing running when
there's no work, and the 80%-of-cores cap holds because the count is
checked against RQ's live worker registry. The alternative is standing
worker services that install and update set up (a systemd unit, a
launchd daemon, a Windows service). That is simpler at run time, but it
is one more service on every platform and holds memory while idle.

**Windows.** RQ's normal worker forks, and Windows can't fork. Windows
workers use RQ's `SpawnWorker` with a timer-based timeout instead of
`SIGALRM`, talking to the Redis in WSL2 (OPS.27). CI's Windows job runs
the queue tests that way.

**Steps, one PR each:**

1. `planetgen.queue.rq`: the Redis connection from `redis.url`, the
   queue names, the burst-worker launcher and the worker command
   (`planetgen.cli.worker`). Tests run against a real Redis
   (`PLANETGEN_TEST_REDIS_URL`, which CI already provides) and skip
   without one.
2. `WorkQueue` runs its tasks as RQ jobs. Each task keeps its
   `task_seed`, `on_done` and the progress bars, and the job tree rows
   (ADM.12) stay in the control database as the history the queue page
   reads. The `work_lease` table and the in-process pool go. Done when
   the same seed gives the same galaxy at 1, 2 and 4 workers.
3. Web jobs are RQ jobs. `jobs.start_job` queues the steps instead of
   spawning `planetgen.cli.job`. Cancel and pause use RQ's stop
   command, and the jobs folder keeps only the downloadable log.
   `planetgen.cli.job` goes.
4. The API moves from the table above: the neighborhood route, sector
   regenerate, wiki uploads, the one-off system page, the estimates
   and, if Boss agrees, single-object generation. Each returns `202` and
   a job id when it doesn't finish quickly, with `GET /api/jobs/<id>`
   for its state.

Live progress and log lines over SSE are ADM.22, built on step 3.
