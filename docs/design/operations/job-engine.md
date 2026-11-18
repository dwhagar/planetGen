# The job engine, as built

OPS.252. Written 2026-11-17 (EST) from the code on main after OPS.234, for stages 0 to 4 and the parts of stages 5 and 6 that are merged. It is the one note for how work runs on a planetGen host: what decides, what is stored where, what a job is, and what is still the old queue. It replaces the queue parts of [`job-graph-and-memory-broker.md`](job-graph-and-memory-broker.md) (DOC.113), [`self-tuning-jobs.md`](self-tuning-jobs.md) (DOC.120), [`one-run-per-job.md`](one-run-per-job.md), [`dynamic-queue-chain.md`](dynamic-queue-chain.md) and the worker, lane and lease parts of [`redis-and-the-work-queue.md`](redis-and-the-work-queue.md); those notes keep their history and each begins with a pointer here. The design reasons, the audit (W1 to W12) and the research are in `/mnt/project-files/reports/job-engine-design-2026-11-13.md` and its two companions; this note does not repeat them, it says what exists.

Where this note and a module's docstring disagree, the docstring of a file that changed more recently is right; fix this note.

## 1. What it is, in ten lines

1. **One supervisor process per host decides** what runs, what waits and what stops (`queue/supervisor.py`). It is taken by a kernel file lock, started on demand by any submit or status poll, and exits after `IDLE_SECONDS` (300 s) with nothing to watch.
2. **Every job of a run is a row** in the control database (`engine_runs`, `engine_jobs`, `engine_attempts`; `engine_items` is in each galaxy's own database). The plan is rows: a stage is cut into rows by an `open` job, and ended by `conclude`, `finish` and `release` jobs with edges between them.
3. **A job's done mark commits with its data** and is fenced by an attempt number: a commit from an attempt that was replaced touches no row, and rolls back.
4. **Each job is a fresh process forked from a warm template** (`queue/template.py`, `queue/jobhost.py`). The template is its parent, so the end of every job, and how it ended, is known at once.
5. **The engine owns the loop** of a command that has `pieces` and `do` (`queue/engine.py`): real counts, a cursor kept with the data, and a safe place to stop between two pieces.
6. **Admission is the memory gate and the fair-share scheduler**, applied by the supervisor to the whole backlog (`queue/dispatch.py` over `memorygate.decide` and `scheduler.judge`).
7. **A job is judged on evidence** (`beatjudge`): its count, processor time, bytes and statement. Only the dead and the proven idle are stopped.
8. **Failures are handled by kind** (`queue/failures.py`): passing ones retry with doubling spread waits, killed ones retry alone with a bigger claim, idle ones once alone, lasting ones **hold** the item and the run (`queue/held.py`) while everything else goes on.
9. **Redis carries live display, wake-ups and the run's record** (`queue/runrecord.py`); the facts are in MySQL. The site's own small jobs (API work, slow pages, tile prefetch) are **light jobs** held in Redis (`queue/lightjobs.py`), not rows.
10. **It is proven** by a simulator under a fake clock (`/mnt/project-files/engine/stage0/enginesim.py`, OPS.212) and kill tests with real processes (OPS.213). The full kill matrix and a soak are OPS.246 and not yet built.

## 2. Where each stage stands

| Stage | What it is | State on main |
|---|---|---|
| 0 Ground | Measure before building (section 8) | Done 2026-11-13. |
| 1 Progress contract | `pieces` and `do`, counts, second-level bars, evidence in the beat, guard tests (OPS.191 to 206) | Done. Every command family has pieces, or a written reason it is single. |
| 2 Supervisor takes the watch | The supervisor, the reconciler, runs started through it, the simulator, kill tests (OPS.207 to 214) | Done (PR 6571). |
| 3 Ledger and facts | Ledger tables, plan as rows, fenced done mark, open/conclude/finish jobs, checkpoints, failure by kind, held runs, resume from rows, attempt history, Queue page from the ledger (OPS.215 to 224) | Done. OPS.225 (the proof) and OPS.261 (delete the old retry paths) are open. |
| 4 Central dispatch, RQ out | Template and forked jobs, dispatch rule, bulk then maintenance and check lanes as forked jobs, light jobs, one memory controller, named resources, controls through the supervisor, Queue page "why jobs wait" (OPS.226 to 234) | Done except OPS.235 (RQ removed) and OPS.236 (the proof). RQ is still the path of **runs with no ledger** until OPS.271 (section 7). |
| 5 Planner | A stored plan for every start path (OPS.238), independence harness (OPS.241), overlap by proof (OPS.242), the engine picks the split (OPS.243), Generate page jobs join the ledger (OPS.244) | OPS.238, 241, 242, 243 done. OPS.244: the one-off system page runs as a light job in the warm supervisor process (slice A); a Generate page job is a ledger run of kind `page` (slice B, section 7); the waiting line as rows (slice C) is open. OPS.239 (single-purpose commands), OPS.240 (sectors one at a time) and OPS.245 (proof) are open. |
| 6 Hardening | Kill matrix, trial, MySQL as a signal, learned thresholds, retire typed timeouts, this note (OPS.246 to 252) | OPS.248 done; OPS.251 slice 1 done; OPS.246, 247, 249, 250 open; this note is OPS.252. |

## 3. The ledger and the plan

**Tables** (`queue/ledger.py`; `control_schema.sql`; migration 0261 for `engine_items`):

| Table | Where | One row per |
|---|---|---|
| `engine_runs` | control database | run: kind, galaxy, state (`waiting`, `running`, `held`, `done`, `canceled`), why held, times |
| `engine_jobs` | control database | job: run, stage key, command, item key, arguments (JSON), claim, class, state (`waiting`, `ready`, `running`, `done`, `failed`, `canceled`), newest attempt number, not-before time, the reason it waits (`waited`) |
| `engine_attempts` | control database | attempt: process id, start, end, outcome (`running`, `done`, `failed`, `killed_memory`, `timed_out`, `stuck`, `replaced`, `lost`), error kind and text, stack file, peak memory, processor seconds |
| `engine_items` | the galaxy's database | job that writes to that galaxy: attempt number, cursor, done time; `engine_piece_results` beside it keeps what finished pieces returned |

**The plan is rows** (`queue/planrows.py`). `record_plan` writes a run's whole step list before its first stage opens (what each step needs and makes: `queue/needs.py`). The **cut** of a stage, `record_cut`, inserts one job row per unit, the stage's `open` row (marked `done`, holding the plan summary) and the end jobs, in one transaction: a cut that stopped half way left nothing, and a cut that finished is found by its `open` row and never made twice. A stage's end is three jobs with edges (`conclude` needs every unit, `finish` needs `conclude`, `release` needs `finish`, and `release` queues the next opening), so a stage advances inside the transaction that finishes a job. The Redis done set and the `advanced` flag are no longer the truth.

A run with **no usable ledger** (no control database, or one with no tables yet: `generate system --output`, the update's math check before its migration) keeps its plan in its Redis record behind the same functions (`planrows.ledgered`). Its units have no rows and no attempts.

**The fence** (`queue/attempts.py`, `ledger.begin_attempt`, `ledger.commit_item`). At hand-out the attempt number is raised in the control database and written to the unit's `engine_items` row in the galaxy before the job starts. Every commit the unit makes through `ctx.commit()` updates that row `WHERE attempt = mine`, in the same transaction as the data; no row touched raises `Replaced` and the unit rolls back and exits. The last commit sets the done time with the last data. The message the job sends the supervisor is only a nudge; the supervisor reads the row.

**Checkpoints** (OPS.219). `ctx.commit()` also stores the number of pieces finished and what they returned. A retry or a resume starts after them (`ledger.resume_point`), for a command that `resumes`; one that cannot resume starts from its first piece and says so (`AGAIN`).

**Resume from rows** (`queue/resumed.py`, OPS.222 part B): a run that carries on from a stopped one carries over every stage the stopped run finished and, in each other stage, marks done the units the earlier rows say are done (their stored results copied) and queues only the rest. An update drains a run to `held` with the reason "stopped for update" (`held.drain_run`); the units stop at their next piece, and the end of the update carries the run on (`held.resume_after_update`) when the generator epoch is the one it was drained under.

**Attempt history** (`queue/attemptlog.py`, OPS.223): each attempt's end, error kind, error text and stack file; a setting runs any restarted job in debug mode with its own log file. The Queue page and the run's Jobs view read them (`queue/jobview.py`).

## 4. The supervisor

One process per host (`queue/supervisor.py`, `cli/supervisor.py`). Its loop takes a clock, a store, a process table and its duties as arguments, so the simulator drives the real loop with fakes. Each turn runs the duties in this order (`cli/supervisor.py`):

| Duty | What it does |
|---|---|
| `Reconciler` | Compares the run records with the rows and processes and corrects (R4: recompute, do not remember): a stage that finished and never queued its next opening gets it; a run whose starter is gone is carried on by the supervisor. Every correction is said in the supervisor's acts. |
| `RunStarts` | Starts runs: a process that wants a run leaves an ask (`queue/runline.py`) and follows the run the supervisor starts for it; a run that finds the machine taken waits in the line. |
| `RunWatches` | The watch on every live run's units (`unitwatch.UnitWatch`, `beatjudge`, `unitstuck.Stuck`) that used to live in the starter. |
| `WorkerStarts` | Starts and sheds RQ workers for the queues that still use them (until OPS.235). |
| `pressure.Controller` | The one memory controller (section 6). |
| the database signal | The same controller with the latency of the engine's own statement as its signal (`queue/dbsignal.py`, OPS.248). |
| `controls.Controls` | Pause, resume, Stop All, cancel, clear and the stuck-unit stop, as requests (section 9). |
| `lightjobs.LightJobs` | The site's light jobs (section 5). |
| `Launches` | The launcher: hands out the ledger's ready jobs (section 5). Holds its hand-outs while the queues are paused by hand. |

If the supervisor dies, the jobs it forked notice their pipe closed, finish the commit they are in and exit at the next safe point (`template.EXIT_SUPERVISOR_GONE`). The next supervisor reads the ledger, waits for processes that still exist and puts back what ran under processes that are gone (`Launcher.recover`). No run fails because a watcher died.

## 5. Jobs, templates and dispatch

**Templates** (`queue/template.py`, `queue/jobhost.py`). A template is a process that has imported the code and built one **warm set** (`jobhost.WARM_SETS`): `tables` (the bright-star sampling tables, 1.5 GB, built once instead of once per worker), `app` (the site's app, for the maintenance lane) and `plain` (nothing). It holds no connection and one thread, so it is safe to fork. It starts when the first job that wants it is asked for, and every job is a fresh fork. The supervisor and the template talk over one `SOCK_SEQPACKET` pair; each job has a pipe of its own to the supervisor for what it is doing and what it is asked (`Link`: `send`, `requests`, `stop_requested`). The supervisor also holds a `pidfd` on each job, so the end of a job is known even if the template died. What a job may run is a function the allow-list names (`queue/allowlist.py`).

**The launcher** (`queue/launcher.py`) is the supervisor's turn that hands out the ledger's jobs. Each turn: (1) settle every job process that ended (the unit writes its own end; a row still `running` after its process is gone is closed `lost`, with the signal named, and handed out once more, else `failed`); (2) read the `ready` rows (the oldest 500) as tickets; (3) `dispatch.dispatch` says which go now, and writes the reason of each that does not go to its row's `waited`; (4) mark each admitted row `running`, then fork it from its lane's template running `forkedjob.run`, which takes the attempt and runs the unit through `units.run_unit_on`.

**Dispatch** (`queue/dispatch.py`) is pure. The site's lanes (interactive and long) are served first against the places kept for them (they do not age against bulk work: ranked together a visitor waited 1,823 s at the 99th percentile in the simulator); every other lane is ranked by `scheduler.rank` (class weight, aging, fair share) and judged by `scheduler.judge`, with backfill. A reservation lasts until the job it was made for runs. Only the top 64 are judged each turn. A claim is the learned figure of the kind (`claims.py`, `memorygate.claim_mb`, `jobprofiles`) raised from the typed one; named resources are tokens (OPS.232: `mysql_write`, `alone`, one per warm set).

**Lanes.** Rows carry a `job_class`: `bulk` (a generation run's units and openings), `maintenance` (the warm-up and the checks), and the site's `interactive` and `long`. The maintenance lane has its own cap (`queue.maintenance_workers` while a generation run is live).

**Light jobs** (`queue/lightjobs.py`, `queue/api_jobs.py`, `queue/pagejobs.py`). An API job, a slow page's answer, a node prefetch and the one-off system page are not rows: they are re-derivable and must add no MySQL write to the site's fast path. Each is a record in Redis that the `LightJobs` duty hands out, oldest first through the same dispatch rule, as a process forked from the `plain` template running `lightjobs.run`. A job says how long it stays wanted and is dropped, if it has not started, when nobody asks after it; its time limit is learned from the kind's own history once it has enough (OPS.251). `api_jobs.submit`, `status`, `wait` and `cancel` keep the shape callers had; `call_and_wait` runs a function in the warm process and returns what it returns.

## 6. Pressure and the database

`queue/pressure.py` is the one controller that replaced the adaptive worker limit, the per-stage caps and the lend rule. Its signal is the kernel's memory pressure (PSI: `some` and `full`) where there is one, else available memory against what the gate keeps back. **Pressed** for `HOLD_AFTER_SECONDS`: nothing of the bulk or maintenance lanes is admitted (the site's lanes go on, and every waiting row says why). **Starved** for `YIELD_AFTER_SECONDS`: the youngest bulk job is asked to stop at its next piece (`Launcher.yield_youngest`, the pipe's `stop` request); its row goes back to `ready`, and the stop does not count against its retry (`ledger.STOPPED`). One at a time, `YIELD_COOLDOWN_SECONDS` apart. The kernel's own kill stays the last line, and a killed job is retried alone with the claim grown 1.5 times. `queue/dbsignal.py` gives the same controller a second look: how long the supervisor's own most frequent statement takes, learned against its own calm figure, so admission slows when the database is the limit.

## 7. What is still the old queue

No stage keeps two paths, but the cutover is not finished. Until these items merge, RQ still runs:

| Still on RQ | Until |
|---|---|
| **Runs with no ledger** (`generate system`, the update's math check): their units are RQ jobs (`units.enqueue_units`, `units.run_unit`), run by RQ workers that `WorkerStarts` starts | OPS.271 (a ledgerless run runs its units without RQ) |
| Non-unit RQ jobs the update guard lists (`updateguard.QueuedJobSource`), Stop All's queue emptying (`stopall._empty_queues`), `runcontrol.stop_units`, the RQ branch of `unitstuck` | OPS.271, then OPS.235 |
| The `rq` dependency, `workers.py`, `redisqueue` queue helpers, `cli/worker.py` | OPS.235 |

A run that is in the ledger never touches these: its jobs are forked, its stops are requests (section 9), and its failures are by kind.

**A Generate page job is a run of the ledger** (OPS.244 slice B, `queue/pagerun.py`). The run's kind is `page`; one row holds what the job is (its id, kind, title, who started it, the database, the form fields, the estimates) and one row each holds a step (`page-step`: label, command line, stages). Only the next step is `ready`; a step that ends well makes the one after it `ready`, a step that fails cancels the rest and holds the run with the reason. The launcher hands a step out like any other row (`launcher._job_of`), begins its attempt with the process it forks, and `forkedjob.run_page_step` runs `web/job_runner.run_step`, which starts the step's command line as a process of its own and follows it: the log (`output.log`, `records.log`, `progress.json`, `memory.json` stay in the job's folder), the stage marks, the memory peak, and the row (a Cancel that changes the row stops the step). The page reads the job from the rows (`pagerun.read`: starting, running, succeeded, failed, canceled, interrupted) and cancels with a request to the supervisor (`controls.stop_engine_run`), or ends the rows itself when no supervisor is alive. A page job is never retried by the launcher; a step cut off by a restart or a kill is the admin's to start again (the page's Retry). A supervisor with a page job on `supervisor.pages_key` does not end while the job runs. The job folder no longer holds `job.json`, `state.json`, `runner.pid` or a `cancel` file, and there is no private queue or burst worker.

## 8. Measured

Stage 0 (2026-11-13, a 4-core 16 GB machine, MySQL 8.0.46 with default durability; scripts in `/mnt/project-files/engine/stage0/`):

* **The ledger costs about 3 ms of MySQL time per job**: insert the plan of a 4,000-job stage 0.82 s; hand-out 1.2 ms (2.2 ms at the 99th percentile); the fenced done mark 1.4 ms more than the data alone; reading back 8,000 rows 21 ms. A replaced attempt committed 0 rows. RQ's path was 14.2 ms per forked job.
* **The installed MySQL account cannot read the Performance Schema**, so the engine does not depend on MySQL's own progress figures: a job reports its own statement (`ctx.sql`), and long work is cut into counted pieces.
* **Forked jobs share the tables**: building them takes 52 to 64 s and 1,540 MB once; a forked job starts at 1 MB private and copies what it reads (120 MB after 416 stars, 464 MB after 8,336). A template and 6 jobs held 2,422 MB where 6 warm workers held 9,240 MB for the tables alone.
* **The simulator** (200 seeds) found three things a first design would have shipped: counting only the count as evidence stopped working jobs in every seed (6.7 a run), cured by processor and byte evidence (none in 200 seeds); ranking the whole backlog with aging starved the site; and a reservation that ended after 300 s let a cut be passed for ever.

**The first one-off system page after the supervisor has been idle** (OPS.244 slice A, 2026-11-17 EST, the test run on a cloud machine): about 95 s from request to page, against the page's 120 s limit. That is a supervisor start plus the template's build; the next requests took seconds. It happens because the supervisor exits after `IDLE_SECONDS` (300 s) and the template with it. Measure it again after OPS.244 slice B and after any change to the idle exit or to what the `plain` template warms.

## 9. Controls

`queue/controls.py` (OPS.234). Pause, resume, Stop All, cancel and clear are **requests to the supervisor** (`controls.ask`: a JSON line on a Redis list, a reply list per request). With no supervisor alive none is started for a request: no job process of the ledger can then be running (a job whose supervisor is gone exits at its next safe point), so the caller does what it always did and `controls.end_rows` settles the ledger itself.

* `pause` and `resume` write and clear the flag the engine reads between two pieces (`queue/pause.py`); while it is set the launcher only settles ends and forks nothing new. A solo hold is not a pause of the hand-outs.
* `stop_jobs` (a canceled or cleared run) and `stop_all` cancel the rows still waiting or ready, ask each job process of the run to stop between two pieces through its pipe, and kill one still alive after `GRACE_SECONDS`; an attempt that was running is closed `lost` with the kind `canceled`. Stop All's reply comes when the last process is gone, so the process scan that follows it (`stopall`, the second line) has nothing of the engine's to find. The scan could not find a forked job in any case: its command line is the template's.
* `stop_unit` (the stuck-unit follower, `unitstuck`) kills the one job process that runs a stuck unit and closes its attempt `stuck` with the reason and the stack; the launcher hands the unit out again, and a second stuck verdict fails the run.
* The update guard drains a generation run (`held.drain_run`, a Redis hold the engine reads between pieces); `planetgen reset-jobs` has the supervisor stop every job process and ends the ledger's rows before it removes what RQ left in Redis.

## 10. Proof, and what is not yet proven

Built: the supervisor simulator (OPS.212: the real loop under a fake clock and store, with kill faults; a finished stage advances exactly once, a run ends or says why, nothing is queued twice), real-Redis and real-MySQL kill tests on the small galaxy (OPS.213: kill the starter, the supervisor, the last unit in a conclusion; the galaxy matches an undisturbed run), the guard tests of stage 1 (no command goes silent beyond a budget, no command sleeps or polls outside the listed primitives), the independence harness (OPS.241, declared reads and writes against a SQL tap) and the all-commands silence run (OPS.262).

Not yet built: the stage-3 proof (OPS.225: two attempts never both commit, a held run resumes to the same result, per-family kill tests), the stage-4 proof (OPS.236: a visitor's job starts within a second under any backlog, no job passed for ever, claims never over budget, kill tests for the template, a job, the supervisor and the kernel), the stage-5 proof (OPS.245: random graphs with crashes anywhere) and the full kill matrix with a soak (OPS.246). Until they merge, "the engine is proven" means the first two lists above and no more.

## 11. Where the code is

| Concern | Modules (`src/planetgen/queue/` unless noted) |
|---|---|
| Ledger, plan, fence | `ledger.py`, `planrows.py`, `attempts.py`, `attemptlog.py`, `needs.py`, `resumed.py`, `held.py`, `failures.py`, `jobview.py`; migration 0261 |
| The engine's loop and commands | `engine.py`, `commands.py` (`CommandSpec`: `pieces`, `do`, `resumes`), `longsql.py`, `unitbeat.py`, `beatjudge.py`, `unitwatch.py`, `unitstuck.py`, `jobbeat.py` |
| Supervisor and its duties | `supervisor.py`, `cli/supervisor.py`, `launcher.py`, `controls.py`, `pressure.py`, `dbsignal.py`, `lightjobs.py` |
| Templates and job processes | `template.py`, `jobhost.py`, `forkedjob.py`, `allowlist.py` |
| Admission and order | `dispatch.py`, `memorygate.py`, `scheduler.py`, `claims.py`, `jobprofiles.py`, `profileclaim.py`, `profiletime.py`, `usage.py` |
| Runs and their Redis record | `runrecord.py`, `runline.py`, `wrapper.py`, `units.py`, `runcontrol.py`, `stopall.py`, `pause.py`, `runview.py`, `runeta.py` |
| The site's jobs | `api_jobs.py`, `pagejobs.py`, `web/jobs.py`, `web/job_runner.py`, `web/system_page.py` |
| Still RQ | `redisqueue.py`, `workers.py`, `requires.py`, `cli/worker.py`, the RQ branches named in section 7 |
