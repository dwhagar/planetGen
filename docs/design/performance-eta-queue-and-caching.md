# Performance, progress estimates, the job queue and caching

Where generation spends its time, how a progress bar or banner should estimate
the time left for work whose speed varies, how the Redis and RQ queue behaves
(measured against the pinned `rq` 2.12.0), and how the page, tile and stamp
caches should be invalidated. It extends the queue audit
([work-queue-audit.md](work-queue-audit.md)), the resumable-run design
([sampling-backfill-and-resume.md](sampling-backfill-and-resume.md) section 7) and the
scheduling note ([ops-scheduling-and-rotation.md](ops-scheduling-and-rotation.md)) and
does not repeat them. The uploaded "Web UX and Job Management Guide.md" and "Web UX
Development Notes.md" are source documents and are not edited; section 7 says which of
their statements are done, replaced or still open.

Informs: PERF.31, PERF.32, PERF.33, PERF.34, PERF.20, PERF.18, UX.3, ADM.15, DB.15 (and the finished PERF.3, PERF.7, PERF.10, PERF.24, PERF.25, ADM.39 they build on)

Status: research, 2026-10-09; decisions marked "Boss" are his, everything else is a recommendation

The measured profile of the scatters and the sector fill is in [generation-performance-study.md](generation-performance-study.md).

Evidence tags: [C] computed or measured in this research (scripts under the research
scratchpad `perf/`, listed in Evidence notes); [S] read in the repository, a downloaded
package or a search result (URLs under Sources); [R] recalled and unconfirmed. The
research environment could read only search-result text, not papers or manuals.
Measurements come from one shared 4-core Linux sandbox (Python 3.12, MariaDB 10.11,
Redis from the distribution, warm file cache): they show orders of magnitude and
ratios, not the production server's speed.

## Decisions already taken

- **UX.3 (Boss, 2026-10-01):** "a warning to all users on the UI when a task is running in the background which is modifying the starmap ... with an ETA until it will be finished to the nearest hour rounding up." The banner reads the RQ job's published progress (ADM.22), not `progress.json`.
- **PERF.33 (Boss, issue #661):** "time remaining on all progress bars should be calculated from this performance metric averaged with the actual performance at the time"; a sub-task "expected to take longer than 15 seconds" gets its own bar.
- **PERF.32 (Boss, issues #661 and #750):** record rates "of every generation ... along with # of workers", show them on the stats page, an admin reset, and "Stats should be deleted on every new version."
- **PERF.31 (Boss, issues #761 and #750, 2026-10-09 04:43Z):** "a timed analysis in full debug from the planning of the galaxy to the galaxy being ready" and "a method to benchmark."
- **PERF.20 (Boss, 2026-10-01 22:17Z):** "Everything should go through the work queue so that we can implement short-term caching where possible"; an item "that needs planning."
- **ADM.15 (Boss, 2026-10-01 22:11Z):** "Ludicrous Speed" up to 95% of the CPU, "periodized so that the web interface still works ... and DB calls still work."
- **Queue (Boss):** Redis and RQ (2026-10-03); short single-object generation queued and waited on (approved 2026-10-07); wiki uploads on the queue (2026-10-08 00:14Z). API edits wait up to 8 s, else answer 202 with `GET /api/jobs/<id>` (PERF.24).
- **DB.15 (Boss, issue #727, 2026-10-09 01:18Z):** the migration bar estimates the time remaining of a single migration.

## Findings in short

1. **The bar's ETA is biased high at the start.** The repository's own `DecayingRate` (`queue/progress_rate.py`), fed a 4-worker run of equal units, reads 2.0 times the true time left at 5% done, 1.6 at 10% and 1.2 at 25% [C], because it starts from the first single completion. A ratio of two decayed sums has no such error.
2. **Weighting units by predicted cost cuts the error from about 50% to about 12%** over six simulated runs and the worst case from 181% to 26% [C]. The predictor already exists: PERF.3 sums `seconds_per_system(kind, density) * expected systems` before a run. The gain holds with a 100% per-unit prediction error (21% against 56%).
3. **Recorded speed must carry the worker count and must not be divided by it again.** One worker measured 0.04 s per system and three workers 0.17 s each (`stats.py`), so a one-worker figure divided by three is 4.3 times too optimistic [C, arithmetic]. The bright-star fill scales 1.87 times at 4 workers [S, report]; the generator alone scales 4.04 times on 4 cores [C].
4. **Every API job pays about 2.4 s and 195 MB before it does anything.** `api_jobs.submit` starts a fresh burst worker per job, and `execute` imports `planetgen.api.common` (2.3 s) before it runs the function [C]. Nothing caps the number of such workers; Boss's batch wiki uploads would start one 195 MB process each.
5. **RQ does not deduplicate by job id, `cancel()` does not stop a running job, and interval retries did not run on a burst worker** [C].
6. **The cache stamp holds one linear `COUNT(*)`** (0.13 s per million placed sectors, 2.7 s at 20 million [C]); it changes at every check during a fill, which empties the page cache for that database each time; and the planetGen version in its base flushes every tile at each release.
7. **In-memory generation is about 5 ms per system** (default configuration) against 16.6 ms per system in a filled sector (a different mix of systems) [C, S], so the database is the target of PERF.31 and PERF.34. Of the profilers, `py-spy` cost nothing measurable, `pyinstrument` doubled the run time and `cProfile` tripled it [C].

## 1. How progress and time left are produced today

| Mechanism | Code | Estimates | Weak point |
|---|---|---|---|
| Bar ETA | `queue/progress_rate.py` `DecayingRate`, one per rich task (`generation/run_common.py`) | units per second, exponential average over time, tau 60 s; ETA `remaining / rate`, never below `(remaining - 1) / rate` | counts units, not cost; start-up bias (2.1); tau fixed whatever the task length |
| Published progress | `queue/progress_file.py` writes `progress.json` (completed, total, rate, `eta_s`, `detail`) at most twice a second; `web/generate_page.py` counts `eta_s` down | the bar's value | one file per web job; the banner (UX.3) does not exist (no template mentions it) |
| Job tree ETA | `queue/work.py` `_roll_up`: `remaining * (mean seconds of done tasks) / workers`, children's ETAs added | the Queue page and `web/queue_page.py` `_remaining_seconds` (ADM.39) | mean of tasks done so far (order bias); a step with no finished task adds nothing, so a Generate job's later steps (plan, then fill) are left out; a task in flight counts as not started |
| Pre-run estimate | `generation/stats.py` `estimate`: per sector `seconds_per_system(kind, density) * expected systems`, summed, divided by `min(workers, sectors)` | the size and time shown before a run (PERF.3) | assumes linear scaling; stored seconds were measured at whatever worker count the earlier run had |
| Stored speed | `generation_stats` rows per `(kind, density half-decade)`, decay 0.05 per sector (PERF.10) | feeds the pre-run estimate | no worker count, no phase breakdown, not deleted on a new version |

`queue.work.timing_by_kind` computes mean, minimum and maximum seconds per kind and
has no caller [S, grep].

## 2. Estimating the time left

### 2.1 What is wrong with the present estimators

The test fed the real `DecayingRate` a 4-worker run of 400 units (3 s mean, lognormal
spread 0.45), one completion at a time like `advance()`, over 200 seeds
(`perf/dr_check.py`) [C]:

| Units done | Median ETA / true time left | P10 | P90 |
|---|---|---|---|
| 5% | 2.03 | 1.61 | 2.37 |
| 10% | 1.63 | 1.40 | 1.80 |
| 25% | 1.19 | 1.11 | 1.26 |
| 50% | 1.04 | 0.97 | 1.12 |
| 90% | 0.98 | 0.86 | 1.07 |

`DecayingRate.last` starts at construction, so the first completion, one task-length
later, gives a rate of one unit per task-length (a quarter of the four-worker rate), and
a 60 s time constant takes minutes to forget it. For a multi-hour run the error is gone
in minutes, but the UX.3 banner is shown from the first second and a Generate job's
steps are often short.

Two structural problems are larger:

- **Counting units when units differ.** A sector's cost scales with its density, and a ring fill crosses densities in a fixed order. When cost falls tenfold along the run, every count-based estimator is wrong by 70% to 180% (2.2).
- **Linear worker scaling.** `estimate` and `_roll_up` divide by the worker count. Measured: 19.7, 34.7 and 36.8 bright stars per second at 1, 2 and 4 workers [S, report]; the generator alone, without a database, scales 1.00, 2.23, 3.19, 4.04 times at 1 to 4 processes and stops at 4 cores (`perf/scale_bench.py`) [C]. The loss is the database and the serialised neighbour-lock section (17.6 of 46.6 ms per bright star), a property of the work and the machine that has to be measured, not modelled.

### 2.2 Estimators compared

Method [C, `perf/eta_sim.py`, 30 seeds per scenario, deterministic]: 400 units on 4
simulated workers, 3 s mean per unit, lognormal spread 0.45 in the real cost, and a
predicted cost known beforehand with lognormal error 0.25 (a per-density model).
Score: mean absolute percentage error of the ETA against the true time left at 10, 25,
50, 75 and 90% done. Scenarios: `iid`; `dense_first` (cost falls tenfold along the
run); `sparse_first` (cost rises); `two_classes` (20% of units cost 8 times more, mixed
at random); `load_spike` (machine 3 times slower for the middle quarter of the time);
`prior_wrong` (equal units, recorded speed 2.5 times too fast).

| Estimator | iid | dense_first | sparse_first | two_classes | load_spike | prior_wrong | Mean | Worst |
|---|---|---|---|---|---|---|---|---|
| Cumulative average | 6 | 181 | 69 | 15 | 21 | 6 | 49.5 | 181 |
| Window of the last 30 s (the shape of rich's estimate [R]) | 9 | 105 | 37 | 22 | 49 | 9 | 38.4 | 105 |
| Time-weighted EWMA of the count rate, tau 60 s (the repository) | 18 | 171 | 48 | 19 | 29 | 18 | 50.5 | 171 |
| tqdm-style EMA of increments, alpha 0.3 | 15 | 84 | 36 | 35 | 50 | 15 | 39.2 | 84 |
| Ratio of decayed sums, count units, tau 60 s | 7 | 132 | 50 | 17 | 29 | 7 | 40.1 | 132 |
| Cost-weighted, EWMA of the cost rate | 17 | 19 | 25 | 26 | 29 | 17 | 22.4 | 29 |
| **Cost-weighted ratio of decayed sums, tau 60 s** | 7 | 7 | 9 | 13 | 29 | 7 | **11.8** | 29 |
| Same, blended with the recorded rate (shrinkage k = 15) | 6 | 6 | 9 | 12 | 26 | 10 | **11.6** | 26 |
| Recorded rate only | 14 | 14 | 16 | 14 | 17 | 140 | 35.7 | 140 |

Cells are percent error (mean over the five checkpoints). Cumulative average,
least-squares slope and the other count-based variants behave like the rows shown.

- Every count-based method fails on `dense_first`, the shape of a real ring fill. No smoothing constant fixes an order effect; only knowing the remaining cost does.
- The ratio of decayed sums beats the EWMA of the rate by the start-up error of 2.1.
- tqdm's EMA reacts to every update (step-to-step jitter 10% to 13% against 1% to 3%), the wrong default when many workers finish together. Its `smoothing` runs from 0 (average speed) to 1 (current speed), default 0.3 [S].
- A recorded rate alone is good when right (14%) and ruinous when the machine or code changed (140%). Blending with weight `n / (n + 15)` on the live value bounded that case at 10%.
- Ninja does the same in spirit: the remaining time comes from each edge's duration in the previous build, edges without history count as unpredictable, a sliding window over the last N edges is used only for the displayed rate, and the code warns that old timings can mislead (a previous build served from a compiler cache) [S, search result on the status printer source].
- Prediction quality matters little (`perf/psig.py`, mean error over five scenarios): with per-unit prediction error 0.25, 0.6 and 1.0 the blend scored 12.0%, 14.1% and 19.3% and the cost-weighted ratio 13.0%, 15.4% and 21.3%, while the repository's EWMA stayed at 56.4%. The model needs the shape of the cost along the run, not accuracy per unit. Where there is no model (migration batches, backfill chunks) the estimator reduces to the count ratio at no cost.
- Window length matters for slow tasks (`perf/tau_test.py`). With 3 s tasks tau 60 s is best (6.1% error, jitter 2.0%; a 15 s window gives 8.2% and 4.5%). With 60 s tasks (a bright-star layer) a 60 s window holds under four completions: error 19.7% and jitter 10.7%, against 12.7% and 3.0% with a 300 s window (20 task-lengths across the 4 workers). On two cost classes at 60 s: 39.5% to 23.1%; on a load spike: 39.8% to 26.5%. Rule: `tau = max(60 s, 20 * mean task wall seconds / workers)`.

### 2.3 Recommended estimator

State per bar: decayed cost completed `N`, decayed time `D`, the predicted cost of each
unit (1 per unit without a model), the recorded rate `R0` in cost units per second for
this kind of work and worker count, and the number of finished units `n`.

```
on a completion at time t, with predicted cost c (c = 1 without a model):
    k   = exp(-(t - t_last) / tau)
    N   = N * k + c            # cost completed, decayed
    D   = D * k + (t - t_last) # time elapsed, decayed
    t_last = t
live = N / D                   # cost units per second, all workers together
w    = n / (n + 15)            # live weight
rate = w * live + (1 - w) * R0 # R0 absent: rate = live; neither: unknown
eta  = remaining_cost / rate
```

- `c` is the predicted cost, not measured seconds; the machine's speed enters through `D`, so the model must be right in shape only.
- `remaining_cost` is the predicted cost of the units not yet finished, the same density sum `stats.estimate` forms. A unit in flight counts as unfinished.
- Keep the repository's floor: the shown ETA does not reach 0 before the last unit completes.
- `R0` is `rate_cost_per_s` of PERF.32's row for the run's `(kind, workers)` (section 6). With no row, show the PERF.3 default-based estimate labelled "from defaults".
- A multi-step job adds the steps' remaining costs; a step not yet started contributes its PERF.3 `Estimate.seconds`, stored on the step node when the job is planned (today it contributes nothing).
- One function serves the terminal bar, the published record, the Queue page, the banner and DB.15's migration bar (units = batches of rows). It is pure and testable with a fake clock.

ETA over true time left for this estimator, 40 seeds [C, `perf/eta_band.py`]:

| Done | Normal runs (P10 / P50 / P90) | Load spike (median) | Recorded rate 2.5 times off (median) |
|---|---|---|---|
| 10% | 0.91 / 1.01 / 1.15 | 0.85 | 1.24 |
| 25% | 0.90 / 1.00 / 1.10 | 0.81 | 1.08 |
| 50% | 0.88 / 0.99 / 1.08 | 1.80 (during the spike) | 1.04 |
| 75% | 0.86 / 0.97 / 1.10 | 1.09 | 1.01 |
| 90% | 0.73 / 0.95 / 1.11 | 1.01 | 0.99 |

For normal runs the middle 80% of estimates lie within 0.88 to 1.15 times the truth. A
wider range is honest only when the machine's speed changes.

### 2.4 What to show

- **Rounding up to the hour (Boss, UX.3).** A naive `ceil(eta / 1 h)` on a 10 hour run rose 9 to 18 times per run and showed fewer hours than remained for 15% to 25% of the time [C, `perf/hour_display.py`]. Showing `ceil(1.15 * eta)`, letting the value fall freely and raising it only when `0.88 * eta` exceeds it cut the rises to 0 to 0.1 per run and the under-promise share to 0.2% to 9.0% (worst on the two-cost-class run). Re-fit the factors from PERF.32's recorded errors.
- **Ranges.** Admin pages: "about 3 to 4 hours" from the same factors. The banner: one rounded hour.
- **Rate unknown.** Before the first completion show the PERF.3 estimate labelled as such ("about 4 h, from earlier runs"), or an indeterminate `<progress>` (no `value`) with "estimating" when no stored speeds exist. Hold the live estimate until `n >= 5` and 20 s have passed. Drop `-:--:--` from web pages.
- **When a bar is worth drawing.** Nielsen's limits (0.1, 1, 10 s) put the need for a percent-done indicator at 10 seconds [S]; Myers (CHI 1985) found people prefer one [S]. Boss's 15 s sub-task rule fits. Bar animations that change perceived duration by about 10% (Harrison, Yeo and Hudson, CHI 2010 [S, press reports]) are not needed; the template's plain `<progress>` is enough.
- **Stalls.** If no unit finishes for more than `max(60 s, 3 * the longest task so far)`, show "waiting, last progress N min ago" and stop counting the ETA down.

### 2.5 One published record

Write one record per run where every web process can read it cheaply: a Redis hash
`planetgen:progress:<run id>` (a `GET` took 120 us over a kept connection [C]) with
`completed`, `total`, `cost_done`, `cost_total`, `rate`, `eta_lo`, `eta_hi`, `phase`,
`detail` and `updated_at`, expiring a day after the run. The banner's endpoint (public,
no secrets) reads it through a 5 s per-process cache; the page polls every 30 s rather
than holding a stream open (the Guide's starvation warning, section 7). The job tree in
the control database stays the history the Queue page reads and stops computing its
own ETA. `progress.json` stays for web-started runs until the jobs folder goes.

UX.3 asks whether command-line runs should appear. They take no lock today. Let
`planetgen galaxy|plan|reset` always publish the record and refresh a Redis key
`planetgen:active:<db>` (60 s expiry), which the banner reads. The file lock stays for
mutual exclusion.

## 3. The job queue

### 3.1 What the queue is today

RQ here is a process launcher more than a shared queue [S, code]:

| Use | Queue | Worker | Admission |
|---|---|---|---|
| Web job (Generate page) | one per job, `web/jobs.py` | one detached burst worker | the `active` file lock admits one job |
| A run's tasks | one per run, `RQExecutor` | up to `worker_count()` burst workers started by the run | the `work_lease` row admits one run's pool at a time |
| API job (`api_jobs.submit`) | one per job, `planetgen-api-<id>` | one detached burst worker | none |

With a queue per job, RQ's ordering features are unused and nothing bounds the number
of API workers. A long neighbourhood request, a wiki batch and a Generate run can all
run at once. Their processes are at lowered priority, but the database they all write
to is not.

### 3.2 Measured cost of an API job [C]

| Item | Measured |
|---|---|
| `submit` + worker start + a no-op function + result | median 2.4 to 2.5 s (5 runs per case, range 2.2 to 3.5 s) |
| The same work in a warm in-process worker | 0.03 s, plus 0.2 s of imports |
| `python -m planetgen.cli.worker --burst` alone, including its job | 0.3 s |
| `import planetgen.api.common` | 2.3 s (`planetgen.db.store` 2.3 s of it; `nltk` 0.75 s, `scipy.stats` 0.5 s, name word lists 1.0 s) |
| Resident memory: bare interpreter / worker module / + `api.common` / + `run_galaxy` | 12 / 33 / 195 / 198 MB |
| `api_jobs.status()` | 2.9 ms per call, one new TCP connection each |

The 2.3 s is `execute` importing `ApiError` at its top, so a no-op job costs as much
wall time as 400 single-system generations, and the 8 s wait of `run_queued` spends
30% of its budget starting a process. The API daemon has 5 threads (`WSGIDaemonProcess
planetgen-api processes=1 threads=5`, `examples/apache/planetgen.conf.example`), and
each waiting edit holds one for the whole spawn, polling Redis about 19 times a second
on a fresh connection. The server has been OOM-killed before
([server-checklist.md](../server-checklist.md)); an uncapped 195 MB process per queued
job is the same failure waiting for a batch.

Recommendations, by value:

1. Do not import the API layer to catch a refusal: define a small `Refused(status)` base in `planetgen.queue` that `ApiError` subclasses, or look the class up lazily in the `except` clause. Expect about 0.4 s fixed cost (an estimate from the 0.3 s bare worker; not measured after the change).
2. Import `nltk` and `scipy.stats` lazily where they are pulled in at module level. This also speeds every CLI start (`planetgen system` took 2.5 s warm) and the tests.
3. Start an API worker only when fewer than `worker_count()` planetGen workers are alive (`rq.Worker.all`), and put API jobs on shared queues (3.4) so one worker serves several jobs in turn.
4. Reuse one connection in `status`/`wait`; poll at 100 ms, backing off to 500 ms after 2 s.

### 3.3 RQ 2.12.0 behaviour that matters (`perf/rq_sem.py`) [C]

| Question | Result |
|---|---|
| Two `enqueue` calls with the same `job_id` while queued | Both accepted: the queue held 2 entries for 1 job. No deduplication. |
| `enqueue` with the id of a finished job | Accepted; the old result is overwritten and the job queued again. |
| Worker on `[hi, lo]` | Strict order: all of `hi`, then `lo`. A steady stream of `hi` starves `lo`. |
| `RoundRobinWorker` on `[a(4), b(2)]` | `a0 b0 a1 b1 a2 a3`. `RandomWorker` also exists. |
| `at_front=True` | Goes to the head of its queue. |
| `job.cancel()` on a queued job | Removed; status `canceled`. |
| `job.cancel()` on a started job | Status `canceled`, **but the job kept running** (worker alive 1 s later). |
| `send_stop_job_command` on a started job | Work horse killed in 0.07 s; status `stopped`. |
| `Retry(max=2, interval=[0, 1])` after a failure on a burst worker | Status `scheduled`, not re-run. Likely cause, inferred and not isolated: burst workers run `with_scheduler=False`. |
| `job.meta` with `save_meta()` | Visible across processes; 200 writes added about 0.14 s (0.7 ms each). |

Consequences:

- **Deduplicating identical work (PERF.20)** needs `SET planetgen:inflight:<digest> <job id> NX EX <ttl>` before enqueuing; return the stored id if the key exists. A deterministic `job_id` is not enough.
- **Cancel.** The web job cancels by file because its runner has child processes ([work-queue-audit.md](work-queue-audit.md)). For tasks inside a run, `cancel()` for queued and `send_stop_job_command` for started ones are both safe: a sector save is one transaction and a backfill chunk commits levels and stars together ([sampling-backfill-and-resume.md](sampling-backfill-and-resume.md) section 7.1), so a killed horse leaves nothing half-written and PERF.30 picks the unit up again. `RQExecutor.shutdown` only waits for running tasks, so a cancel of a 60 s layer task waits up to 60 s; "cancel now" would end it in under a second.
- **Retries.** Keep interval-less `Retry(max=1)` for dead workers. A task's own exception already comes back as a value (`work.py`, `("error", exc, 0)`), so only a dead worker is retried, which is right. Do not add intervals without a worker running the scheduler.

### 3.4 Fairness and priorities

The aim is not fairness between users (only admins start bulk work) but that an edit
or a short estimate is not stuck behind hours of bulk fill, and that bulk work leaves
the site and database usable (PERF.34, ADM.15).

| Queue | Work | Workers |
|---|---|---|
| `planetgen-interactive` | single-object generation, estimates, one wiki page, settle of a few sectors, neighbourhood requests under a size limit | 1 reserved worker serving only this queue |
| `planetgen-bulk` | a run's sector, scatter and backfill tasks, Generate page runs | up to `worker_count()`, serving `[interactive, bulk]` in that order |

- The reserved worker bounds an edit's wait to its own start-up; without it an edit waits for the current bulk task (a bright-star layer takes about 60 s or more).
- Bulk workers take `[interactive, bulk]`, so interactive jobs go first; strict order cannot starve bulk work because interactive jobs are rare.
- The `work_lease` and the web `active` lock keep their job: one bulk run at a time, a second waits.
- **ADM.15.** The protections Boss asks for are the right ones, with one gap: lowered priority protects the web server's CPU but not the database, which is not niced and is the real bottleneck (5.1). Keep one connection per worker, cap total worker connections below the server's limit, set Ludicrous Speed to 95% of the cores while keeping the reserved interactive worker, and re-run PERF.34's page-time test under it.

## 4. Caching

### 4.1 What exists

| Cache | Where | Invalidation | Notes |
|---|---|---|---|
| Page cache (PERF.2, PERF.25) | `web/lib/pagecache.py`, in memory per web process, `cachetools.TTLCache` | any API write in this process; the content stamp every 15 s (a new stamp drops that database's entries); age limit 300 s | not shared between processes |
| Tile and stage cache | `web/lib/tilecache.py`, disk `<dir>/<db>/<generation>/` | the stamp from `GET /api/galaxy/changes` every 60 s; only changed tiles deleted; `full` starts a new generation; a busy fill keeps the old generation up to 600 s; 200 MB budget | the browser keeps its own copy by generation |
| Opening-view warm-up | `web/warmup.py`, `cli/warm_map.py` | run by `update.sh` | only the opening tiles and stage |

### 4.2 The stamp

`galaxy_content_state` (`db/query.py`) is `base` (hash of the galaxy shape, bright-star
seed, **`__version__`** and naming key), `COUNT(*)` of placed sectors, and four maxima
(sector id, `sectors.modified_at`, system id, bright-star id).

- **Cost.** The maxima are index lookups (8 ms). `COUNT(*) ... WHERE center_x_pc IS NOT NULL` scans `idx_sectors_center` and grows linearly [C, MariaDB 10.11, `perf/stamp_bench.sh`]:

  | Placed sectors | `COUNT(*)` warm | Maxima |
  |---|---|---|
  | 1 million | 0.13 s | 8 ms |
  | 5 million | 0.60 to 0.70 s | 7 ms |
  | 20 million | 2.7 to 2.8 s | 8 ms |

  The check runs per web process per database every 15 to 60 s. The count exists to
  detect deletion. Replace it with a counter only deletion changes (a `galaxy_epoch`
  column incremented by the delete routines), leaving the four maxima.
- **The version in the base.** Each release changes `__version__`, so every tile and cached page is dropped on update. Releases are frequent: 26 to 123 a day on recent days [C, `CHANGELOG.md` headings per date], so a development server never keeps a cache; a production server loses it at each `update.sh`, partly covered by the warm-up. A `TILE_FORMAT` constant bumped only when the tile builder or page serialisation changes would keep the cache across unrelated releases. The risk is forgetting the bump; a golden test removes it (build one tile and one page body from a fixture galaxy, hash the sorted key structure, compare with a constant, and fail with "bump `TILE_FORMAT`"). Do this once PERF.31 shows what a cold rebuild costs; if a full rebuild takes under a minute, leave it.
- **During a fill the stamp changes at every check**, so the page cache drops that database's entries every 15 s and serves almost nothing when the database is busiest. The tile cache already has `busy` handling (PERF.34); the page cache needs it too: while `busy` is reported, serve entries up to `max_age_seconds` and invalidate only on an edit in this process or by age.
- **Several web processes** learn of each other's writes only through the stamp. A cheaper shared signal: writers (the API, generation tasks after each chunk commit) `INCR planetgen:epoch:<db>` in Redis, and each process compares it with its own copy (120 us [C]) before using its cache. Keep the content stamp as the fallback when Redis is down.

### 4.3 Where a cached body should live [C, `perf/cache_cost.py`], microseconds

| Body | In-process `TTLCache` hit | Redis `GET` | File read | `json.loads` |
|---|---|---|---|---|
| 1 KB | 0.4 | 120 | 11 | 4 |
| 100 KB | 0.6 | 663 | 57 | 184 |
| 1 MB | 0.6 | 3,283 | 405 | 1,890 |

Redis suits small shared facts (epoch counters, the progress record, in-flight keys, a
job's result summary) and not page or tile bodies: a megabyte costs 3.3 ms against 0.4
ms from a file. Keep bodies in the per-process cache and the tile directory, and use
Redis only for the epoch that says whether they are still good.

### 4.4 PERF.20: what the queue can and cannot cache

- **Most queued work must not be cached.** Regenerate, re-rolling a body, neighbourhood and sector fills create new random content or change the galaxy (`regenerate_sector` draws a fresh seed). The job record, kept a day for polling, is all there is to keep.
- **What is repeatable and costly:** the Generate page's estimate (`--estimate-only`, up to 120 s), tile and stage builds, the stats page and search aggregates, wiki renders. Key them by `(kind, arguments digest, epoch, speed-stats version)`.
- **Routing reads through RQ is slower than computing them.** A queued job costs 2.4 s before it starts (3.2); a cache read costs microseconds. Cache in front of the work, in the web or API process, and use the queue only on a miss for work long enough to need it.
- **Coalescing concurrent identical requests is the real win.** `fetch_tiles` has no single-flight: when the cache is cold or just invalidated, every visitor runs the same slow query, which fits the PERF.34 symptom. Take `SET planetgen:build:<key> NX EX 30` before building; the loser polls for a few seconds, then builds itself.
- **Expiry stampedes.** For entries that expire by age and are costly (the stats page), recompute shortly before expiry with a probability that rises toward it (XFetch, Vattani, Chierichetti and Lowenstein, PVLDB 8(8), 2015 [S]; its beta defaults to 1 [S]). The epoch plus single-flight covers the rest.

Done for PERF.20: this plan, with build work filed as separate items in the handoff.

## 5. Where generation spends its time

### 5.1 Measured

One star system in memory (`StarSystem(system_config=SystemConfig())`, default
configuration, 600 systems after warm-up) [C, `perf/prof_bench.py`]:

| Profiler | ms per system | Overhead |
|---|---|---|
| none | 4.1 to 5.4 | |
| `py-spy record -r 100` (out of process) | 5.1 | not measurable here |
| `pyinstrument` 5.1.3, 1 ms interval | 12.1 | about +124% |
| `cProfile` | 14.6 | about +190% to +260% |

Under `cProfile`, self time by file: `physics/position.py` 22%, standard library 22%,
`util/checks.py` (the `finite_domain` decorator) 17%, `physics/planets.py` 14%,
`util/draw.py` 7%. `finite_domain` ran 500 times per system; speeding its leaf
generator gained about 3%, inside the noise (`perf/checks_bench.py`), so it is not a
lever. `pyinstrument` puts moon generation (`generate_moons`, a `Planet` per moon) at
about half the run and `SpatialPosition3D` construction at about 12%. These are
inclusive sampled times: enough to point at moons and position objects, not to rank small
functions.

In the database path the generator is the minority: a bright-star sector costs 46.6 ms,
6.6 of it `generate_sector`; a normal sector with 205 systems costs 3.4 s, about 16.6
ms per system against 5 ms of generation [S, bright-star timing report]. The report's
three levers (the serialised neighbour-lock steps, a second uid pass, one registry
upsert per sector) are the PERF.31 candidates that matter and are not repeated here.
Because the generator alone scales to all cores, extra workers buy little on a
database-bound run (1.87 times at 4); the cheap gain is fewer round trips per sector.

### 5.2 Tools for Python 3.9 and later

| Tool | Version | Python | Use here |
|---|---|---|---|
| `cProfile` | standard library | all | exact call counts; heavy on small functions; per process, so wrap a task |
| `pyinstrument` | 5.1.3, 2026-07-29, BSD, wheels cp39 to cp314 [S, PyPI] | 3.8+ | readable call tree for one command; signal mode samples the main thread only, so profile a task function inside the worker |
| `py-spy` | 0.4.2, 2026-04-24, MIT, `py2.py3-none` wheels for Linux, macOS, Windows x64 [S, PyPI] | any | launches or attaches from outside; `--subprocesses` follows forked RQ work horses; worked here in launch mode [C]; may need elevated rights on macOS and hardened Linux [R] |
| `memray` 1.20.0 | PyPI | 3.9+ | only if memory (195 MB per worker) becomes the question |
| `sys.monitoring` | standard library | 3.12+, not 3.9 | cannot be the project's method while 3.9 is supported |

None is a runtime dependency. `planetgen benchmark` itself uses only wall-clock timers.

### 5.3 PERF.31: the benchmark

`planetgen benchmark [--sectors N] [--profile]` plans a small fixed-seed galaxy in a
scratch database and runs the real phases (plan, scatter per layer, sector fills,
backfill, population, paths, checks), recording per phase and sub-phase:

- wall seconds, from a `perf.phase("name")` context manager (always on, microseconds per call), summed per run;
- units done, so rates per second fall out and go into PERF.32 with the worker count;
- the database's share, from `Com_*` and `Innodb_rows_*` counters (or `performance_schema` statement digests) read before and after;
- with `--profile`, a `pyinstrument` profile of one task per kind.

The sub-phases of a sector fill come from the timing report's table. Run it at 1, 2 and
4 workers and two densities so the scaling and cost-model inputs of sections 2 and 6
come from one run, and once with a page-request thread running, which is the test
PERF.34 describes.

### 5.4 PERF.34: suspects, best supported first

Status: PERF.34 shipped in PR #811 after this section was written: the Galaxy Map keeps its tile cache while a fill changes the database, and serves the cache when the freshness check times out. The suspects below were ranked before that change; re-measure the page-time test on current main before acting on any of them.

1. A cold or just-invalidated cache with no single-flight (4.4).
2. The page cache emptying at every stamp check during a fill (4.2).
3. The stamp's linear `COUNT(*)` once millions of sectors are placed (2.7 s at 20 million).
4. Five API threads held by `run_queued` waits of 2.4 s or more (3.2).
5. Database contention: worker commits and the serialised neighbour-lock section. Not measured under page load here.
6. Connection exhaustion. Not measured; the cap in 3.4 prevents it.

Items 1 to 4 can be fixed without waiting for the benchmark.

## 6. PERF.32: what to record

Extend `generation_stats` rather than add a table:

| Column | Meaning |
|---|---|
| `version_key` | `galaxy/version_key.py`'s 22 hex digits (release, Python, OS, architecture) |
| `kind` | `sector`, `scatter`, `backfill`, `phenomena`, `population`, `migration`, and PERF.31's phase names |
| `workers` | worker count of the run |
| `density_bucket` | the existing half-decade bucket, 0 for kinds without density |
| `n`, `sum_wall`, `sum_units`, `sum_cost` | counts and sums, so decayed means and their ratio can be formed |
| `rate_cost_per_s` | decayed throughput of the whole pool at this worker count |
| `updated_at` | |

- Record **pool throughput** at `workers`, not only per-task seconds, and let the estimate read `(kind, workers)` and interpolate between neighbouring worker counts instead of dividing by the count. With no row for the count, scale the nearest one by the measured curve once two rows exist and show a range until then.
- `stats.estimate` and the live estimator (2.3) read the same row; `R0` is its `rate_cost_per_s`.
- **Deleted on every new version (Boss).** Delete rows with another `version_key` on the first read or write that sees one, in one transaction, and log how many went. The key includes Python, OS and architecture, so a Python upgrade also resets, which is right. The cost: the first run after each update has no recorded rate, when PERF.33 wants one; on a development server with 26 to 123 releases a day [C] the table never fills. See the first open question in the handoff.
- The stats page shows the table; an admin reset deletes it; benchmark rows carry a `bench:` kind prefix so they never feed a live ETA. Stored speeds keep surviving a galaxy reset ("they describe the machine, not the galaxy", `stats.py`); only a version change and the admin reset delete.

### Built (PERF.32)

Control schema v12 (`generation_stats`, `control_schema.sql`). As built:

- The key is `(kind, workers, bucket)`; `version_key` is a column, not part of the key, because every row of another key is deleted by the first write of a run (`GenerationStats.flush`, once, logged at debug level) and ignored on read. A table from before v12 is dropped and recreated (`store._drop_old_generation_stats`), as the version change would have deleted its rows anyway.
- Kinds recorded: `sector` (a sector's fill), `scatter` (a bright-star layer) and `phenomena` (a phenomenon-scatter layer), each with the run's worker count (`run_common._worker_count`). The sums and `rate_cost_per_s` columns of the sketch above are not added: the rows still hold the decayed per-task and per-system seconds, and PERF.33's estimator is where a pool-throughput column belongs.
- `seconds_per_system(kind, density, workers)` reads the rows for that worker count; with none, it blends the two neighbouring counts by distance, or takes the nearest when only one side exists. `estimate` passes the run's count, so a four-worker run no longer borrows one-worker times.
- `BENCH_PREFIX` (`bench:`) is reserved for PERF.31's benchmark rows; no live estimate asks for such a kind.
- The Stats page lists the rows with their worker counts and has a Reset stats button (`POST /api/admin/generation-stats/reset`, audited as `generation_stats.reset`).
- The first open question stands as Boss defaulted it: delete on every new version, no flagged prior.

## 7. The uploaded guides, checked against the code

| Statement | Status |
|---|---|
| Job Management Guide: broker-less `ProcessPoolExecutor` with SQLite state, `future.cancel()` | Replaced by Redis and RQ (Boss, PERF.24); marked superseded in [api-design-standards.md](api-design-standards.md). Its point that cancelling must terminate the worker is right: RQ's `cancel()` alone does not stop a running job (3.3). |
| Guide: SSE route `/api/jobs/<id>/stream` fed by an in-memory queue per job | Built differently (a page route reading the log by byte offset, ADM.22; see api-design-standards.md); an in-memory queue cannot reach a separate worker process. |
| Guide: WSGI worker starvation during log streaming | Still applies to any held connection, including a banner; poll every 30 s instead (2.5). |
| Development Notes section 4: native `<progress>` with ARIA; delete bespoke DOM updates in `generatefolds.js` and `progressRate.py` | Native `<progress>` is in `partials/job_status.html`. `progressRate.py` is now `queue/progress_rate.py`, the rate source this document replaces. `generatefolds.js` and `generatejobs.js` still exist. |
| Development Notes section 5: RQ, a broker to deploy and monitor | Chosen; OPS.21 and OPS.27 cover the broker. |
| "Faster background execution reduces ... job-status timeouts" | For single-object work the queue adds 2.4 s until the start-up cost in 3.2 is removed. |

## Evidence notes

Computed [C] in the research scratchpad `perf/` (copied to `/mnt/project-files/research/scripts/perf/`; re-runnable; the ETA simulator is
deterministic by seed): `eta_sim.py`, `eta_band.py`, `dr_check.py`, `tau_test.py`,
`psig.py`, `hour_display.py`; `spawn_bench.py`, `phases*.py`, `redis_cost.py`,
`cache_cost.py`, `rq_sem.py` (with `rq_sem_out.txt`); `stamp_bench.sh` (MariaDB
10.11.14, 2 GB buffer pool, synthetic `sectors` table with the production index);
`prof_bench.py`, `prof_share.py`, `checks_bench.py`, `scale_bench.py`. The simulated
runs assume a cost model (prediction noise 0.25, real-cost noise 0.45); measure the
real prediction error from PERF.32's data and re-fit the band factors of 2.4. The 2.3 s
import and 195 MB figures are warm-cache values; a cold server pays more, and Windows
`SpawnWorker` re-imports in a fresh interpreter, so its fixed cost is probably higher.
Release counts per day come from the `CHANGELOG.md` headings.

Recalled [R], to verify when access is allowed: that rich's speed estimate is a 30 s
window of samples (the table row is that shape, not rich's code); py-spy's ptrace
requirements on macOS and hardened Linux; the roughly 10% figure of Harrison et al.
(press reports only); Myers' preference finding (abstract only).

Not measured: the production server; page times under a fill with the cache changes;
the effect of the import fix; Windows and macOS start-up; connection exhaustion under
a fill.

## Sources

- `rq` 2.12.0 wheel (PyPI): `rq/worker/base.py`, `rq/command.py`, `rq/job.py`, `RoundRobinWorker`, `RandomWorker`, `Retry`; behaviour tested locally (3.3).
- PyPI JSON pages for versions, dates, licences, Python support and wheel tags: https://pypi.org/project/py-spy/, https://pypi.org/project/pyinstrument/, https://pypi.org/project/memray/, https://pypi.org/project/rq/, https://pypi.org/project/cachetools/
- tqdm `smoothing`: https://pypi.org/project/tqdm/4.28.1 (PyPI copy of the docstring)
- Ninja status printer: https://fuchsia.googlesource.com/third_party/ninja/+/refs/heads/main/src/status_printer.cc and `.../status_printer.h` (search results only)
- RQ 1.8.0 release notes (`RoundRobinWorker`, `RandomWorker`): https://newreleases.io/project/pypi/rq/release/1.8.0
- Nielsen's response-time limits: https://www.nngroup.com/videos/3-response-time-limits-interaction-design/ and https://www.uxtigers.com/post/progress-indicators
- Myers, CHI 1985: https://www.cs.cmu.edu/~bam/papers/percentdoneCHI85.pdf (in search results; not opened)
- Harrison, Yeo and Hudson, CHI 2010: https://www.figlab.com/research/2010/faster-progress-bars
- Vattani, Chierichetti and Lowenstein, PVLDB 8(8), 2015: https://iris.uniroma1.it/handle/11573/877847 and https://docs.rs/xfetch/
- py-spy and pyinstrument overviews: https://pypi.org/project/py-spy and https://pypi.org/project/pyinstrument/0.13.2/ (the project's own overhead benchmark: sampling cheapest, `cProfile` heavier; percentages differ between versions)
- Project files: the design documents linked above, `docs/server-checklist.md`, `/mnt/project-files/bright-star-timing/report.md`, `/mnt/project-files/galaxy-oom/apache-recovery.md`, `/mnt/project-files/site-load-fix/server-checklist.md` (the original of the repository checklist).
