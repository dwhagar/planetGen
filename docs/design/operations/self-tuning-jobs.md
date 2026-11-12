# Self-tuning RAM, timing, worker counts and splits (DOC.120)

Status: design note. OPS.155 (the recording, step 1) is built: `queue/jobprofiles.py`, `queue/usage.py`, control schema v22. The readers are not. Asked by Boss, 2026-11-12, 2:50 PM EST: can the system keep stored values for RAM, process times, worker counts and splits, and work them out itself from the machine's memory and performance? Boss added at 2:51 PM EST: not just Stage 17, every job queue. At 2:52 PM EST he said to start by recording memory per unit and per item together with time per unit and per item, so recording is the first piece and ships alone (see "Record first").

**Answer: yes, and about half of it already exists.** The gate already learns a peak per job kind, the ETA code already stores times and rates per stage, and the adaptive limit already moves the bulk worker count with memory. What is missing is that the learning is thin (one number, kept only in Redis, never used for the cut) and that the starting figures are typed constants. This note joins the pieces into one stored table per job kind and says what each consumer reads from it.

## What exists today (reuse, REUSE RULE)

| Piece | Where | What it already does |
|---|---|---|
| Memory claim and learned peak | `queue/memorygate.py` (`learn`, `claim_mb`, `grow`, `forget`), `queue/claims.py` (`claim_for`) | Every job of every queue and lane takes a claim. The claim is the typed table figure raised to the largest peak seen for the kind, plus 10 percent. A job killed for memory retries once at 1.5 times. Kept in the Redis hash `planetgen:gate:learned` as one number per kind (`unit:<command>`, `cut:<command>`, `api:<ref>`, `fn:<name>`). |
| Measured peak per job | `queue/units.py` (`peak_resident_mb`, `resident_mb`), the worker heartbeat (`memorygate.refresh`) | Pss of the worker and its horse, shared pages counted once. This is the sample the table would store. |
| Stage times and cost per worker second | `queue/stagetimes.py` | Running average per stage key, newest run weighs half. Redis only. |
| Items per second per command | `queue/unitrates.py`, `queue/progress_rate.py` (`blended_rate`) | Running average per command, folded in when a run ends `done`. Redis only. |
| Longest quiet gap per command | `queue/commandgaps.py` | Feeds the heartbeat judge. |
| Worker count follows memory | `queue/adaptive.py` (OPS.148, PR 6147) | Lowers the bulk limit after a 30 s memory wait or under pressure, raises it after 120 s of quiet when one more claim fits. Short-lived state, 1 hour, list of recent changes. |
| Fair share between lanes | `queue/scheduler.py` (OPS.130) | Charges cores times seconds per kind. |
| Cut by weight | `queue/units.py` (`static_chunks`, `coarsen`, `MAX_UNITS`, `STATIC_UNITS`) | The cut depends on the rows alone, never on the worker count. |
| Per-stage unit ceiling | `util/settings.py` `queue.unit_ceiling_seconds` (300 s) and `CommandSpec.unit_seconds` (23 commands declare one) | A test holds each command's declared longest unit to the ceiling. |

## What is hardcoded today

1. **Memory per worker.** `workers.STAGE_WORKER_MEMORY_MB` (the 2,100 / 700 / 3,050 MB classes, measured once on one galaxy), `claims.UNIT_BY_NATURE_MB`, `claims.API_CLAIMS`, `claims.FUNCTION_CLAIMS`, `claims.CUT_MB` (1,500), `claims.SUPERVISOR_MB` (700), `claims.DEFAULT_MB` (1,000). The learned peak can only raise a claim above these, never lower it, so a figure that was too high on this machine stays too high.
2. **Worker counts.** `workers.STAGE_MEMORY_SHARE` (0.5), `HOST_MEMORY_SHARE` (0.6), `WORKER_MEMORY_MB` (300), four per core, `cores + 4` floor, the interactive and long lane caps in settings, and the adaptive edges (30 s, 20 s, 120 s, 0.75, 1.25; OPS.151 is already open to tune them by hand).
3. **Gate edges.** `RESERVE_MIN_MB` 1,024 and `RESERVE_SHARE` 0.08, `CPU_SHARE`, `LOAD_LIMIT`, `MAX_WAIT_SECONDS` 300, `GROW_FACTOR` 1.5, `LEARN_MARGIN` 1.10, `KEEP_MB_PER_PLACE` 700.
4. **Split sizes.** `STATIC_UNITS` 512, `MAX_UNITS` 4,000, `MAX_SHELLS_PER_CHUNK` 64, `BACKFILL_CHUNK_SECTORS` 200 and `BACKFILL_UNIT_CELLS` 10,000, `LINK_UNIT_SECTORS` 500, `SETTLE_UNIT_SECTORS` 100, `MAP_LUMINOSITY_BATCH_IDS` 200,000, `HABITABILITY_BATCH_SYSTEMS` 500, `uncharted_index.UNIT_CELLS` 300,000 (OPS.153), the sky and bright-star batches, `web/warmup.py` `NODE_CHUNK` 64 and `OBJECT_CHUNK` 250, `tileprefetch.CHUNK` 16, plus the typed `unit_seconds` per command. Each was picked for one galaxy size on one machine. The Stage 17 index at about 10k rows per second in one process is the example: the split was fixed by hand after a 40 hour run.
5. **Times.** The 300 s unit ceiling and `roughtime` (the "not timed yet" guess).

## The scheme

### One stored profile per job kind

A table `job_profiles` in the **control database** (not Redis, so it survives a Redis flush and a restart, and a wiped galaxy database does not lose it). Key: job kind (the same `unit:<command>`, `cut:<command>`, `api:<ref>`, `fn:<name>` strings the gate uses) and a **size bucket** (rows or cells in the unit, in powers of two). Per key it keeps:

* peak MB: the median and the 95th percentile of the last N samples (N about 30), not only the maximum, plus the count and the largest ever seen;
* seconds per unit and items per second, same window;
* how it ended: done, killed for memory, timed out, stuck, so a kill is data;
* the machine it was measured on: total RAM, cores, MySQL buffer pool size, and the code version (epoch). A profile from another machine shape or an older epoch is used at half trust and replaced by the first fresh samples.

The existing writers feed it (reuse): `memorygate.release` already knows the peak; `unitrates.note` and `stagetimes.record` already know rate and wall time; `units.horse_killed` already knows a memory kill. No new sampling is added. The Redis hashes stay as the fast cache the running code reads.

### Record first (Boss, 2026-11-12, 2:52 PM EST)

Step 1 of the build changes nothing about how jobs run. For every finished unit, of every command and queue, it stores one sample: the kind, the items the unit did (`unitbeat` already counts them), the unit's own worked seconds (without gate and lock waits, as `unitrates.note` already measures), the peak MB (`peak_resident_mb`), and from those **seconds per item** and **MB per item**. The profile keeps per kind both the per-unit figures and the per-item figures, because the two answer different questions: MB per unit sizes a claim for a unit of a known size, MB per item says how memory grows with the split, which is what a later split choice needs. MB per item is read as the slope between size buckets, not peak divided by items alone, since a unit also holds a fixed base (the worker with its tables is about 1.4 GB before the first item). So each sample also carries the unit's size bucket, and the profile fits `peak = base + slope x items` once a kind has samples from at least two buckets; until then only the per-unit figure is used. Jobs that are not units (API jobs, page jobs, the cut) store their peak and seconds with items left empty.

### What each consumer reads

1. **Gate reservation (claim).** `memorygate.claim_mb` becomes: the profile's 95th percentile times a margin that shrinks with more samples (1.25 at 3 samples down to 1.10 at 30), else the typed figure. The typed figure stops being a floor once the profile has enough samples, so an over-generous constant is corrected downward on this machine. A memory kill still wins: the retry at 1.5 times, and the kill is recorded against the bucket so the next run starts higher.
2. **Worker counts.** `workers.stage_worker_cap` is worked out from the claim the gate would now make, not from `STAGE_WORKER_MEMORY_MB`. The share of the host a stage may take stays a setting. `adaptive.py` keeps its role (react to waits within a run) but its starting limit for a run is the last run's settled limit for the same stage kind, kept in the profile, instead of the configured cap. Its edges are tried against recorded waits (OPS.151) rather than typed.
3. **Splits.** At plan time (the planner writes the cut into the run record, so it is fixed for the run), a stage asks the profile for a target unit size: the number of items that took about `target_seconds` at the recorded rate, and whose memory fits `claim x the worker count the stage may have`. The size is clamped between a floor and a ceiling kept next to the typed constant (the constant becomes the **default and the clamp**, so it is never removed). `coarsen` and `MAX_UNITS` still apply last.
4. **Times.** Unit expected time, the "not timed yet" estimate and the heartbeat judge's expected rate read the profile before `roughtime`. The 300 s ceiling stays a setting; the profile only picks the size that lands under it.
5. **Host edges** (reserve, load limit, wait limit): these stay typed. They are safety edges, not workload measurements. The only automatic part is already there (the adaptive lowering under memory pressure).

### Fallback when no history exists

In this order, per field: profile with enough samples (at least 3) on this machine shape, then profile from another shape at half trust, then the typed constant that is in the code today. A fresh install, a new job kind, a wiped control database or a new epoch therefore behaves exactly as the build does now, and the first run measures. `forget` clears one kind or all. A setting `queue.self_tuning` (default on once built, off restores today's behaviour exactly) and the existing `PLANETGEN_MEMORY_GATE=0` switch both turn it off. A reading that is missing or implausible (zero seconds, a peak below `MIN_CLAIM_MB`) is dropped, never stored.

### Safety rules

* **Output must not change.** `work.task_seed` keys a task, not a worker, so worker count and order are safe. A changed **cut** is only safe for a stage whose tasks are seeded per item, not per unit index. Every stage gets a flag, `cut_neutral`, set only after a test shows the stage's golden is byte-identical at two different unit sizes. A stage without the flag keeps its typed split and only tunes memory and workers. This is the one real risk of the scheme and the reason splits come last.
* **Frozen per run.** Splits and the claim a run started with are written into the run record at plan time (plan then execute, DOC.113), so a restarted unit uses the same cut and a run is reproducible from its record.
* **Slow to move.** A learned figure changes by at most 25 percent per run, and only after the unit finished `done`. A run that was stopped or failed adds nothing (same rule as `unitrates`).
* **Learning never blocks a job.** Every read and write is wrapped like `adaptive.adapt`: on any failure the typed figure stands and the job goes on.
* **Visible.** The Queue page and admin show, per kind: typed figure, learned figure, samples, last change. A change is logged with the old and new value and the run that caused it. Boss can pin or clear any value by hand.

## Order of work

1. `job_profiles` table (one row per sample, rolled up by a small job so it cannot grow without bound), migration, and the writers fed from the existing release, rate and stage-time paths, storing memory and time per unit and per item as under "Record first". No reader yet, so nothing changes in a run. Test: a small run fills it for every command it ran. This ships first and alone so the data is gathering while the readers are built.
2. Gate claim reads the profile with the fallback chain (memory only). Test: no history gives today's claims; history lowers a claim and a kill raises it.
3. Worker caps (`stage_worker_cap`) and the adaptive starting limit read the profile. Test: caps equal today's with no history.
4. Status and admin display, pin and clear (extends OPS.150).
5. Planner reads the profile for split targets, **one stage at a time**, each behind its `cut_neutral` golden test. Start with the stages that already scale by row count: Stage 17 uncharted index, block lights, link, settle, warm-up chunks.
6. Times: expected seconds and the heartbeat judge read the profile.
7. Retire the typed tables only where the profile has run on Boss's server for a full galaxy; the constants stay as the default and clamp (no shims, but no deletion of the fallback either).

## Open questions for Boss

1. Table in the control database (recommended) or Redis only? Redis only is simpler and loses everything on a flush.
2. Should a profile from a different machine shape count at half trust (recommended, so a rebuilt lane machine is not cold), or start from the typed figures?
3. Splits last and only for stages shown cut-neutral (recommended), or tune memory and workers only and leave splits typed?

## Items to file

Designations to be given by TODO. Suggested split: DOC (this note); OPS profile table and writers; OPS gate claim from profile; OPS worker caps and adaptive start from profile; OPS Queue page and admin display with pin and clear; PERF per-stage cut-neutral golden test; PERF split targets from profile (one item per stage kind, starting with Stage 17); OPS times from profile; OPS retire typed tables after a full-galaxy run. Related open items: OPS.149, OPS.150, OPS.151, OPS.154.
