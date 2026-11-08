# Every run as a graph of single small jobs, built from what Boss asks for

Written 2026-11-08, 3:35 PM EST by the job-queue design thread (DOC.113), for Boss. Design only: no code was changed. All times are EST. Evidence tags: [S] read in the repository on main (c406c9a), [R] read in an earlier report, [A] my arithmetic or judgement, [?] not verified, to be measured before a build. Item numbers (OPS.122 to OPS.127) are TODO's, filed in docs PR 5750.

## 1. Answer first

1. **You are right about how it works today.** `planetgen galaxy` is a fixed, hand-ordered list of 13 stages (`generation/stages.py` `galaxy_stages`), turned into steps one after the other (`commands.build_steps`). Each stage is cut into chunks of many sectors, and a stage must be wholly finished before the next is cut. The cut, the preparation, the closing work and the start of the next stage are not jobs at all: they run inside the last chunk's work process. That hidden work is where the three kills of 10-31 and 11-08 happened. [S]
2. **The fix is two phases for every task: plan, then execute (your 3:29 PM rule).** First a planner turns your request into an ordered **sequence of commands**, written down as data before anything runs. You say what you want (a galaxy of this shape, with these layers, or one scatter, or one sector); the planner works out which commands that takes, drops the ones whose result is already stored, and puts the rest in order, each with what it needs first. Then an executor runs that sequence and does nothing else: it never decides what to run, only when it may start (its predecessors are done, memory is free). The sequence is the contract: by default entry 5 starts after entry 4 is done. There is no fixed stage list in code. Letting independent entries overlap is a later, optional step (section 7, step 7), only where shown safe.
3. **Every entry in the sequence is a command that does one thing (your 3:32 PM rule).** A job, whether a Generate page job, an API job or a command line run, is nothing more than such a sequence. Today a command bundles several things behind hooks: the finish of the phenomena scatter draws the special phenomena, writes them, builds five indexes, stamps an epoch and records a seed in one go; the population pass scans, prunes, dates and founds in one go. Those become separate commands in the sequence, each with its own name, result and memory claim (section 3.0).
4. **Every job has a memory claim, and one broker admits jobs by the live free memory.** The claim starts from the measured figure for that kind of job and is corrected from what jobs of that kind actually used on this galaxy. A job that would not fit waits. A job that outgrows its claim is stopped by us, put back with a bigger claim, and run again, before the kernel has to kill MySQL. Small jobs make this cheap: the loss is a few seconds, not the 2 hours the last kill cost.
5. **A fair scheduler shares the machine across all kinds of work (your 3:44 PM EST rule).** Site-serving jobs have reserved places and always go first; the rest share the host by weight across job kinds (generation, warm-up, checks, API jobs); and waiting time raises priority, so nothing starves. The warm-up or any other command can never block the site, and every task that can run gets a turn. Small jobs are what make this work: a turn comes up every few seconds, not every hour (section 4.5).
6. **Jobs can spawn jobs.** An opening job reads the database and spawns the item jobs; a closing job is held until they are all done. Spawning is windowed (a rolling few thousand at a time), so a galaxy with tens of millions of sectors never puts tens of millions of jobs in Redis.
7. **It is built in small steps, each one a whole replacement with no old path left beside it, and each compared byte for byte against the run before it.** The first step (cut and close work as jobs of their own, with a memory claim) is the one that stops the kills, and it changes no output.
8. **Two decisions are yours** (section 8): (1) whether "single small command" means strictly one queue job per sector, or one command per sector that a worker may take several of in one go when each takes far less time than the cost of starting a job. I recommend the second; the first is a setting away. (2) Strict sequence or overlap; I recommend strict sequence for now.

## 2. How a galaxy is built today, and where the grouping is

All [S] unless marked.

**The stage list is fixed in code.** `galaxy_stages(args)` returns, in this order: the sector stages (start and neighborhood, or one sectors stage, picked by flags), `link`, `mass`, `luminosity`, the three phenomena stages (`phenomena`, `phenomena-index`, `nebula-draw`), `block-lights`, `backfill`, `population`, `territories`, `settle`, `run-blocks`. Order is the list order. A stage the request does not need is dropped or marked skipped by an `if` in code, and a run carried on after a stop works by tagging stages `finished_by` (`finish_run.py`).

**The only dependency data is four lines.** `commands.REQUIRES` names four edges (phenomena index after phenomena scatter, territories after population, link and settle after the three sector commands). Those edges are only used to refuse a wrong list (`validate_order`) and to write a "released" line. Nothing runs two stages at once; the doc says so (`redis-and-the-work-queue.md` 5.13: "making two stages run together is a later step").

**Grouping, layer by layer:**

| Where | What is grouped | Evidence |
|---|---|---|
| Stage barrier | Every unit of stage N must be done before any unit of N+1 exists. The unit that finishes the set queues the next stage. | `units._end_stage`, `wrapper.advance` |
| Chunks | A unit is many items. Sector stages: up to 4,000 chunks of consecutive sectors (`SECTOR_UNITS_MAX`). Link: 500 sectors a unit. Settle: 100 a unit. Backfill: at least 10,000 sky cells a unit, at most 2,000 units. Cut is fixed when the stage opens, by `static_chunks`, from the weights, never from memory or load. | `run_galaxy.expand_generate_sectors`, `units.static_chunks` |
| "Whole" commands | One job doing a large job: plan skeleton, finish phenomena, run population, run territories, update block totals, clean nebula field, reset. 12 commands are `whole`. | `commands.REGISTRY` |
| Hidden steps that are not jobs | `prepare`, the cluster-field build, `expand` (the cut), `finish`, `conclude`, `finish_stage`, `advance`. They run in the work process of the last unit of the stage before. | `wrapper._open`, `units._end_stage` |
| Followers | The process that started the run holds the display objects and some in-process stage work between stages. | `generation/followed.py` |
| Memory | A table of one figure per command kind (`workers.STAGE_WORKER_MEMORY_MB`, 700 / 2,100 / 3,050 MB) and a cap of half the host. Nothing reads the host while a run goes. The cut is outside the table. | `workers.py`, [R] `memory-kills-and-dynamic-workers-2026-11-08.md` |

**What that cost.** The 11-08 kills: 15 workers held 15,128 MB of 15,992 MB, the cut of the next stage ran in a work process for 1 hour 54 minutes with no log line, and the kernel killed it. [R] The job that was killed was named as the stage before, because the cut had no job of its own. [R] `perf252-opening-oom-2026-11-08.md`.

**What already works and is kept.** Every command is already a registered, pure, restartable command with JSON in and JSON out (`CommandSpec`). Every draw already comes from a seed of the thing it draws, not from the worker or the chunk (`work.task_seed`, GEN.39), so results do not depend on worker count. A unit that runs twice counts once. The run record in Redis is the only store of a run, and a second process can follow it. Stage statistics, ETA, Stop All, pause and the update guard are all built on it. This design reuses all of that; it does not start over.

**The earlier plan this continues.** Boss, 2026-10-17, 4:12 PM EDT (`redis-and-the-work-queue.md` 5.12): every part of generate is an individual command, a wrapper makes the list of commands and puts them into the queue individually, no command knows the order, "a chunk of a chunked stage is its own command". Built as PERF.106 to PERF.128 and PERF.172. PERF.247 (5.13) added "a job can require another job". PERF.254 (memory governor; the Dask assessment `perf254-dask-assessment-2026-11-08.md` is kept as the record) added the memory half. This design finishes the thought: the wrapper was meant to *select* commands from the request, but selection today is a fixed list with `if`s, and the "individual command" stopped at the chunk.

A note on a reference in the brief: the 11-06 report "Review phases and trains for efficiency" (`phase-and-train-efficiency-2026-11-06.md`) is about build phases and merge trains (the test gate, sharding it), not the job queue. Nothing in it changes this design.

## 3. The job model

### 3.0 The one-thing rule

Boss, 2026-11-08, 3:32 PM EST: all jobs are executed as a sequence of singular commands that each do one thing. This is the governing principle of 2026-10-17, 4:18 PM EDT ("each command does one and only one thing and does it exceptionally well") taken to its end.

**The rule.**

1. A command has one verb and one named effect: *clear the phenomena scatter*, *cut the mass scatter into shells*, *scatter one shell of massive stars*, *store the mass limit*, *build one index*. Its name says it, and its registry description needs no "and" or comma list.
2. A command has **no hooks**. `prepare`, `expand`, `finish` and `conclude` are removed from `CommandSpec`. Each of them is a command in the sequence, in its place, or it is not a command's business at all.
3. **Bookkeeping is the executor's, not a command's.** Timings, speed statistics, progress, the run's log lines and the stage row are written by the executor from what every job reports. `finish_scatter` and `finish_scatter_phenomena`, which only record speeds and metrics, disappear into the executor.
4. A **job** (Generate page job, API job, `planetgen galaxy`, `warm`, `check-db`) is a plan: a sequence of such commands. The single-verb command line runs are plans of one entry.
5. A guard test fails a registered command that has a hook, a second effect (a registry field `effects` names exactly one), or a description that lists several actions. It sits beside `test_one_wrapper_run`.

**Where the code bundles several things today** (read from the functions' own docstrings [S]; each split is a candidate for the sequence, with Boss's veto):

| Today | What it does in one go | Becomes, in sequence |
|---|---|---|
| `prepare_scatter_phenomena` | empties the earlier scatter, then takes off the indexes | `clear-phenomena` (exists as a command already), `drop-phenomenon-indexes` |
| `finish_phenomena_unit` (the `finish-phenomena` command) | draws the special phenomena (nucleus, hypervelocity stars), writes them, builds the table's indexes one at a time, stamps the hypervelocity stars' epoch, records the scatter's seed and mass cut | `draw-special-phenomena`, `write-special-phenomena`, `build-index` (one per index), `stamp-hypervelocity-epoch`, `record-phenomenon-scatter` |
| `conclude_scatter_phenomena` | counts the phenomena per class and stores them | `store-phenomenon-class-counts` |
| `prepare_scatter_pass`, `conclude_scatter_pass` | empties bright stars (or the luminosity pass); stores threshold, seed and mass limit; stores the generator epoch | `clear-bright-stars` (or `clear-luminosity-pass`), `store-scatter-limits`, `store-scatter-epoch` |
| `prepare_draw_nebula_field`, `conclude_draw_nebula_field` | checks the skeleton and empties the table; links the stored nebulae to their rows and makes the nebula block totals | `clear-nebula-field` (exists), `link-nebulae-to-field`, `make-nebula-block-totals` |
| `conclude_field` (`warm-cluster-field`) | makes the field from rows and parts, and publishes the file | `make-cluster-field`, `publish-cluster-field-file` |
| `_ensure_cluster_field`, called by the wrapper before a stage that reads the galaxy's shape | makes the file if absent | `make-cluster-field`, planned before the commands that need it (a declared fact: `cluster-field`) |
| `population_unit` (`run-population`) | scans the life worlds, drops species with no civilization, dates the civilizations, founds or dissolves polities | `scan-life-worlds`, `drop-civilization-less-species`, `date-civilizations`, `update-polities` |
| `conclude_rescore_habitability` | makes the blocks' habitable totals follow the sector counts | `make-habitable-block-totals` |
| `expand_*` (every cut) | reads the database and decides the items | `cut-<command>`: one thing, makes the item list |
| `generate_sectors_unit` | draws a sector and saves it in one transaction, and builds the log summary | stays **one command**, `generate-sector` (draw and save are one effect: making a sector exist; splitting them would move a whole sector through Redis). The summary lines are bookkeeping and move to the executor. Boss may overrule. |
| `refresh_block_lights_unit`, link, settle, one shell of a scatter, one backfill cell | one effect each | unchanged, but as items (section 3.5) |

The commands that have no hook today (checks, `reset-galaxy`, `plan-skeleton`, the warm families) are mostly already singular; `plan-skeleton` and `reset-galaxy` are one effect each. A few `check-*` and `warm-*` commands may also bundle (the audit in step 0 lists them).

**What a plan looks like then**, for `planetgen scatter --mass` (stored first, then run):

```
1  clear-massive-stars          (only if a re-draw was asked for)
2  cut-mass-scatter             makes: the list of shells
3  scatter-mass-shell  x N      one per shell, spawned in windows, run together
4  store-scatter-limits         needs: all of 3
5  store-scatter-epoch          needs: 4
6  make-block-lights  x blocks  needs: 5
```

Every line is one thing, can be run alone (`planetgen command run`), has its own claim and its own log line, and is named in a kill message if it is the one that dies.

### 3.1 A job

A job is one registered command applied to one item. It is a record in Redis:

| Field | Meaning |
|---|---|
| `id` | Hierarchical and stable: `<run>/<parent path>/<kind>:<item key>`, for example `r7/sectors/generate:3.12.401`. The item key is the thing worked on (an address, a shell, a sky cell), never a chunk number or a worker. |
| `kind`, `args` | A command in `REGISTRY` and its JSON arguments, as now. |
| `seed key` | The item key. The draw is `work.task_seed(run_seed, item key)`, as now. |
| `needs` | What must be done first (section 3.3). |
| `makes` | The facts this job's result stores (section 3.3). |
| `claim` | `mem_mb`, `est_seconds`, and named tokens (`mysql_write`, `star_tables`, `alone`) (section 5). |
| `parent`, `group` | Who spawned it, and the group (the old "stage") its bar and ETA belong to. |
| `state` | `planned`, `held`, `ready`, `admitted`, `running`, `done`, `failed`, `canceled`. |
| `result` | JSON, kept only as long as a dependant or the run's summary needs it. |

A command stays one thing and knows nothing about order. The five hooks of today (`prepare`, `expand`, `unit`, `finish`, `conclude`) are not three job kinds of one command, as an earlier draft of this note had it. They become **commands of their own in the sequence** (section 3.0): the cut is `cut-<command>`, a clearing is `clear-<thing>`, the unit is the item command, and a conclusion is one `store-...` or `make-...` command per effect. The cut therefore stops being something a worker does on the side: it is a command with a name, a claim, a log line and a place in the Queue page, and when it is killed the kill message names it.

### 3.2 Plan first: the planner writes the command sequence

**Plan, then execute, for every task** (Boss, 2026-11-08, 3:29 PM EST: "all tasks should be planned with commands in sequence and then should execute the command sequence"). This is the shape of the job spec already designed in `one-run-per-job.md` (PERF.203: a JSON list of `{command, argv}`, one run per Generate job), carried through to every kind of run.

* **The plan is data and has no side effects.** `plan(request, stored state)` returns the ordered list of command entries: `position`, `command`, `args`, `needs`, `makes`, `claim`, `group` (its label for bars). It is stored in the run record before the first job is queued, so it can be printed (`--plan-only`, and the Generate page shows it before you press start), estimated (the size and time estimate come from it, not from a separate path), and compared with what ran.
* **One planner for every start path.** The Generate page, the API's jobs, `planetgen galaxy`, `scatter`, `sectors`, `warm`, `check-db`, `sector`, `system`, and the jobs `update.sh` starts all build their request and call the same function. The page and the API cannot drift, which is what PERF.203 set out to guarantee.
* **The executor runs the sequence and nothing else.** It reads entry by entry, releases an entry when what it `needs` is done and memory is admitted (section 4), and records the result. A command still knows nothing about order. A running command's **open** job may expand into its item commands; those are appended to the plan record under their entry as they are written (in windows), so the record ends as the complete list of what ran, and the dry run shows them as "about N items".
* **A plan is re-plannable.** If a run stops, plan again from the stored facts: finished entries are not planned, the rest are. That replaces "finish a stopped run" as a special path.

Input to the planner: what you are doing, as a **goal**, not as a stage list.

* "Make a galaxy of this shape and these layers" (the Generate page, `planetgen galaxy`).
* "Scatter the bright stars" (`planetgen scatter --luminosity`).
* "Generate this block", "regenerate this sector", "re-score habitability", "warm the site", "check the database".

The planner (replacing `galaxy_stages`, `for_command`, `KIND_OF_STAGE`, `build_steps` and the `if`s that skip stages) does three things.

1. **Names the facts the goal needs.** Every command declares what it `needs` and what it `makes`, as named facts. My reading of the galaxy, from the code comments and the existing `REQUIRES`; each pair must be proven by the data-touch check in section 7 before any overlap is allowed [?]:

   | Command | Makes | Needs |
   |---|---|---|
   | plan-skeleton | `skeleton` | nothing |
   | generate-start-sector, generate-neighborhood, generate-sectors | `sector(address)` | `skeleton` |
   | link-neighbors | `links(address)` | `sector(address)` and `sector` of each neighbor |
   | settle-sector-paths | `path(address)` | `sector(address)` (the existing REQUIRES has all three sector commands; per sector is enough) |
   | scatter-mass | `stars:mass` and the mass limit | `skeleton` |
   | scatter-luminosity (and band) | `stars:bright` | `stars:mass` (it works at the stored mass limit) |
   | scatter-phenomena | `phenomena:shells` | `skeleton` |
   | finish-phenomena | `phenomena:index` | `phenomena:shells` (existing REQUIRES) |
   | draw-nebula-field | `nebula` | `skeleton`, `phenomena:index` [?] |
   | scatter-backfill | `stars:backfill(cell)` | `sector` of the cell's sectors, `stars:mass` |
   | refresh-block-lights, refresh-run-blocks | `blocks` | every fact that writes a block total (PERF.251 already wants it last) |
   | run-population | `population` | `sector` (all) |
   | run-territories | `territories` | `population` (existing REQUIRES) |

2. **Subtracts what already exists.** A fact already stored for this galaxy (same seed, same build) needs no job. This is what `skip`, `finished_by`, `_refuse_stored_scatter` and "run again after a stop" each do by hand today. One rule replaces them: *a job whose `makes` are all present is not planned.* A run carried on after a stop is just a new plan over the facts that were saved.

3. **Emits the sequence, lazily below the command.** The plan holds the commands. Each command's items (a range of sectors, a shell list, a block) are described compactly and written by the command's open job in windows (section 3.4), not stored up front.

The order of a run is therefore a result of the planner, not a list someone typed. Add a command with `needs`/`makes` and it fits; ask for less and fewer jobs exist. The Generate page, the API and the command line all send the same thing (a goal and its settings) to the planner. The "Stage 4 of 13" bars become **groups**: each command in the plan is a group, bars and ETAs are the sum over the group's jobs, so `stagetimes`, the ETA and the Queue page keep their meaning.

### 3.3 Dependencies

* **Counters, not RQ `depends_on` for wide fan-in.** RQ's `depends_on` is quadratic when one job waits on N: 2.3 s, 10.1 s, 45.4 s for 200, 500, 1,000 units. [R, measured] The linear way, 7 ms a unit, is an atomic counter in Redis. So each held job has a `waiting` count; each fact has a list of dependants; a finishing job runs one Lua script that decrements them and moves any that reach zero to `ready`. That is the done-set barrier of 5.1, made general.
* **Two shapes of edge.** *Narrow*: link of sector S waits for S and its (at most 8) neighbors; settle of S waits for S. *Wide*: backfill waits for the generated sectors of its cell, or population waits for all sectors; a wide edge is one group counter, and it is the old barrier, now only where it is real.
* **Failure.** A failed job fails its dependants with the line "required job X failed", as `requires.job_failed` does now, and the run ends with the job named. Retry once (as now); the second failure stops the run.
* **Cancel, pause, Stop All.** Cancel is by id prefix (`<run>/...`), which is why ids are hierarchical. Pause stops admission only. Stop All and the update guard read the same run record.

### 3.4 Spawning and the size of the queue

A job may spawn jobs: `spawn(children)`. The only ways a job can create work:

* **Cut commands** (`cut-<command>`, today's `expand`) read the database and spawn the item commands. A cut does that one thing and nothing else.
* **Windowed spawners.** The cut of a command with a huge item count does not write all items. It writes its cursor (`next = 0, total = N`) and the first window (for example 4 times the workers, never more than 2,000 queued). Each finishing item runs the same Lua script that releases dependants; if outstanding items fall under the low-water mark and the cursor is not at the end, it also asks the spawner to write the next window (a tiny continuation job, or the same script if the cut is a pure range). Redis then holds at most one window per command whatever the galaxy size. [A] Measured numbers: 1,831 bytes per queued job, 7,229 with the result kept [R]; with `result_ttl=0` as now and a window of 2,000 that is under 4 MB.
* **Joins.** The `store-...` command that follows needs the group counter to reach `total`; the counter counts items done, including ones not yet spawned, so the join is correct while the window rolls.
* **Depth is small and fixed by design:** run, command (cut, items, then the `store-` and `make-` commands that follow). A spawner may not spawn a spawner of its own kind without its claim being counted, so a bug cannot multiply work.

### 3.5 How small is "small"

The item is the smallest piece that has a meaning of its own, a seed of its own and a result of its own:

| Command | Item today (chunk) | Item in the new model |
|---|---|---|
| generate sectors | up to thousands of consecutive sectors | one sector |
| link | 500 sectors | one sector |
| settle | 100 sectors | one sector |
| scatter mass / bright / band / phenomena | a run of shells | one shell (or one shell slab if a shell is too small to be worth a job) [?] |
| backfill | at least 10,000 sky cells | one 50 pc sky cell |
| block lights | one of 90 units | one block |
| nebula field | a slab of cells | one cell slab |
| skeleton, reset, checks, finish phenomena | one job | unchanged: they are one thing |
| population, territories, run blocks | one job | split by sector later, once the sectors are items (section 7, step 8) |

Cost of a job: an RQ fork-per-job worker runs a no-op job in 14.2 ms (70 a second); the done-set design measured 7 ms a unit end to end. [R] That is nothing next to a sector (seconds; a generate unit is declared 20 s for a chunk). It is a lot next to a link or a settle of one sector [?, not measured]. That is why section 8 asks the question: the command is always one item; the only choice is whether a worker may take several tiny ones in one go. Because every item has its own seed and result, packing changes no output and can be chosen at run time from measured speed, never from the cut.

## 4. Memory and fair share: the broker and the scheduler

### 4.1 The claim

Each job carries `mem_mb`. It starts as the larger of

* the figure in `STAGE_WORKER_MEMORY_MB` for the command (700, 2,100 or 3,050 MB; cut and `store-`/`make-` commands get figures of their own, which is new), and
* the learned figure: the 90th percentile of `units.peak_resident_mb` (proportional size, shared pages once) for that command on this galaxy, kept in Redis and refreshed by every finished job. The stage report already records the biggest worker per stage [S]; this keeps the whole curve and ties it to galaxy size (sectors, shells), so a bigger galaxy raises the claim *before* a job runs.

A packed job (section 3.5) claims `base + per item x items` where per-item is learned.

### 4.2 Admission

The broker is one function that runs atomically over a Redis hash of live claims:

```
available = MemAvailable (and the cgroup limit, whichever is lower)
          - reserve                         (MySQL, Redis, web: the size they have now, not a guess)
          - sum over running jobs of max(claim, last observed size) - observed size     (headroom still to be used)
admit job if claim <= available, and its tokens are free
```

* Workers do not take whatever RQ gives them. A ready job goes to a lane by memory class (lean, tables, sector, alone); a worker of that lane asks the broker before it runs the job it was given, and if refused hands it back to the front of the ready list. [A] RQ cannot pop conditionally, so the honest options are this hand-back, or our own ready list in Redis with a `claim_next(free_mb)` script. I recommend the own ready list once jobs carry claims, because it also gives best-fit and age priority; RQ stays the worker runtime and the interactive and long lanes are untouched.
* **Token slots** for other shared resources: `mysql_write` (a cap on concurrent heavy writers, section 6), `star_tables` (the 1.4 GB tables are built once per worker), `alone` (nothing else runs; this is what the cut of the stage needs, and what PERF.226 imitates by shedding workers).
* **Aging.** A job that has waited long gets priority over smaller later ones so a big claim is never starved by small ones.

### 4.3 Adjusting by load

Every 2 to 5 seconds each running job reports its proportional size (psutil, or `/proc/<pid>/smaps_rollup`). The broker then:

1. keeps `max(claim, observed)` for the sum above, so a job that grew uses up the room it really takes;
2. when available falls under the reserve: stops admitting, asks the newest *idle* workers to leave (the `shed_bulk` call, now by load);
3. if it keeps falling: **stops the running job with the most growth that has the least done**, puts it back with `claim = observed x 1.3`, and logs "job X was put back: it held N MB against a claim of M MB". Items are small and restartable and a unit that runs twice counts once, so this is cheap, and it is a stop *we* choose rather than a kernel kill of whatever it picked (MySQL, last time);
4. a job whose claim alone exceeds what the host can give runs `alone`; if even that is not enough it fails with a message naming the job and the MB, because no governor shrinks a pass that wants 57 GB. [R] The remedy then is the size of the job, which is the argument for items.

Everything here is code in the program. Nothing asks you to change a systemd drop-in, Apache limit or MySQL setting.

### 4.4 Relation to PERF.254, and Dask is out

* **PERF.254 (rewritten by TODO, 2026-11-08, after Boss's decision)** is the first slice: RQ stays the job runner, and our own psutil reader plus a Redis memory gate admits and sheds workers from the live free memory. Boss's rule is that **all jobs of every kind** run under that gate (generation units, warm-up, checks, API jobs), not only the stages of a galaxy run. That is this broker: the claim, the live reading and the admission rule are the same code, and this note's per-job claims are what the gate takes in place of a per-stage figure once OPS.122 lands.
* **Dask is dropped** (OPS.125 withdrawn). The graph, the planner, the claims and the broker live in our Redis, because the run record is what the Queue page, the update guard, Stop All, pause and finishing a stopped run all read. Lane 22's assessment (`perf254-dask-assessment-2026-11-08.md`) stays as the record of why it was considered. Nothing in this design needs it, and no step below depends on it.

### 4.5 Fair share: one scheduler for all job kinds (OPS.130)

Boss, 2026-11-08, 3:44 PM EST: the job system uses a task scheduler that shares system time fairly, so the warm-up or any other command never blocks the site, and every task that can run gets a turn (fair share across job kinds, priority and reserved places for site-serving work, aging so nothing starves).

**What exists today** [S]: a fixed slicing of the host's worker places. Three lanes (`bulk`, `interactive`, `long`) and a maintenance lane of three queues taken in turn by one worker (`workers.lane_of`, PERF.245); `reserved_places` keeps `queue.interactive_workers + queue.long_workers` places out of the bulk and maintenance lanes' reach; `maintenance_cap` gives upkeep `queue.maintenance_workers` while a generation run is live; a pause flag holds the bulk and maintenance lanes only; the site's workers run at a lower nice value. It is a good first cut (PERF.150 to PERF.153, PERF.245) but it is a **count of places per lane, set by hand**, it knows nothing about memory, and a single unit of a lane can hold its place for an hour.

**The scheduler.** The same single admission point as the broker (section 4.2), which now answers two questions for every free place: *which class of work goes next*, then *which job of that class fits*.

1. **Classes.** Every job kind belongs to a class: `site` (jobs a visitor or an admin is waiting on: the interactive and long API jobs, page builds), `run` (generation), `upkeep` (warm-up, checks, database repair). A class has a **weight**, a **reserved number of places** and a **priority**. They are settings (`queue.share.*`, `queue.reserved.*`), not code constants and not operating-system settings.
2. **Reserved places, first.** `site` jobs have places no other class can take, so a warm-up that has used every other place still leaves a visitor's search a worker at once. A reserved place with no `site` work in it is lent to the other classes and taken back as soon as a `site` job is ready; because lent jobs are small, taking it back costs seconds, not the end of a unit of an hour (today's lend rule, PERF.150, with the same idea).
3. **Fair share, for the rest.** Between `run` and `upkeep`, and between job kinds inside a class (a run of units and another run, the warm-up families), places go by **weighted fair queuing**: each kind keeps a running total of what it has used (claimed MB times seconds, and places times seconds), and the next free place goes to the ready kind that is furthest behind its weight. One run with ten million ready items does not shut out a warm-up with two jobs, and the warm-up cannot hold the run back past its weight.
4. **Aging.** The longer a ready job has waited, the higher its effective priority, up to a cap; a job with a large claim that has waited past a limit gets the next room that opens (the other classes are not admitted until it fits), so a big claim is never starved by small ones (also 4.2).
5. **Two resources, one rule.** Fairness is taken over places (cores) and over memory, because either can be the one that is short. A `mysql_write` token (section 6) is shared the same way.
6. **No preemption of running work.** A job is never cut short for fairness. What makes the share fair is that jobs are short (a few seconds to a minute at most); the longest turn anyone waits is the longest job in the class ahead. Jobs longer than the unit-length ceiling (`queue.unit_ceiling_seconds`) are already refused.
7. **Pause and Stop All** stay what they are, a hold on admission by class. The admin's worker count (`wanted`, "Ludicrous Speed") becomes the total number of places the scheduler shares.

**Replaces:** `reserved_places`, `maintenance_cap`, `maintenance_share`, the lend rule and the per-lane worker counts, in one change, with no old path beside them. **Related:** OPS.121, the warm-up priority work of PR 5582 (the rule that two warm-ups never run at once stays), PERF.150 to PERF.153, PERF.245. **Tests:** a fake clock and a fake host; a run, a warm-up and a check all ready at once must each get their weight within a bounded number of admissions; a `site` job waits no longer than one job's time when its place is reserved; a large claim is admitted within its aging limit while small jobs keep arriving.

## 5. How a galaxy run looks

Boss asks for a random-start galaxy of a given shape with all layers. This section shows the full design once overlap (steps 6 and 7) is allowed; with strict sequence (Decision 2, option 1) each command below simply starts after the one before it is done, and the items inside a command still run together.

1. The planner reads the request and the stored state. Nothing is stored, so it plans the whole goal.
2. The plan is stored as the sequence of commands, with root run `r7`. Only `plan-skeleton` is ready.
3. `plan-skeleton` finishes, its fact `skeleton` appears, and the cuts of generate-start-sector, scatter-mass and scatter-phenomena (all need only the skeleton) become ready. They run beside each other if the broker admits them. [Today they run one after another; running them together is allowed only after the check in section 7, step 7.]
4. `cut-start-sector` makes one sector item and `generate-sector` runs it; then `cut-neighborhood` makes the neighborhood's items in a window.
5. As each sector finishes, the link job for it waits only for its neighbors; settle of it waits only for itself. No barrier at the end of the sectors.
6. Backfill waits per sky cell for the sectors of that cell and for the mass scatter.
7. Block lights and the final block totals wait on every fact that writes a block total. Population waits for all sectors; territories for population.
8. The run is done when every command of the plan is done and the executor has written the run summary. Stop at any time and plan again: the stored facts are subtracted and the rest is planned.

## 6. Risks and how each is held

| Risk | Why it matters | How it is held |
|---|---|---|
| **Byte-identical output** | You need the same galaxy for the same seed, at any worker count and any packing. | Seeds already come from the item (GEN.39, `task_seed`). **Unaudited:** (a) any unit whose draws depend on chunk edges (a run of shells that shares a stream), (b) any row id taken in completion order (`id_blocks`, auto-increment), (c) any unit that reads rows another unit of the same stage writes. The audit is step 0 and lists them. Every migration step is compared with the run before it through the golden and the census and row hashes keyed by content (address, shell), with workers 1, 4 and 15, with packing off and on, and with a forced job kill and put back in the middle. A step that fails the comparison does not merge. |
| **MySQL write contention** | One sector is one transaction today. One-sector jobs raise the commit rate, link takes a neighbor lock in turn, and each job opens a connection. | The `mysql_write` token caps heavy writers; a worker keeps its pool between jobs; the lock for linking becomes per-pair, which is what narrow edges give (link of S only waits for its own neighbors). Measure commits a second and lock waits before and after, on the host's 16 GB, in step 0. The binlog and disk growth with many small commits is a known trap on lane machines. |
| **Queue size** | Tens of millions of sectors would be tens of millions of jobs. | Windowed spawning (3.4), compact group descriptions in the plan, `result_ttl=0`, and counters in place of rows for the Queue page. Redis never holds more than a window per command, about 4 MB at 2,000. |
| **Scheduler complexity and bugs** | A broker that is wrong stalls a run or kills it slowly. | Fall back to the table figure and the old cap when the broker cannot read the host; bound every wait; log every decision (admit, refuse, shed, put back) in the run's log; a fake host with memory curves in the tests; a stall detector that fails the run with the waiting jobs named after a limit. |
| **Starvation and deadlock** | Big claims behind small ones; circular needs; one kind of job taking every place. | Fair share and aging (4.5); `alone` jobs drain others; the planner rejects a cycle at plan time. |
| **Fairness against throughput** | Reserved places sit idle when no site work is waiting, and sharing slows a generation run while a warm-up runs. | Idle reserved places are lent and taken back in seconds; the weights are settings and the log shows each class's share, so a slower run is a number to tune, not a surprise. |
| **Retry loops** | A job that always outgrows its claim could be put back forever. | Two put-backs, then it runs alone, then it fails with its size in the message. |
| **Observability** | A run of millions of jobs cannot be listed row by row. | The Queue page shows groups with counts and the slowest and the held; a job page by id; the run log keeps the same lines. |
| **Progress, ETA, stage statistics** | They are keyed by stage. | A group is the old stage: same keys, same `stagetimes`, bars from group counters. |
| **Redis becomes the graph** | The run record already lives there. | `appendonly yes` and the ACL set up under SEC.39 and the install notes; the run stays safe to rebuild: the plan is a pure function of the goal and the stored facts, so a lost graph is re-planned, not lost. |
| **Replace, no shim** | Your standing rule. | Each step deletes the thing it replaces in the same change (the stage list, `REQUIRES`, `finished_by`, the hooks inside the horse). No switch stays beside it. |

## 7. Migration, one step at a time

Every step: one change, old way deleted, compared byte for byte with the step before, tested on the toy galaxy and then on the 16 GB host. "For Boss" is what you would see.

| Step | What changes | For Boss | Output |
|---|---|---|---|
| 0. Measure and audit | (a) per-job overhead on the real worker with a real one-sector job, (b) the cut of the backfill and of the other stages on a full galaxy (Lane 13's measurement, PERF.252), (c) commits a second and lock waits, (d) the three byte-identity audits of section 6. | A report. | none |
| 1. Hooks become commands, one thing each | `prepare`, `expand`, `finish` and `conclude` are removed from `CommandSpec`; each becomes a command in the sequence (`clear-*`, `cut-*`, `store-*`, `make-*`), with its own claim, log lines and `alone` where needed; the speed statistics move to the executor. Bundles are split as in section 3.0, one command per PR. Stages still run in the same order. `_end_stage` and the work done inside the last unit's horse are deleted. The guard test of rule 5 goes in at the end. | Kill messages name the real job; the cut shows in the Queue page and the log with its memory. This is the change that removes the hiding place of the 11-08 kills. | unchanged |
| 2. The broker | Claims, live reading of memory and per-job sizes, admission, idle shedding, put-back with a bigger claim, `alone`. `STAGE_WORKER_MEMORY_MB` becomes the starting figure. (This is PERF.254 step 1.) | Memory is read while the run goes; the log says why a job waited or was put back. | unchanged |
| 2b. Fair share | Classes, reserved places, weighted fair queuing and aging in the one admission point (OPS.130). Replaces `reserved_places`, `maintenance_cap`, `maintenance_share` and the lend rule. | A warm-up never starves the site, and a run, a warm-up and a check all make progress together; the log says which class got each place. | unchanged |
| 3. Job records and counters replace the barrier | Run record gets job ids, `needs`/`makes`, group counters and the Lua release. Stages become groups. `REQUIRES` is replaced by declared facts. Still one group at a time. | Same bars; the Queue page lists held jobs by group. | identical by construction |
| 4. The planner, plan then execute | The request becomes a goal; the planner writes the command sequence (stored, printable, estimated) and subtracts stored facts; the executor only runs it. This takes over PERF.203's job spec for every start path. `galaxy_stages`, `for_command`, `KIND_OF_STAGE`, `build_steps`, `finished_by` and the stage `skip`s are deleted. The Generate page and the API send goals. | The run is built from what you ask for. "Finish a stopped run" is just planning again. | identical |
| 5. Items as jobs, windowed | Sector generation first (the biggest), then scatter shells, backfill cells, link, settle, block lights, nebula slabs. Windowed spawners; optional packing (section 8). One command per change, each compared with its chunked self. | Progress per sector; a kill loses one sector, not a chunk. | identical |
| 6. Narrow edges | Link and settle per sector, backfill per cell, so they start while sectors are still being generated. | Shorter runs (no idle workers at barriers); the link lock is per pair. | identical rows, a different order of writing: proven by content hashes |
| 7. Overlap of independent commands | Allowed only for a pair of commands shown not to touch the other's data (the `makes`/`needs` table plus a test that fails a command that reads a fact it did not declare). First candidates: the skeleton-only group (sectors, mass scatter, phenomena). | Fewer idle minutes; more memory at once, which the broker meters. | identical |
| 8. Split the whole commands | Population into its four effects first (step 1), then by sector; territories, run blocks by block, once the sectors are items. | Population gets bars and a memory claim per piece. | identical |

Rough effort [A, not measured]: steps 1 and 2 are days (queue code, no output change; step 2 is PERF.254 step 1); step 3 is about one to two weeks; step 4 about a week; step 5 two to three weeks across the commands (one lane, sector generation first); steps 6 to 8 depend on what step 0 finds. Lane 22 estimated 2,500 to 3,500 of the 7,273 queue lines are in the path of a change to how jobs are submitted [R]; this design touches about the same lines (`wrapper`, `units`, `workers`, `runrecord`, `commands`, `stages`, `finish_run`, `followed`), so the cost is of that order. The test suite has 36 files that touch the queue.

Done when: a full galaxy on a 16 GB host runs past stage 13 with no kernel kill, with groups overlapping where proven safe, output equal to the chunked run, and the cut of every command seen in the log as a job with its memory.

## 8. Decisions

**Decision 2 (yours, small): strict sequence or overlap.** I read your 3:29 PM rule as strict sequence: each command starts after the one before it is done, as today, with the plan written first. Option 1 (recommended): strict sequence now; allow overlap of independent commands only later, one proven pair at a time, if you ask for it. Option 2: allow overlap from the start wherever the declared `needs` permit. Option 1 changes least and keeps every step byte-identical by construction.

**Decision 1 (yours): how small is one job.**

1. **One command per item, and a worker may take several tiny items in one go** (recommended). Sectors, shells, sky cells and blocks are always their own items with their own seed, result and claim. A worker takes more than one only when an item costs far less than starting a job (links and settles of one sector, probably). The number is chosen at run time from measured speed, never from the cut, and `queue.min_job_seconds = 0` makes it strictly one job per item.
2. **Strictly one queue job per item, always.** The simplest to explain. Costs about 14 ms of job start per item [R], which is nothing for a sector and a lot for a link; the queue would carry more jobs, still bounded by the window.

My recommendation is 1 for the reason in 3.5; nothing in the design prevents 2.

Smaller choices I took, which you can overrule: the graph lives in our Redis and Dask is out; PERF.254 step 1 is the first slice, not a separate project; the planner takes a goal and subtracts stored facts, so there is no "skip" flag anywhere.

## 9. Items (numbers from TODO, docs PRs 5750 and after; all proposed, build after Boss approves this note)

| Item | What | Owner proposed | Steps of section 7 |
|---|---|---|---|
| DOC.113 | This design note | this thread | |
| PERF.254 | Ships first. Rewritten by TODO: RQ stays; our own psutil reader plus a Redis memory gate admits and sheds workers; Boss's rule is that all jobs of every kind run under it | Lane 22 | 2 (first form) |
| OPS.122 | Job model: job record, declared dependencies, subjob spawning, per-job memory claim, replacing the stage barrier | Lane 22 | 3 |
| OPS.123 | Hooks (`prepare`, `expand`, `finish`, `conclude`) become single-purpose commands in the sequence | Lane 13 | 1 |
| OPS.124 | Memory broker: per-job claim, measured and adjusted by load; the gate of PERF.254 grows into it | Lane 22, after PERF.254 | 2 |
| OPS.126 | Sector generation as single-sector jobs fed by a spawner | Lane 13 | 5 (sectors first) |
| OPS.127 | Migration stage by stage, each compared byte for byte | Lanes 22 and 13 | 0, 5 to 8 |
| OPS.128 | Run planner: every task planned first as a stored command sequence, for every start path; takes over PERF.203's job spec | to be set by TODO | 4 |
| OPS.129 | Split the bundled commands (the phenomena finish into about five, run-population into four, and the rest of section 3.0) | to be set by TODO | 1 and 8 |
| OPS.130 | Fair task scheduler: classes with weights and reserved places for site-serving work, weighted fair sharing across job kinds, aging; replaces the hand-set lane counts (Boss, 3:44 PM EST) | Lane 22, after PERF.254 | 2b |
| OPS.125 | Optional Dask engine: **dropped** | | |

Not covered by an item yet: the narrow-edge and overlap steps (6 and 7). They would be filed after step 5 shows what the measurements say, and only if Boss asks for overlap (Decision 2).

## 10. What I did not verify

* The default galaxy's sector count and the time of one sector job; the commit rate at one sector a job; per-job cost for a link or a settle.
* The `makes`/`needs` table beyond the four edges the code already declares and the comments I read.
* Whether any current unit draws from a stream shared across its chunk, or takes row ids in completion order.
* The size of the cut of each stage on a full galaxy (Lane 13's open measurement).

## Sources (all in the repository at c406c9a unless noted)

`src/planetgen/generation/stages.py`, `generation/followed.py`, `generation/run_galaxy.py` (`expand_generate_sectors`, `generate_sectors_unit`), `queue/commands.py`, `queue/wrapper.py`, `queue/units.py`, `queue/workers.py`, `queue/requires.py`, `queue/runrecord.py`, `queue/pause.py`, `cli/worker.py`; `docs/design/operations/redis-and-the-work-queue.md` (sections 3, 5.1, 5.12, 5.13); `docs/TODO.md` (PERF.254, PERF.252, PERF.226); `/mnt/project-files/reports/perf254-dask-assessment-2026-11-08.md`, `memory-kills-and-dynamic-workers-2026-11-08.md`, `perf252-opening-oom-2026-11-08.md`, `phase-and-train-efficiency-2026-11-06.md`.
