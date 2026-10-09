# Generation performance study: where the time goes in the scatter and the sector fill

Profiles of the three generation phases on a scaled-down galaxy: the bright-star
scatter, the phenomenon scatter (GEN.100) and the sector fill, with the database path
split out, the word-salad share measured, ranked bottlenecks, and what each fix is
expected to save. It answers PERF.31 for the phases it covers and corrects one claim in
PERF.18. It extends [performance-eta-queue-and-caching.md](performance-eta-queue-and-caching.md)
and the bright-star report (`/mnt/project-files/bright-star-timing/report.md`, PR #657) and
does not repeat them.

Informs: PERF.31, PERF.18, PERF.39, PERF.34, GEN.100, GEN.72, GEN.64, GEN.69, GEN.126

Status: research, 2026-10-09; no generator code was changed. Everything marked "recommend" is a proposal for Boss.

Evidence tags: [C] measured or computed in this research (scripts under
`/mnt/project-files/research/scripts/genperf/`, outputs under
`/mnt/project-files/research/genperf-results/`); [S] read in the repository; [R] recalled
and unconfirmed; [B] Boss's words.

## Method and limits

- One 4-core Linux container, MariaDB 10.11 and Redis on the same box, Python 3.13,
  planetGen 8.0.711 (the project's minimum is Python 3.9, which is slower) [C].
- A quarter-scale galaxy (`--disk-scale-length-pc 650 --disk-scale-height-pc 75
  --bulge-scale-radius-pc 395`): 509 layers, outer ring 940, 4.6e9 expected systems,
  420,840 bright stars, 18,241,331 phenomenon rows [C]. The default galaxy is about 64
  times the volume [C].
- Wall times come from timers and from the sampling profiler (pyinstrument, 1 ms), not from
  cProfile, whose per-call cost inflates code made of many small calls (the word-salad
  share reads 7% under cProfile and 4% in real time).
- `--debug` logging was measured on the scatter and costs nothing there (36 s against 39 s)
  [C].
- What a quarter-scale run cannot tell is in "What this does not predict".

## The answer in short

1. **The bright-star scatter is dominated by per-job start-up, not by drawing stars.** A forked RQ
   work horse re-imports `planetgen.generation.run_plan` (2.4 s wall, 2.9 s CPU, of which
   `nltk` 0.87 s) and rebuilds the band tables (`star_population._bright_table`, 0.19 s per
   population) for every layer job. With the modules imported and the tables warmed before
   the fork, the same 420,840 stars take **39.5 s instead of 322 s** (8.2 times), identical rows [C].
   At this scale 95% of the CPU was start-up.
2. **The phenomenon scatter writes 43 times as many rows as the bright-star scatter** (18.2M
   against 0.42M here; about 1.17e9 against 26.9M at the default scale [C]) and each row is
   138 bytes, so about 161 GB at the default scale. Its per-row cost is small (32 to 37 CPU-us);
   its volume is the problem.
3. **The sector fill costs 20 to 27 ms per system** (3.4 to 7 s per sector) in one process.
   Python generation is about half of that, the database save about 41%, and SQL
   execution (with pymysql's client-side escaping) 31% of the wall [C]. The statement count is
   already low (1.3 to 2.2 per system, PERF.13 batching); the cost is the volume of rows
   (about 19 non-system rows per system) and the second pass that updates every one.
4. **Word salad is 4% of the sector fill and about four fifths of that is wasted.** A call
   costs 93 to 96 us warm and about 10 happen per system, 77% of them in
   `RoguePlanet.__init__`; but a placed phenomenon is named by its object ID (GEN.64) and the
   generated name is thrown away [C][S]. Stars and systems still need theirs. The scatters
   do not call it.
5. **I cannot map "10 hours for the bright stars" to one phase without his log.** The numbers
   that could produce it are below ("Where ten hours could come from").

**Correction (same day).** The default-scale phenomenon count above (1.06e9 rows, 146 GB) came from
`phenomenon_scatter.layer_expected`, which runs 10% low; 64 times the 18.24M rows actually drawn is
1.17e9 rows and 161 GB. [phenomenon-scatter-mass-cut.md](phenomenon-scatter-mass-cut.md) has the corrected
counts, traces the notes' 1.6e8 and recommends a mass cut that takes the table to about 2.7e5 rows.

**A figure not to rely on.** Bugfixes lane 1's "30 s of 84 s in `reserve_system_names`" was a contaminated
measurement: it profiled while another run was saving into the same database, so it measured lock waiting.
Alone, name reservation was 1.7 s of a 17.8 s dense core sector (about 870 names), about 10% (reported by the
TODO thread; PERF.49 starts by re-measuring). The profiles in this note were taken with nothing else writing
to the database.

## Findings by phase

### Bright-star scatter (`scatter_bright_stars`, one RQ job per layer)

| Run (420,840 stars, 509 layers, 4 workers) | Wall | CPU | Note |
|---|---|---|---|
| Default `planetgen worker` | 322 s | 1123 s user + 106 s sys | 2.2 CPU-s per layer job |
| Modules imported and tables warmed before the fork | 39.5 s | 110 s | same stars |

Per star, warm: draw 40 to 49 us, insert 38 us (executemany, batches of 10,000, commit per
batch). In `scatter_layer` the largest cost is `population_densities` (941 rings times 32
angle bins, about 30% of a layer), then `_sample_one_bright` and the cell check
`_point_in_cell` with `sector_address_at` [C].

`cli/worker.py` imports only `planetgen.queue.redisqueue` and `work`; RQ forks a horse per
job, and the horse imports the task module and fills the `lru_cache` from scratch every
time [S]. The same fixed tax applies to every queue job (fill tasks, API jobs; PERF.39 measured
2.4 s and 195 MB per API job). **PERF.18's remark that the CPU gain is small** holds for the
backfill, whose per-cell draw is 0.06 ms, but not for the scatter, where start-up is the cost.

### Phenomenon scatter (`scatter_phenomena`, GEN.100)

18,241,331 rows in 169.9 s with the preload (3.43M black holes, 14.8M neutron stars, 1,765
hypervelocity stars, 445 planetary nebulae, 306 supernova remnants) [C]. Per row, warm: draw
13.7 us, insert 18.6 us, total 32.3 us in layer 0 (95,258 rows). Expected rows per
system 3.6e-3 against 9e-5 bright stars, so 40 times the rows [C].

Insert split on 60,000 rows [C]: pymysql `executemany` 13.8 us per row; one hand-formatted
multi-row INSERT 10.1 us per row (client formatting 3.4 against 1.9 us); the same
`executemany` into a table without the address index 13.2 us per row. So the client's
escaping is about a quarter of the insert and the index does not matter at this size.

Row layout: `kind` VARCHAR(24), three BIGINT positions, three nullable DOUBLE velocities,
an 8-byte seed per row, ring/layer/slot columns; 1559 MB data plus 807 MB index for 17.9M
rows (138 bytes a row) [C]. Most rows are neutron stars and black holes whose velocity columns
are NULL or derivable.

### Sector fill (`generate_and_save_sector_at`, one process, normal density)

Five runs of 5 to 6 sectors picked as `fill_bench.py` picks them (weighted toward sectors
that hold bright stars, so denser than average): **20.3, 21.0, 21.8, 22.9, 26.6 ms per
system** [C]. The bright-star report measured 16.6 ms per system on an earlier version.

The largest run (1,569 systems, 35.9 s) [C]:

| Part | Share of wall |
|---|---|
| SQL execution, including pymysql's client-side escaping | 31% (11.1 s) |
| Word salad | 4.2% (15,593 calls, 1.50 s, 96 us each) |

The sampled profile of another run (391 systems, 16.4 s) [C]:

| Part | Seconds | Share |
|---|---|---|
| Generate the sector's systems | 8.2 | 50% |
| of which planets and moons (`Planet.__init__`, moons 2.5 s, orbital motion 1.5 s, spin axis 0.5 s) | 3.7 | 22% |
| of which rogue planets (salad name 0.94 s under the sampler) | 1.4 | 8% |
| of which placing bodies in the galaxy (`place_in_galaxy`, `add_system`) | 2.3 | 14% |
| of which `FillContext.bright_share` building the band tables (once per process) | 0.7 | 4% |
| Save the sector (`insert_sector`) | 6.8 | 41% |
| of which writing the held INSERTs (flush, triggered by the registry upsert) | 2.7 | 16% |
| of which unique IDs (`assign_uids`: SELECT back, UPDATE every row) | 1.05 | 6% |
| of which nearest-system links (`_add_sector_to_nearest`) and neighbour search for the location | 1.6 | 10% |
| of which building the rows in Python (`insert_star_system`) | 0.7 | 4% |
| Containment / sector stats / nebulae | rest | |

Statements per system (server counters): 1.3 to 4.4 in total across the runs (the
smallest sectors cost the most per system), 0.2 to 0.5 inserts, 0.5 to 1.5 selects, 0.5 to
2.2 updates [C]. Rows written per system in one 449-system run: 25 moons,
9 planets, 8 rogue planets, 1.3 stars, 0.7 comets; the uid pass then updated 11,933 moon rows,
4,353 planet rows and 6,958 rogue-planet rows (594 moon statements alone) [C].
Per row, the INSERT costs 70 to 90 us (moons 87, planets 85, rogue planets 70) and the
later UPDATE about 38 us.

## Ranked bottlenecks

Ranked by the share of the phase they take, at the scale measured.

| # | Where | Share | Fix | Saving |
|---|---|---|---|---|
| 1 | Scatter jobs: import and table build per forked job | 88% of scatter wall | Pre-import and warm in the worker | up to 8 times faster [C] |
| 2 | Phenomenon rows: 1.17e9 at the default scale, 161 GB | whole of GEN.100 | Compact the rows or derive them on demand | 3 times smaller; or no rows [C][R] |
| 3 | Sector save: row INSERTs, 70 to 90 us each, about 19 per system | 16 to 25% of fill | Fewer columns written per row; the escape cost; see below | 5 to 10% |
| 4 | Sector save: unique-ID pass updates every row just written | 6 to 9% of fill | Compute the uid in Python while the row is built | 6 to 9% |
| 5 | Sector generation: planets and moons, orbital position updates | 22% of fill | Fewer coordinate resyncs; make the `finite_domain` checks optional in bulk fills | 5 to 8% |
| 6 | Neighbour search and nearest links per sector | 10% of fill | Move to a batch pass; cache coordinates | 4 to 10% |
| 7 | Word salad for names that are discarded | 4% of fill | Name lazily | 4% |
| 8 | Band tables built in a cold process | 4% of a cold fill | Preload with #1 | 4% (per process) |
| 9 | pymysql client escaping of values | 3 to 4% of fill, 25% of a scatter insert | C driver or a formatter for numeric-only tables | 3 to 4% |

Savings are shares of the phase, not additive with the dependent rows; the sampler and the
timer runs differ by a few points.

## Recommendations

### Cheap wins (small code, no change in what is generated)

- **A. Pre-import and warm the worker before it forks.** In `cli/worker.py` import the
  generation modules and fill `_bright_table` for the four populations at the default
  threshold. Expected: scatter 322 s to 40 s here; at the scale of the bright-star report
  (74.6 CPU-min) about 20 to 38 CPU-min saved. The same fixed 2 to 3 s goes from every fill
  and API job (PERF.39). Risk: after `update.sh` a long-lived worker would run old code, because
  today the horse re-imports on each job; the update script has to restart the workers, or
  the horse must compare `planetgen.__version__` with the one the parent loaded and exit when
  they differ. Tables for a non-default luminosity threshold stay cold (cache them on disk
  keyed by version and threshold if that matters).
  **Built (PERF.42).** `cli/worker.py` `warm()` imports the task modules
  (`WARM_MODULES`) and, when the worker serves the generation queue, builds
  `_bright_table` for the four populations at every floor a scatter or
  backfill uses (the galaxy-wide floor and the four backfill tier floors).
  The per-job API workers (one queue each) only import. A worker that has
  pre-imported would keep serving an old release after an update, so before
  it waits for the next job it compares the release in `_version.py` on disk
  with the one it loaded and, when they differ, deregisters and re-executes
  itself (POSIX; a spawned horse on Windows is a fresh interpreter anyway).
  Measured on this build box, a forked horse's start-up (importing
  `run_plan` and building those tables) went from 2.6 s to under 1 ms; the
  one-off warm-up costs 2.7 s per worker. The tables are the same as a cold
  build's (a test compares them), so a seeded run gives the same rows. A
  floor outside that list stays cold.
- **B. Do not generate names that will be discarded.** `RoguePlanet`, `Comet`, `Nebula`,
  `AsteroidField`, `SupernovaRemnant`, `Quasar` and the compact remnants call
  `generate_phoneme_salad_name` in their constructors, but a phenomenon placed in a
  galaxy sector gets its object ID as name unless `name_given` (`store._reserve_phenomenon_name`,
  `_sector_object_ids`) [S]. Make the salad name lazy: draw it on first read when the object
  is not placed. Saves about 3.5% of a fill. Risk: anything that reads `.name` before saving must still
  get a name; keep `name_given` as is. I could not reproduce the hang Boss saw on rogue
  planets: warm, a name is 96 us; the likely causes are the cold NLTK load (0.87 s) in every
  forked job and the retry loop, both before the current checks.
  Built in PERF.43: `names.wordsalad.LazySaladName` keeps one draw at construction
  (the name's own seed) and draws the salad on first read, so an object's other draws
  don't depend on whether its name is read. 3,000 rogue planets build in 0.15 s
  instead of 0.60 s. The one draw in place of many changes a seed's galaxy once.
- **C. Write each row's unique ID with the row.** `assign_uids` selects every new row back
  and updates it. The uid is `derived_uid(seed, kind, parent uid, rank)` where rank is the
  row's order among its siblings by id, and ids come from `id_blocks` before the INSERT
  [S]. The builder knows all of them in memory. Removes the SELECTs and every UPDATE (about
  1.8 s of 11 s SQL in the 1,045-system run). Risk: the id order and the generation order must
  match, which `insert_sector` already relies on; GEN.69's idempotence ("a row that has one keeps
  it") must still hold for a re-save.
  **Built (PERF.44).** `Connection.execute` asks a `_UidIssuer` (set by `insert_sector`) for
  the `uid` of each INSERT into a table that has one, and adds it to the statement; the issuer
  counts ranks per parent as rows go in, which is their id order. A bright-sweep system keeps its
  position ID (GEN.72) as before. If a row's parent was not one the issuer saw, `insert_sector`
  runs `assign_uids` for the rest as it did. The IDs are identical (14,000 rows across 30 sectors
  compared with the SELECT-and-UPDATE pass; a test clears every ID and checks `assign_uids` finds
  the same). The pass it replaces cost 0.026 s per sector at about 9 systems a sector, more in
  a denser one.
- **D. Build INSERT text without pymysql's per-value escape for numeric-only tables**
  (`bright_stars`, `phenomenon_scatter`): 13.8 to 10.1 us per row, 27%. Risk: floats and
  ints only plus a `kind` drawn from a fixed list; any string value must still go through
  the escaper. Low priority after A.

### Design changes

- **E. Stop storing every phenomenon row.** Two options. (1) Compact the row: `kind` as
  TINYINT, positions as cell-relative 32-bit integers (the ring/layer/slot columns already
  give the cell), velocities computed when the object is built, no per-row seed (derive it from
  the cell and the row's index). About 3 times smaller (to about 45 bytes) [R]. (2) Do not
  store the numerous kinds (neutron stars and black holes are 99.9% of the rows) at all:
  derive them per cell from the existing per-unit seeds (GEN.39, `util/draw.py`) when the
  sector is filled or queried, and keep rows only for the rare kinds. This removes the
  scatter's writes, the 161 GB and the index, at the price of a cell-level count the map and
  the settle step must get from the generator instead of the table. Needs Boss's decision;
  I recommend (2) with (1) as the fallback.
- **F. Move the nearest-system and containment work out of the per-sector save.** A single
  later pass (the GEN.126 settle step already runs after a plan) would take the neighbour
  locks off the critical path, which is what limits the fill to 1.87 times at 2 to 4
  workers [S, report]. Risk: queries made mid-fill see incomplete links, which they already
  can.
  **Built (PERF.45).** A `galaxy` run saves its sectors with `link_neighbors=False` and
  `run_galaxy.link_after_run` links them once at the end (`store.link_sector_neighbors`):
  containment, each object's nearest systems and quadrant, then the new systems merged into the
  lists of the already-linked sectors around, 200 sectors per transaction under the neighbour
  lock, retried on a deadlock. It runs after a cancelled or failed run too (for the sectors that
  were saved). The first sector of a molecular cloud still takes the lock while it saves, since
  it stores the cloud. Sectors made on demand (API, `ensure_sector_generated`) link as they are
  saved. Measured by decision rule "whichever is more efficient": 4 workers, 400 sectors in ring
  2000, 52.0 s linking as saved against 46.9 s linking at the end (about 10% faster; the
  link pass is serial, 150 sectors took 3 s). The two runs' nearest-system lists and
  containment, compared system by system, are identical; so are single-worker runs compared
  row by row. If a run dies without linking, the orbit update script (`planetgen.cli.orbits`)
  links everything again.
- **G. Cut the repeated coordinate work in planet and moon generation.**
  `SpatialPosition3D._sync` ran 83,725 times in four sectors and `update_orbital_position`
  15,481 times; set the position once after the orbit is final.
  **Built (PERF.46, position half).** `SpatialPosition3D` now marks its derived coordinates
  and sector address out of date when a body moves and works them out when first read, and
  `carry_anchors`/`carry_sector_center`/`carry_star_center` read the coordinates they keep
  straight from the stored truth. `place_system` asks only a planet for its place (its moons
  need it). Coordinates are identical (a hash of 60 seeded systems, 1,666 bodies, before and
  after) and `place_system` took 0.035 s instead of 0.049 s. The `finite_domain` wrapper
  ran 202,622 times (about 3 to 4% of the fill); an environment switch to skip it during bulk
  fills loses a safety net, so only if Boss accepts that.

### What not to spend time on

- `is_name_valid` itself: 24 us, a compiled regex over the 279 forbidden words saves 1.3 us [C].
- Statement count or round trips: already 1.3 to 2.2 per system in a normal sector, not the limit.
- `--debug`: free on the scatter.
- Dropping the address index during bulk insert: no measurable gain at this size (13.8 against 13.2 us) [C]; at 1e9 rows it may matter, see below.

## Boss's question: how to cut the database cost without losing anything

The calls are not chatty; what costs is **how many rows are written and how many times**.
In order of safety:
1. **A** removes the per-job start-up (scatter and every queued job).
2. **C** writes each row once instead of insert-then-update, removing 6 to 9% of a fill
   with the same data.
3. **B** and **G** are Python-side and remove another 8 to 12%.
4. **F** takes the neighbour work off the lock so more workers help.
5. **E** is the large one: the phenomenon table.
Nothing here drops a column or a row that the maps, search or API read today; E(2) changes
where the numerous rows come from, not what they are.

## What this does not predict

- **Index and buffer-pool behaviour at 1e9 rows.** Quarter-scale tables fit in memory; the
  per-row insert cost rises once the index no longer fits the InnoDB buffer pool, and 161 GB
  is far outside a 4-core box's pool. Expect the phenomenon scatter to run several times
  slower per row than 32 us. Measure on the real server with `innodb_buffer_pool_size`
  recorded.
- **Contention among workers.** Measured here with 4 workers on 4 cores that also host
  MariaDB; the fill's neighbour locks gave 1.87 times at 2 to 4 workers in the earlier report.
- **Python version and hardware.** Python 3.13 here, 3.9 minimum for the project; the server
  may be slower or faster.
- **The galaxy mix.** The sectors were picked toward bright-star sectors; an average sector
  has fewer moons per system.
- **Disk and commit latency.** Commits are about 0.05 to 0.1 per system here; a spinning disk
  or a durable `innodb_flush_log_at_trx_commit=1` on a slow volume moves that share.

## Where ten hours could come from

Counts at the default scale, from the quarter-scale run times 64 [C]:

| Phase | Rows or stars | Warm CPU | At 4 workers |
|---|---|---|---|
| Bright-star scatter | 26.9M | 39 CPU-min | about 10 min |
| Phenomenon scatter | 1.17e9 | about 10.4 CPU-h | about 2.6 h |
| Bright-star fill (46.6 ms per star) | 26.9M | 14.5 days single worker | about 7.7 days at 1.87 times |

The report's figure of 24 minutes for the scatter does not reach ten hours. The phenomenon
scatter could: 9.4 CPU-hours before the per-job tax, the cost of a growing 161 GB index and
the database sharing the cores. If 10 hours was `planetgen plan` with the phenomenon scatter, A and E are the
fix. The notes record 1.6e8 neutron-star and black-hole rows; my count at the default
scale is 1.17e9 rows in total, so the two disagree by 7.3 times and one of them comes
from a different configuration. If 10 hours was the bright-star **fill**, only a smaller
galaxy fits and the per-star fill cost (46.6 ms, of which 2.6 ms is generation) is the target
(C, F). The run's log (or `planetgen plan` output with the phase times) would settle it.

## Proposed items (the TODO thread assigns IDs)

1. Pre-import generation modules and warm `_bright_table` in `cli/worker.py` before the fork;
   restart workers on update or compare versions in the horse (A). Gives PERF.18 and PERF.39 their measured number.
2. Lazy word-salad names for phenomena named by object ID (B).
3. Compute and write uids with the row, drop the post-insert pass (C).
4. Decide the phenomenon rows: derive on demand or compact (E); needs Boss's choice.
5. Nearest-system and containment as one later pass (F).
6. Orbital position updates once per body; optional `finite_domain` skip in bulk fills (G).
7. Numeric-only INSERT formatter or C driver for bulk tables (D), low priority.
8. PERF.31 benchmark to include the three phases above and the worker start-up cost, and to
   record `innodb_buffer_pool_size` and the table sizes with each run so the quarter-scale
   numbers can be compared with a server run.

## Evidence notes

Scripts in `/mnt/project-files/research/scripts/genperf/`: `common.py` (load the skeleton),
`scatter_profile.py`, `scatter_split.py` (draw against insert per star), `galaxy_totals.py`
(expected systems, bright stars and phenomenon rows), `phenom_profile.py`, `insert_split.py`
(executemany against hand-formatted and index-free), `fill_profile.py` (cProfile, server
statement counts, per-statement time, salad time), `fill_sample.py` (pyinstrument). Output of the
statement table: `/mnt/project-files/research/genperf-results/fill_statements.txt`. The preload
was an experiment in a `sitecustomize.py` outside the repository; no repository file was changed.
The bright-star report's earlier numbers (24 min, 74.6 CPU-min, 46.6 ms per star) are [S] from
`/mnt/project-files/bright-star-timing/report.md` and were not re-measured.
