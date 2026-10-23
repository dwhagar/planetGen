# Generation performance: every step timed, and what to speed up next

**Status:** reference (measured at release 8.0.1930 for PERF.146; the numbers are repeatable with `planetgen benchmark --suite`).

The earlier [generation-performance-study.md](generation-performance-study.md) profiled three phases on a scaled-down galaxy.
This note covers every generation step of a whole galaxy, from the plan to the cache warm-up, at four galaxy sizes,
with the suite that measures them, the fixes that came out of the first round, and a ranked list of what is left.

## How to run it

```sh
planetgen benchmark --suite --size small --workers 3 --step-timeout 90 --output suite.json
```

It builds a scratch galaxy (`planetgen_suite_*`, dropped at the end), runs each step as its own child command, and
measures the child's wall time, the machine's CPU (from `/proc/stat`, less the database and Redis servers), the MySQL
server's CPU, the peak memory of the `planetgen` processes, the MySQL counters (questions, rows read, inserted and
updated), the stages the step reports, and what it made, so the time per star or per cell can be read off. `--repeats N`
runs the whole suite N times and prints the fastest, median and slowest run of each step. `--steps a,b` runs only some
(a step needs the ones before it to have run). `--step-timeout MINUTES` stops a step that runs away and records it as
not finished. `docs/cli.md` has the options. A test (`src/tests/test_benchmark_suite.py`) runs the first step on the
`tiny` size and checks the JSON's shape.

The sizes grow by the galaxy's mass (`--mass`): `tiny` 2e6, `small` 3e7, `medium` 4e9 (the size a seed draws), `large`
5e10 solar masses. **Even `tiny` is not tiny**: the galaxy's 157 globular clusters and 11 satellites exist at every
mass, and the 6.4e9 systems of the satellites make millions of phenomena and hundreds of thousands of band stars whatever
the galaxy's mass. The tiny size is for checking that the suite runs, not for timing.

The numbers below are from a 4-CPU machine with a 128 MB InnoDB buffer pool, MySQL 8.0.46, three workers. They are
about the shape, not the seconds: another box will differ, and the run-to-run spread of one step is 4% to 15%.

## The profile recipe

`py-spy` (`pip install py-spy`; it is not a project requirement) attaches to a running worker; the horses the worker forks are its
children, so `--subprocesses` follows each unit:

```sh
W=$(pgrep -f 'planetgen.cli.worker' | head -1)
py-spy record --pid "$W" --subprocesses --duration 60 --format raw -o worker.raw
```

The supervisor (the `planetgen scatter ...` process) is the one that expands units serially; profile it with
`--pid` alone. The raw output is `stack count` lines; summing the counts by frame, once as inclusive (every frame of
the stack) and once as self (the last frame), gives the tables used here. Profile a step only while nothing else runs,
and read the step's **steady state**, not its first minute (a worker warms its caches before its first fork).

## What a full run costs (`small`, 3e7 solar masses, release 8.0.1930)

Total 25.7 minutes. One run, the final code; the warm-up failed (see the ranked list), so its time is a lower bound.

| Step | Wall | Share | Client CPU | Server CPU | Peak memory | Questions | Rows read | Rows written | Made | Per unit |
|---|---|---|---|---|---|---|---|---|---|---|
| Plan the galaxy | 17.7 s | 1.2% | 20 s | 1 s | 485 MB | 264 | 11 | 721 | 327 components | 54.25 ms |
| Scatter the massive stars | 2.0 min | 7.9% | 157 s | 57 s | 1188 MB | 78,080 | 428,151 | 164,735 | 24,558 stars | 4.95 ms |
| Scatter the bright stars | 2.8 min | 10.7% | 288 s | 55 s | 1341 MB | 76,219 | 425,331 | 80,975 | 14,037 stars | 11.78 ms |
| Scatter the phenomena | 4.3 min | 16.9% | 407 s | 234 s | 1267 MB | 109,609 | 3,401,724 | 3,313,979 | 3,108,512 phenomena | 0.08 ms |
| Draw the nebula field | 39.8 s | 2.6% | 95 s | 10 s | 1179 MB | 11,306 | 73,113 | 29,158 | 1,338 clouds | 29.78 ms |
| Add the dimmer band of bright stars | 5.2 min | 20.1% | 573 s | 205 s | 1320 MB | 422,282 | 3,930,516 | 2,397,360 | 486,761 stars | 0.64 ms |
| Make the nebula block totals | 4.3 s | 0.3% | 4 s | 1 s | 253 MB | 163 | 29,800 | 27,801 | 24,294 blocks | 0.18 ms |
| Fill sectors (6 sectors) | 46.4 s | 3.0% | 117 s | 3 s | 1362 MB | 3,202 | 24,155 | 1,619 | 6 systems | 7,735 ms |
| Run the population pass | 4.4 s | 0.3% | 5 s | 0 s | 169 MB | 140 | 105 | 8 | 0 polities | |
| Re-score the habitability | 14.3 s | 0.9% | 15 s | 1 s | 455 MB | 529 | 511 | 35 | 106 planets and moons | 135 ms |
| Fill the map luminosity | 75.7 s | 4.9% | 45 s | 162 s | 1034 MB | 1,549 | 6,217,120 | 3,094,263 | 3,108,512 phenomena | 0.02 ms |
| Warm the caches (failed) | 8.0 min | 31.3% | 120 s | 395 s | 523 MB | 42,920 | 95,888,396 | 134 | 6 sectors | |

Worker start-up is 1.5 s a process (1.1 s of imports, 0.1 s for the first connection). Every step also pays a fixed
2 to 3 s for the command's own start.

**Run-to-run spread** (`tiny`, three runs): plan 4%, massive stars 15%, bright stars 11%, phenomena 10%, nebula field
8%, band 10%, nebula block totals 10%, fill 14%, population 13%, habitability 4%, map luminosity 9%. The warm-up did not
finish in any run, so it has no spread.

## What the first round fixed

All of these leave every drawn value as it was (the generator epoch is untouched); each has a test that compares the
new path with the old one bit for bit.

| Fix | Before | After |
|---|---|---|
| A queue worker builds the galaxy's component field and the sampling tables **once, before it forks a job**. The horse is a fresh process per unit, so each unit was rebuilding the field from about 22,000 stored clusters and satellites. | massive-star scatter at the default size over 2 hours (stopped) | about 3 minutes |
| The mean system mass at an age (`mean_system_mass_at_age_solar`) works out the remnant states and the turn-off mass once, and the lifetime terms once per metallicity. | 1.0 | 4.4 times faster |
| The scatter's band tiers are looked up **together**: one pass over the density arrays for every place a component's stars land in, and one for the marked cells of a shell, instead of one band at a time. The longitudes use Python's `math` trig, so the floats are the loop's own. | 86% of a bright-star unit | the bright-star step 5.2 min to 2.8 min at `small` |
| The giant branch's metallicity terms (`a27 .. a32`, `b9 .. b13`) are worked out once per metallicity. | 0.36 s per 20,000 branches | 0.07 s |

Whole `small` suite: 30.4 minutes with the first two fixes, 25.7 with all four. The spread of one run is 4% to 15%, so
smaller differences are not shown by one run.

## Ranked: what to speed up next

Each cost's share is of the `small` run above. "Value-neutral" means the drawn values stay as they are; a change that
alters a drawn value needs its own item, a generator epoch bump and Boss's word, and goes to a later version.

1. **Warm the caches, 31% and it fails.** (a) The deep-tiles phase builds 1,044 tiles in 320 s and reads 95.9 million
   rows doing it (the whole step, about 92,000 rows a tile): a tile reads far more than it draws. (b) The search check
   (`/search`) answers `202` (the page is slow enough to be queued as a job, see
   [slow-page-answers.md](slow-page-answers.md)) on a galaxy of 3.1 million phenomena, and the warm-up counts a `202`
   as a failure, so the step stops after its fourth phase. Either the warm-up should wait for the job or the search
   should not need one. Value-neutral. A bug and a cost, filed as both.
2. **The band scatter, 20%.** Steady state, per unit: about 53% drawing each star's member tree
   (`scatter_companions.star_trees` into `multistar.place`, the stability checks), about 19% writing it
   (`executemany`), the rest drawing the stars. The tree draw is the generator itself, so a change to it changes values;
   the writes are not: the step makes 348,000 `UPDATE` statements (one per row) that could be batched, and the row
   writes could go in larger groups. The batching is value-neutral; its saving is not measured.
3. **The phenomenon scatter, 17%, for 3.1 million rows in a 3e7 galaxy.** 3.3 million rows written, 234 s of server CPU
   (the database is the cost, not the Python). Ideas: larger insert batches, fewer secondary indexes while the rows go
   in (add them after), and the question of why a dwarf galaxy holds 3.1 million phenomena (the satellites'
   6.4e9 systems at every galaxy mass: a modelling question, not a speed one).
4. **The scatter supervisors, 8% + 11%.** The massive and bright-star steps use 1.3 and 1.7 of the machine's 4 cores:
   the supervisor expands the 65,536 shells serially, so the workers idle. Value-neutral; expanding in the workers or
   in larger groups would use the cores.
5. **Fill the map luminosity, 5%.** 3.1 million single-row updates and 6.2 million rows read, 162 s of server CPU
   against 45 s of client CPU. One set-based `UPDATE ... JOIN` per block would do it. Value-neutral.
6. **Tier lookups in the far-reach scatter** (`bright_stars._reach_spots`) still ask `tier_at` per candidate; the
   component pass was batched in this round. Small at `small`, grows with the reach.
7. **The nebula field, 2.6%,** for about 1,400 clouds (30 ms each), and **`clear-scatter --nebula-field`** (11 to 15 s):
   both small; the suite runs the clear untimed.
8. **The plan, 1.2%** (18 s): the field build and storing the 22,000 components. Fine at this size; check at `large`.
9. **The fill, 3%** at 7.7 seconds a system: too few systems at this size to say anything; the fill needs a size of
   its own to be measured (the next version of the suite should fill a fixed count of sectors).

## The plan stage at the default size (PERF.163)

Boss's plan of a default galaxy (4e9 solar masses, 22,267 clusters and satellites, 1,089 shells) took 264 s against a
3.7 s skeleton at release 8.0.1961. The plan log now ends with a table of where the unit's time went
(`util/step_timer.py`), costliest first. On the 4-CPU test box, release 8.0.1960 (before PERF.146) took about 2 minutes
for the whole `plan` command, the release after PERF.146 about 30 s, and the release with this table about 25 s (three
runs each, on a machine that is not Boss's, so the seconds differ and the ratios are the finding). The table of the last run, 12.3 s
of the unit:

| Step | Seconds | Share |
|---|---|---|
| The component field (the mean system mass of each of 22,267 clusters and satellites) | 4.5 | 36.7% |
| Storing the clusters and satellites (22,267 rows) | 3.1 | 25.3% |
| Storing the 1,089 shells | 1.1 | 8.8% |
| The report lines (pins, halo, outline) | 1.0 | 8.4% |
| The young open clusters (22,099 drawn) | 0.7 | 5.6% |
| The outline (the skeleton's shells; the "3.7 s" of Boss's log) | 0.7 | 5.6% |
| The galaxy's shape | 0.6 | 5.1% |
| The globular clusters, the satellites, the rest | 0.5 | 4% |

The rest of the 25 s is outside the unit: the command's start-up and argument checks (about 3 s, 0.6 s of it
`plan_shape` run once for the check and once for the unit), the start-up math check (2 s), and starting the queue and
its workers (about 6 s).

What was wrong: the cache of `mean_system_mass_at_age_solar` held 512 ages while a plan asks for 22,267 distinct ones (every
open cluster rounds to its own age), so the log line that prints the clusters' system counts paid the whole bill a
second time (4.4 s here, about 3 times that on Boss's box). It now holds 65,536 ages; the search for the mass of a
given lifetime also reads the metallicity's terms once instead of once per bisection step. Neither changes a value
(tests compare with the plain loop bit for bit). Together with the 4.4 times faster mean system mass from PERF.146 this
brings the plan to about a fifth of its time at 8.0.1960.

What is left, ranked: (1) the mean system mass of each age costs 0.2 ms, mostly the 50-step bisection for the turn-off
mass (`youngest_progenitor_sol`): a table of the turn-off mass by age, or a root finder with fewer steps, would change
the last bits of the answer, so it needs its own item and a generator epoch; (2) storing the 22,267 components takes 3 s,
about 0.14 ms a row; (3) the plan's start (the `plan_shape` check run twice, the math check, the queue's workers) is
about 10 s of the 25 s. GEN.624 will cut the satellites' systems; the component field's cost depends on the number of
components, not their size, so it will not shrink.

## What `medium` and `large` show (PERF.164)

Same 4-CPU box, three workers, release 8.0.1973, seed `...BEEF`, one run each. The container restarted during the
medium run, so it ran in three pieces over the same scratch databases (`--resume`); the steps that were run twice
differ by 10% to 30% between the two runs (the plan 29 s and 35 s, the band scatter 10.5 and 14.3 min), which is the
spread to read these by. The medium numbers below are from the complete pass of each step.

| Step | `small` (3e7) | `medium` (4e9) | `large` (5e10) |
|---|---|---|---|
| Plan | 17.7 s, 327 components | 29 s, 22,267 components | 43 s, 9,306 components |
| Massive stars | 2.0 min, 24,558 stars | 3.6 min, 68,499 stars | 4.0 min, 69,407 stars |
| Bright stars | 2.8 min, 14,037 stars | 8.8 min, 899,474 stars | 11.2 min, 1,307,239 stars |
| Phenomena | 4.3 min, 3.1 million | 14.1 min, 6.7 million | over 60 min: stopped at the cap with 66.0 million rows stored |
| Nebula field | 40 s, 1,338 clouds | 4.0 min, 8,166 clouds | not run |
| Dimmer band | 5.2 min, 486,761 stars | 10.5 min, 800,953 stars | not run |
| Block totals | 4.3 s | 18 s, 153,815 blocks | not run |
| Fill (6 sectors at small, 12 at medium) | 46 s | 75 s, 36 systems | not run |
| Population, re-score | 19 s | 15 s | not run |
| Fill the map luminosity | 76 s | 2.6 min | not run |
| Warm the caches | 8.0 min, fails | 27.8 min, fails | not run |
| **Total without the warm-up** | **17.7 min** | **46 min** | |
| **Total with the warm-up** | 25.7 min | 74 min | |

A galaxy 130 times heavier (small to medium) takes under 3 times as long: the cost is the count of things drawn
(components, shells, phenomena), and the clusters and satellites exist at every mass. Medium to large (12 times the
mass) adds 11% to 48% to the first three steps (the plan 29 to 43 s, the massive stars 3.6 to 4.0 min, the bright
stars 8.8 to 11.2 min). The per-star times therefore do not carry from one size to another: the
massive-star scatter costs 4.95 ms a star at small and 3.1 ms at medium, the band 0.64 ms and 0.79 ms. A real
default-size run is the `medium` column; its steps in a real run are the same commands, so these times are what a
real run takes on this box, plus the whole sector fill (the benchmark fills only 6 or 12 sectors).

The warm-up is the largest single cost at both sizes and still fails: at medium its deep-tiles phase builds 1,665
tiles in 948 s, and the `/search` check answers `202` again (see the ranked list, item 1). At large the phenomenon scatter
was still running when the 60-minute cap stopped it, with 66.0 million rows stored (ten times the medium galaxy's 6.7
million by the time it was stopped, and not done), so the steps after it were not run: the large galaxy's run is at least
1.3 hours on this box (16 minutes for the first three steps, then the phenomena), with the phenomena the first thing to speed up (ranked list, item 3). The cap is the
suite's `--step-timeout`, 240 minutes by default; the run used 60.

## Findings sent to the TODO thread

The warm-up's `202` failure (1a), the deep-tile row reads (1b), the band scatter's 348,000 single-row updates (2), the
supervisor's serial expansion (4) and the map-luminosity update (5) are filed as items. The bright-star unit's
`tier_at` per candidate in the reach scatter (6) is filed with them.
