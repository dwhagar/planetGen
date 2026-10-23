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

## What `medium` and `large` show

See the section at the end of this note once those runs finish (the numbers are added with the run date).

## Findings sent to the TODO thread

The warm-up's `202` failure (1a), the deep-tile row reads (1b), the band scatter's 348,000 single-row updates (2), the
supervisor's serial expansion (4) and the map-luminosity update (5) are filed as items. The bright-star unit's
`tier_at` per candidate in the reach scatter (6) is filed with them.
