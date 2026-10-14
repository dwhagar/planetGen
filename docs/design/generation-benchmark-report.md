# Generation benchmark report (PERF.31)

Measured with `planetgen benchmark --workers 1,2 --profile` on the lane's
build box (4 cores, 16 GB, MariaDB 10.11 on the same box, Redis local).
Galaxy: seed `...BEEF`, plan to ring 12 (a small galaxy: 2,041 layers), then
`galaxy --ring 12 --limit 2`. Two runs of the 1-worker case agreed within 5%.

Re-run with `python -m planetgen.cli.generate benchmark --workers 1,2 --profile --output report.json`
(about 15 min). It makes and drops its own scratch databases.

## Where a 1-worker run goes (138 s)

| Stage | Seconds | Share |
|---|---|---|
| Plan: scatter the phenomena | 23.7 | 17% |
| Plan: scatter the massive stars | 10.6 | 8% |
| Plan: scatter the bright stars | 8.2 | 6% |
| Fill: generate the 2 sectors | 48.4 | 35% |
| Fill: save the sector paths (settle) | 33.4 | 24% |
| Fill: link neighbours | 4.5 | 3% |

The database ran at about 980 statements per second (311 inserts/s, 180 updates/s, 179 commits/s).

## Findings, largest first

1. **Two workers make the plan scatter 12 times slower** (plan 47 s -> 594 s;
   bright-star scatter 8 s -> 282 s, massive stars 11 s -> 206 s, phenomena 24 s -> 101 s).
   The server answered 542,000 statements against 135,000, so the extra
   cost is queue and database round trips per layer (2,041 small layers),
   not computation. The fill itself does speed up with two workers
   (generate 48 s -> 30 s). Cause not yet isolated; the poll constants in
   `queue/work.py` and a task per layer are the places to look.
2. **Save the sector paths costs 33 s for 2 sectors** (24% of the run and
   flat with workers). In the profile `settle_after_run` is 80 s of 227 s
   and `integrate_path` (17,284 calls) is 54 s.
3. **Generating a sector**: `nearest_neighbors` makes 1.8 million
   `distance_to` calls for 1,916 systems (41 s of 227 s in the profile); a
   spatial grid would cut this by a large factor.
4. **Inserting a sector** is 40 s of the profile (`insert_sector`), 17% of it.
5. The scatter stages scale with layers, not objects: bright-star scatter
   visits 2,041 layers and places 0 objects (8 s).

The profiler slows the run 2.5 times, so read its seconds as shares of the
run, not as wall time.

## Not measured

Default-scale galaxy; MySQL 8.4; more than 2 workers. The profile is of
the fill only, not the plan.


## Save the sector paths, after PERF.80

Measured by re-saving the paths of the benchmark galaxy's two sectors
(`sector_paths.settle_sectors`, 17,284 bodies: 1,916 star systems and the
rogue planets and comets; the two sectors and their neighbours) on the
kept scratch database, 4 cores, MariaDB 10.11, an otherwise idle machine:

| | Seconds |
|---|---|
| Before (main) | 44.7 |
| After, 1 worker | 20.0 |
| After, 2 workers | 13.3 |
| After, 4 workers | 10.8 |

The saved paths are the same to the last bit: a hash over every path and knot
of the 17,284 bodies is identical before and after, at every worker count.

What changed:

- **The cell test.** `integrate_path` asked `geometry.sector_address_at` (a
  ring, a layer and a slot, each worked out afresh, `ring_sector_count` three
  times over) at every step and 40 more times to bisect where the path leaves
  the cell: 43% of its time. `geometry.sector_cell_test` repeats the same
  arithmetic with the cell's numbers worked out once and stops at the first of
  ring, layer and slot that is not the cell's; a test holds it equal to
  `sector_address_at == address` on thousands of points including the cell's
  faces. `ring_sector_count` and `ring_master_count` are cached.
- **Smaller costs.** The bisection builds only the position of the Hermite
  curve it tests; the integrator unpacks its masses to `(x, y, z, mu)` once
  per path instead of computing `G * mass` at every step for every mass.
- **One statement batch per sector.** A body's path took a lookup, a delete, an
  insert and a knot insert; a sector's paths are now saved with one lookup per
  500 bodies, one insert of the path rows and one of the knots.
- **Worker processes.** `galaxy` passes its worker count to `settle_sectors`, which hands
  a crowded sector's bodies (250 per task) to that many spawned processes;
  loading, saving and the commit stay in the calling process, and the paths
  and their order do not depend on the count. What remains at 4 workers
  (about 10 s) is that serial part: the sector and neighbour query and the
  inserts.

## After PERF.79 (scatter layers in chunks)

The plan's three scatter passes now send their layers to the workers in at most 8 chunks per worker instead of one task per layer (`WorkQueue.submit_each`). Same benchmark (ring 12, seed ...BEEF, 2 sectors), machine otherwise idle:

| workers | plan | phenomena | massive | bright | objects |
|---|---|---|---|---|---|
| 1 | 38.5 s | 18.4 s | 8.8 s | 6.7 s | 52,517 |
| 2 | 24.2 s | 8.2 s | 6.0 s | 5.4 s | 52,517 |

Before the change two workers took 593.9 s for the plan (about 12 times slower than one). The object counts match between one and two workers, so the rows are the same. The scatter is still one job on the RQ queue; the chunks are tasks of that job.

## Writing the phenomenon rows (PERF.73)

At the default scale the phenomena pass writes about 68 million rows, so the write is most of its time. Measured on 1,000,000 synthetic rows of the real shape (40 layers of 25,000, same machine, MariaDB 10.11, default buffer pool, one connection, commit per statement batch):

| way of writing | rows per second | notes |
|---|---|---|
| driver `executemany`, 10,000 per statement (before) | 37,000 | |
| same, 1,000 or 50,000 per statement, or one commit per layer | 37,000 to 40,000 | batch and commit size do not matter |
| same, `unique_checks`/`foreign_key_checks` off | 37,000 | the table has neither |
| multi-row `INSERT` text built in the store instead of the driver's per-value escaping | 49,000 | the driver's escaping alone is about 9 s of the 27 s |
| built text and the two secondary indexes dropped, then built once at the end | 60,000 (79,000 for the write, 4 s to build the indexes) | **what the scatter does now** |
| `LOAD DATA LOCAL INFILE`, indexes dropped | 89,000 | needs `local_infile` on the server and the client; not used |

So the scatter now writes 1.6 times faster than before on this table, and the 68 million rows go from about 31 minutes of database time to about 19 plus the index build. The rows are the same (a test stores numpy numbers, `NULL`s and quoted names both ways). A scatter that dies half way leaves the indexes off but also no recorded seed, so no sector fill uses the table until the next scatter clears and rebuilds it. `LOAD DATA` stays an option if the host's server allows it.

## PERF.87 and PERF.90: placement as arrays

`SpaceSector._random_position` tests a candidate against all massive neighbours
with one NumPy comparison; points within 1e-9 of a boundary fall back to
`math.dist`, so accept and reject decisions match the old loop. Two tests prove
it (a fill compared with the old loop, and a borderline check). Benchmark,
2 sectors at ring 12, 1 worker, idle machine: the sector fill went from
68.0 s to 64.0 s (the plan stages are untouched and moved by noise, 40 to 45 s).
`_ring_centerline_edges` joins `ring_master_count` and `ring_sector_count` in
the cache.

## PERF.89: population densities as arrays

`density.population_densities_array` gives the four population densities for
every angle bin of a ring in one call; `bright_stars._ring_bins` and
`layer_expected_stars` use it. NumPy's `exp`, `log`, `cos` and `arctan2` are not
bit-for-bit the maths library's: about 12% of values differ in the last digit,
so the tests hold the array form to 1e-12 relative of the one-point form rather
than to the last digit (Math Research asked for the last digit; flagged). A
bin weight moving by 1e-16 can change a star's placement only if a random draw
lands within 1e-16 of a bin edge. Benchmark bright-star stage, 2 sectors at
ring 12: 8.5 s to 7.5 s (the other plan stages moved by up to 5 s between runs,
so this is within noise at benchmark scale; the default-scale galaxy has far
more bins).

## DB.24: the nebula field table

`planetgen plan` writes `nebula_field` inside the phenomena stage (one task per
first-index slab of 50 pc field cells, through the work queue). Benchmark
(2 sectors, ring 12): the phenomena stage stays at 18.8 s on one worker and
12.3 s on two, as before the table, so the build is small at that scale. The
table has no object IDs for clouds no sector has stored: working an ID out costs
a nebula shape (`interior_centroid` over the fit grid, about 112 ms), which for
every drawn cloud would be about 30 hours of CPU at the default scale. A save
sets `uid` and `stored` on the row of each cloud it stores.
