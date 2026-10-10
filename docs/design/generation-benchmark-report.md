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

## After PERF.79 (scatter layers in chunks)

The plan's three scatter passes now send their layers to the workers in at most 8 chunks per worker instead of one task per layer (`WorkQueue.submit_each`). Same benchmark (ring 12, seed ...BEEF, 2 sectors), machine otherwise idle:

| workers | plan | phenomena | massive | bright | objects |
|---|---|---|---|---|---|
| 1 | 38.5 s | 18.4 s | 8.8 s | 6.7 s | 52,517 |
| 2 | 24.2 s | 8.2 s | 6.0 s | 5.4 s | 52,517 |

Before the change two workers took 593.9 s for the plan (about 12 times slower than one). The object counts match between one and two workers, so the rows are the same. The scatter is still one job on the RQ queue; the chunks are tasks of that job.
