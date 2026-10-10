# Keeping slow database reads inside the time limit (PERF.71)

Boss (2026-10-10 20:50Z): "is there a TODO item to maybe batch or give partial responses to queries that take too long? ... research how to optimize DB calls to large data sets to avoid timeouts and give systems a chance to respond. We could put recall passes into the queue too if they take longer than 10 seconds, then keep the query open until it finishes? Or split it up? I'm not sure, we need research."

Informs: PERF.68 (Galaxy Map tile queries), PERF.69 (Planets, Moons and Phenomena counts), PERF.70 (Systems sort by Sector or Octant), and the new items in section 7. Builds on PERF.19 and PERF.24 (the queue), PERF.34 and PERF.38 (tile cache under a fill), PERF.36 and PERF.64 (candidate-cell guard, stored counts, the busy page).

Status: research, 2026-10-10; nothing here is built. Evidence tags: [S] seen in the repo, [C] measured or computed in this research, [R] recalled and not confirmed (listed at the end).

## Summary

- **Do not make "queue it and keep the request open" the general answer.** It does not make the query faster (it runs the same SQL on the same busy database), it holds one of the API daemon's 5 threads for the whole wait, a queued job needs 2.4 s just to start a worker, and the polling route is admin-only [S, C: sections 2 and 3]. Keep it for the few operations that are meant to run long and are started by an admin (estimates, the deep database check, exports).
- **Every page that timed out in these measurements is slow because it reads far more rows than it shows, not because the work is irreducibly big.** Each one has a cheap form: a keyset page instead of `OFFSET` (8.5 s to 2 ms), a stored per-sector count (18.6 s to 1 ms), an index for the map's scattered points (31.7 s to 0.28 s), a capped count (6.0 s to 30 ms) [C].
- **So the order of preference is: read fewer rows, then answer in pieces, then show what is ready, and only last queue.** That fits what the code already does (stored counts, single-flight tiles, the busy page) and extends it.
- **A fill makes every read 1.2 to 2 times slower** on the test box (a page at 0.13 s took 1.08 s; a sector list at 17 s took 29 s) [C]. A query that needs 6 s idle is a timeout during a fill. The target for any page query is therefore about 1 s idle, which every proposed form meets.
- **Partial answers are worth building for three places:** the table pages (rows first, filter counts after), the search page (each result panel on its own time limit) and the Galaxy Map tiles (a tile without its slowest piece, filled in later). Section 5 says what the visitor sees.
- Eleven build items for the TODO thread are in section 7. Three questions for Boss with defaults are in section 8.

## 1. What happens today

[S] unless marked.

| Piece | Behaviour |
|---|---|
| Statement time limit | `mysql.statement_timeout_seconds`, default 10. The web interface's and API's read-only connections set MariaDB `max_statement_time`; a stopped statement raises error 1969 (3024 on MySQL). |
| Page that hits it | 504 page, "The database is busy", reloads itself every 10 s, `Retry-After: 10` (PERF.64). The API answers `QUERY_TIMEOUT` with `Retry-After`. The Galaxy Map asks again after 2, 4, 8 ... 30 s. |
| Table pages | `static/datatable.js` shows a virtual list of 50-row pages and fetches `/table/<name>?offset=&limit=50` as the visitor scrolls. The first page also asks for the filter-menu counts (`facets=1`) in the same request. `web/tables.py` turns any API error, a 504 included, into a 502 `{"error"}` with no `Retry-After`, and `datatable.js` throws on a non-200 answer without retrying. |
| Stored counts | `db/countcache.py` (PERF.64): totals and filter-menu counts for the Systems and Sectors tables are served from the last count; one background thread makes the next. Planets, Moons and Phenomena still count in the request (PERF.69). |
| Tiles | `web/lib/tilecache.py` builds missing tiles one request at a time (Redis claim, `singleflight.py`) and keeps the old generation while a fill runs. |
| Queue | RQ on Redis. `api_jobs.submit` starts a detached burst worker per waiting job (up to `worker_count()`); `run_queued` waits 8 s and answers `202` with a job id; `GET /api/jobs/<id>` gives the state. The route needs `require_admin(scope="read")`. Results are kept a day. |
| Web serving | Apache `WSGIDaemonProcess planetgen-api processes=1 threads=5 request-timeout=60`. |
| Earlier measurements | A queued job costs 2.4 s before it starts (median of 5; the import of `planetgen.api.common` is 2.3 s of it); a cache hit costs microseconds; the earlier design concluded "routing reads through RQ is slower than computing them" and put the queue only behind cache misses of work long enough to need it ([performance-eta-queue-and-caching.md](performance-eta-queue-and-caching.md) 3.2 and 4.4). The reserved `planetgen-interactive` worker proposed in its 3.4 was not built. |

## 2. Measurements

Setup [C]: MariaDB 10.11 on a 4-core, 16 GB box, InnoDB buffer pool 1 GB. Database `gp3`: 2,000,000 synthetic single-star systems over 150,000 synthetic sectors, cloned from one real system (about 3 GB of data and indexes, so it does not fit the pool, as a real galaxy would not). `gp2`: the earlier quarter-scale galaxy (17,793,451 `phenomenon_scatter` rows, 415,084 bright stars), migrated to schema 81. Each case is the median of 3 runs of the project's own query function (`query.list_systems` and the others), no stored-count cache, no statement limit. "Loaded" is a fill stand-in: 3 threads inserting committed 40-row batches (about 18,000 rows/s) and 4 CPU burners at nice 10, the way a fill competes for the box. Scripts and raw output: `/mnt/project-files/research/db-timeouts/scripts/` (`make_synthetic.py`, `bench_pages.py`, `tile_parts.py`, `fill_load.py`, `alternatives.py`) and `alternatives.log`.

### 2.1 Pages as they are

| Page query | Idle | Loaded | Why |
|---|---|---|---|
| Systems, total (no filter) | 0.25 s | 0.46 s | index scan |
| Systems, total with a star-type filter | 6.0 s | not run | joins `stars` |
| Systems, page 1 by name | 0.004 s | 0.016 s | index |
| Systems, rows 20,000 to 20,050 by name | 0.13 s | 1.08 s | `OFFSET` reads and drops 20,000 rows |
| Systems, offset 500,000 | 5.8 s | 9.8 s | same, 500,000 rows |
| Systems, offset 1,500,000 | 12.1 s | 15.6 s | same, 1,500,000 rows (7.3 to 8.5 s on warm repeats) |
| Systems, page 1 sorted by Sector (PERF.70) | 4.8 s | 6.3 s | sorts the whole join |
| Systems, page 1 sorted by Octant (PERF.70) | 0.9 to 4.1 s | 4.8 s | sorts 2,000,000 rows (`quadrant IS NULL` first) |
| Sectors, any sort, page 1 | 17.2 s | 29.0 s | counts all 2,000,000 systems per request to give each sector its count, then sorts |
| Systems search, a rare name | 2.8 s | not run | `search()` panel |
| Systems search, a common word (333,000 matches) | 3.7 s | not run | same |
| Systems search, a 2-letter word | 8.4 s | not run | same |
| Phenomena, total / page 1 / offset 1,000,000 | 0.08 / 0.15 / 0.42 s | not run | `phenomenon_scatter_classes` and the union view (fast) |
| Galaxy Map tile, level 2 / 5 / 8 (centre) | 33.5 / 32.8 / 7.4 s | 0.34 s at level 8 *after the index below* | the scattered-points query reads 8.6 million rows (EXPLAIN: range on `idx_phenomenon_scatter_address`, filesort) |
| Galaxy Map tile, level 10 / 12 | 0.11 / 0.01 s | not run | narrow address band |
| Same tiles, bright-star part only | 0.21 to 0.24 s | not run | already index-served |

The window to the limit is narrow: a page at 4 to 6 s idle is a timeout under a fill. The Sectors list is already over the limit at 2,000,000 systems idle. These grow linearly with the table: 10,000,000 systems means about 85 s for the Sectors list and 60 s at offset 1,500,000 [C by scaling, not run].

### 2.2 The alternatives

| Alternative | Before | After | Cost |
|---|---|---|---|
| Keyset page (`WHERE name > ? OR (name = ? AND id > ?)`), Systems by name, after row 1,500,000 | 8.5 s (`OFFSET`) | 0.002 s | Cannot jump to row N, only forward or to a name or letter ("jump to F" is 0.001 s). **The row-constructor form `(name, id) > (?, ?)` took 5.0 s on MariaDB 10.11**: it does not use the index range. |
| Capped count (`SELECT COUNT(*) FROM (SELECT 1 ... LIMIT 10001)`), star-type filter | 6.0 s | 0.030 s | The answer is "10,000 or more" above the cap. |
| Estimate from table statistics | 0.22 s (exact) | 0.000 s | Within tens of percent; already the PERF.64 fallback for an unfiltered table. |
| Stored per-sector system count (summary table `sector_id, n`, index on `n`) | 17 to 18.6 s | 0.001 s (by name), 0.23 s (by systems) | A refresh is one `GROUP BY` on the sector index: 0.42 s to 0.78 s for 2,000,000 systems, so a background refresh like `countcache` is cheap; or maintain it where systems are written. |
| Stored sector-name sort column on systems, index `(sector_sort, name, id)` | 4.8 s | 0.000 s (1.06 s with the `IS NULL` term) | One-off fill of 2,000,000 rows took 202 s as an `UPDATE ... JOIN`; index build 5.9 s; must be kept in step with sector renames. |
| Index `(quadrant, name, id)` for the Octant sort | 0.9 to 4.1 s | 0.000 s (1.07 s with the `IS NULL` term) | 6 s to build; the `ORDER BY quadrant IS NULL, ...` form defeats the index, so the order must be written to match it. |
| Index `(kind, subtype, mass_solar)` on `phenomenon_scatter` for the coarse tile points | 31.7 s at level 5 | 0.28 s (level 2: 0.05 s; level 8: 0.26 s) | 29 s to build on 17.8 million rows, about 0.6 GB. This scratch galaxy predates the scatter's subtypes, so a current galaxy returns more rows than this one; the plan (an index range on the few big classes, not an address range over millions) is what matters. |
| Statement limit | a 17 s query | stopped at 10.001 s with error 1969 | The limit is reliable and prompt; what to do after it is the open question. |

### 2.3 Queued read, from earlier measurements [S]

Starting the worker costs 2.4 s (median), mostly imports; a warm worker adds about 0.03 s; `status()` is 2.9 ms; the poll loop spends one of the API's 5 threads for the whole wait. Nothing in these measurements is specific to a read, and none was re-run here.

## 3. The four options, weighed

### 3.1 Queue a read that passes 10 s, keep the request open or poll

For: the visitor gets an answer instead of an error for work that cannot be made cheap; the machinery exists (`run_queued`, `GET /api/jobs/<id>`).

Against, from the numbers above:

- The query still runs on the same database, which a fill has made 1.2 to 2 times slower. Moving it to a worker only moves the wait.
- **Keeping the request open costs an API thread each.** The daemon has 5. Five visitors waiting on queued Sectors lists would leave the API unable to answer, and `request-timeout=60` kills the daemon if they wait longer than a minute.
- **Polling needs a public route.** `/api/jobs/<id>` is admin-only; ordinary visitors would need their own, with its own rate limit and a way to keep one visitor's job id from reading another's results.
- Each job pays 2.4 s to start a worker, and no reserved worker exists; during a fill a read job also competes with the fill's own workers for cores.
- A result goes through Redis (small pages are fine; a tile is a megabyte and costs 3.3 ms to read back, 4.3.3 of the earlier note).
- It hides a query that should not be this slow. The three worst cases above fall to a few milliseconds with an index or a stored count.

Verdict: **not for page reads.** Use it for operations an admin starts and expects to take long (Generate estimates, the deep database check DB.21, exports): they already show a progress or "working" state and already live behind admin routes.

### 3.2 Split one query into pages or key ranges

Already the model for tables (50-row pages). What is wrong is the *kind* of paging: `OFFSET` makes page *n* cost *n* times as much (0.004 s, 0.13 s, 5.8 s, 12.1 s at rows 0, 20,000, 500,000 and 1,500,000). Keyset paging makes every page cost the same (2 ms). The price is that the scroll bar cannot jump straight to row 1,500,000 by number; it jumps by value (a letter, a sector) or lands on an estimate and fills forward. For a list sorted by name that is what a reader wants anyway.

Tile pieces: a tile is eight independent queries (placed sectors, planned slots, filled cells, clouds, bright stars, generated stars, points, scattered points). Time per piece shows one piece (scattered points) holds 31.4 of 31.6 s at level 2. Splitting by piece, each under its own budget, is cheap to build and gives partial answers (section 3.3).

Verdict: **yes**, for table paging (keyset) and tiles (per-piece budgets).

### 3.3 Partial or streamed answers

- **Tables.** Rows and counts are separate costs: the stored count cache already separates them for totals, but the first page still asks for the filter-menu counts in the same request. Fetch the rows first; fetch `facets=1` in a second request; paint the rows at once and fill the menus when they arrive.
- **Search.** `search()` builds one panel per object type (sector, system, star, planet, moon), each with its own query. Today one slow panel times out the whole page. Run each panel under its own statement limit and fetch panels separately; a panel that times out says "still looking, try again" and the others show.
- **Tiles.** Give each piece a budget (say 3 s), build the tile with whatever finished, mark it `incomplete` in the tile cache so it is rebuilt after a short wait and is not served as final, and let the existing retry fetch the missing piece. A map tile missing its black-hole points for a few seconds is better than the "Took too long" page Boss saw.
- **Streaming** (chunked JSON, server-sent events) is not worth it here: every case above is "one page of 50 rows" or "one tile", small enough to return whole once cheap.

Verdict: **yes** for the three places above.

### 3.4 Stored counts, summary tables and indexes

This is where almost all of the win is. The cases:

- A *count* of everything: stored (done for Systems and Sectors, PERF.69 for the rest).
- A *count with a filter nobody has used*: capped count, "10,000 or more" (30 ms against 6 s).
- A *per-row aggregate* shown in a list (systems per sector): a stored column or summary table, refreshed in the background, not recomputed per request.
- A *sort across a join*: a stored sort key, or an index that matches the exact `ORDER BY` including its `NULL` handling.
- A *geometric read* (map tiles): an index that matches the filter (class), not the address range; later a per-level summary of the largest objects.

Verdict: **do first.**

## 4. Which option for which page

| Page or call | First | Then | Queue? |
|---|---|---|---|
| Systems list | keyset paging; stored sector sort key and an `(quadrant, name, id)` index (PERF.70) | rows before filter counts; capped counts for filters | no |
| Sectors list | stored per-sector system count (summary table) | rows before counts | no |
| Planets, Moons, Phenomena lists | stored counts (PERF.69); keyset paging | rows before counts | no |
| Galaxy Map tiles | index on `phenomenon_scatter(kind, subtype, mass_solar)` (PERF.68); a test that no piece reads more rows than the tile needs | per-piece budgets and `incomplete` tiles | no |
| Search | per-panel statement limits and separate fetches | prefix match rather than `%x%` for short words (the 2-letter case, 8.4 s) | no |
| Counts and facets | stored (done) | capped counts where nothing is stored | no |
| `GET /api/*` list routes | same indexes and keyset (`after=` next to `offset=`) | `Retry-After` on `QUERY_TIMEOUT` (done) | no |
| Generate estimate, deep database check (DB.21), exports | already admin and long-running | progress state | **yes**, with a reserved warm worker |
| Neighbourhood and sector fills | already queued (PERF.24) | | yes |

## 5. What the visitor sees, and what changes under a fill

- **A page that is cheap again** shows nothing special; that is the goal for every row in section 4 except the last two.
- **Rows before menus.** A table shows its first 50 rows immediately and its filter menus a moment later, with "..." in place of the counts until they arrive.
- **A panel or tile that timed out** says so where it is (the search panel: "This took too long; try again or narrow the search"; the map: the tile shows without its points and fills in on the next retry). The whole page does not turn into the 504 page.
- **The 504 page stays** for a page whose main query timed out; it already reloads itself.
- **Table JSON:** `/table/<name>` should return the 504 with `Retry-After` (today it returns 502 without it) and `datatable.js` should retry with the Galaxy Map's 2, 4, 8 ... 30 s backoff, showing the rows already loaded and a "database busy, trying again" line.
- **While a fill runs:** keep the 10 s limit. Under the stand-in fill every measured cheap form stayed under 0.5 s, so no longer limit is needed; a longer limit would only let slow queries make the fill's own writes wait longer. Stored counts keep serving the last count during a fill (PERF.64), and tiles keep their last generation (PERF.34). The summary tables (per-sector counts) refresh in the background, not sooner than every 30 s, like the count cache.

## 6. What the recommendation costs

- Keyset paging changes `datatable.js` and the list functions (an `after` parameter next to `offset`); the scrollbar-by-row-number becomes jump-by-value.
- Stored columns and summary tables need writers to keep them right: sector renames (sort key), system adds and deletes (counts). The refresh-by-`GROUP BY` form needs no writer changes and cost 0.4 to 0.8 s in the test; start with that.
- Migrations on a big database take real time (an index on 17.8 million rows: 29 s; the stored sector sort column: 202 s), so they go through the migration runner's batch and progress reporting.
- The per-piece tile budget adds an `incomplete` state to the tile cache, which is the one part with new moving pieces.

## 7. Items to file (for the TODO thread)

Numbers are placeholders; the TODO thread assigns IDs.

1. **Keyset paging for the data tables.** `after=<sort key, id>` next to `offset=`; `datatable.js` pages forward by key and jumps by value; use the `OR` form, never a row constructor, and a test on MariaDB and MySQL for the plan. Measured 8.5 s to 2 ms at row 1,500,000.
2. **Stored per-sector system count.** Summary table refreshed in the background (0.8 s per 2,000,000 systems) and served by the Sectors list; covers the sorts by systems and density. 17 s to 1 ms.
3. **Systems sort by Sector and Octant (refines PERF.70).** Index `(quadrant, name, id)` with an `ORDER BY` written to match; a stored sector-name sort column with `(sector_sort, name, id)`, filled in batches by a migration, kept in step on sector rename. 4.8 s to 0 s.
4. **Scattered-points index for map tiles (feeds PERF.68).** `phenomenon_scatter (kind, subtype, mass_solar)`; a test that no tile piece examines more rows than its limit allows; a second look at the neutron-star and non-coarse path on a current galaxy. 31.7 s to 0.28 s on the scratch galaxy.
5. **Per-piece tile budgets and `incomplete` tiles.** Each of the eight tile queries under its own statement limit; the tile is served with what finished and rebuilt after a short wait.
6. **Capped counts.** Where no stored count exists for a filter, count to 10,001 and show "10,000 or more" (30 ms against 6.0 s); extends PERF.64's fallback.
7. **Rows first, menus second.** `datatable.js` fetches the first rows without `facets=1`, then the menus in a second request.
8. **Table JSON busy answer and retry.** `/table/<name>` returns 504 with `Retry-After` for a statement timeout; `datatable.js` retries with the map's backoff and says the database is busy.
9. **Search panels on their own limits.** Each result panel fetched separately under its own statement limit; short-word searches use a prefix match (2-letter search was 8.4 s).
10. **Reserved warm worker for long admin operations**, and the poll pattern for them (Generate estimates, DB.21, exports). This is the `planetgen-interactive` queue from the earlier design (3.4), built only for work that is meant to take long. No public polling route.
11. **A reusable big-galaxy query budget test.** Move `make_synthetic.py` into the test tools (2,000,000 systems in about 2 minutes), run `EXPLAIN` on every page and API list query, and fail when a query examines more rows than its page needs; this is the test PERF.68 and PERF.70 each ask for, built once.

## 8. Questions for Boss

1. **Is jump-by-value acceptable in place of jump-by-row-number** on the big tables (scroll bar position becomes an estimate; a "go to letter" box)? Default: yes.
2. **Show "10,000 or more" for filtered counts nobody has stored?** Default: yes; the exact figure appears once the background count finishes.
3. **Should a tile or search panel that timed out show partly, with a retry, rather than the full busy page?** Default: yes.

## Limits of this study

- One 4-core box, one 1 GB buffer pool, MariaDB 10.11; MySQL 8.4 may plan differently (the row-constructor case is the one to check first).
- The 2,000,000 systems are clones of one real system with synthetic names and sectors, so selectivities (star types, sector sizes) are uniform; real data is lumpier. Counts of rows examined scale with table size, times are warm-cache medians, and a cold cache would be slower.
- The scatter table is from an earlier galaxy (no subtype for most black holes, no masses), so the tile numbers show the shape of the problem, not the exact current figure.
- The fill stand-in (inserts and burners) is lighter than a real fill (the real one also runs the generators and 4 or more workers); the 1.2 to 2 times slowdown is a floor.
- The default galaxy is about 64 times the quarter-scale one in volume; row counts and times scale accordingly [R from the generation-performance study].
- The queue costs in 1 and 2.3 are from the earlier PERF.20 measurements, not re-run.

## Evidence notes

- [S] `src/planetgen/util/settings.py` (statement limit), `web/errors.py`, `web/app.py` (504 handler), `web/tables.py`, `static/datatable.js`, `db/query.py` (`list_systems`, `list_sectors`, `count_*`, `galaxy_tiles`, `galaxy_scattered_points_in_box`), `db/countcache.py`, `web/lib/singleflight.py`, `queue/api_jobs.py`, `api/routes.py` (`/jobs/<id>`), `examples/apache/planetgen.conf.example`.
- [C] everything in section 2; scripts and outputs in the shared folder.
- [R] "default galaxy about 64 times quarter scale" from the generation-performance study; the linear scaling to 10,000,000 systems.
