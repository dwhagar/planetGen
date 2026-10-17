# Slow page requests as jobs (PERF.105)

Informs: PERF.105, UX.118 (the loading state and progress bar), PERF.104 (the page cache), PERF.86 (Redis rules), PERF.99
Status: built (Core 2, 2026-10-17): the mechanism. No page uses it yet; UX.118 puts the phenomena lists and the search on it. Nothing here changes a table; no migration.

## The problem

A page whose answer reads a big table (the whole list of one phenomenon kind, a wide search) builds it inside the web request. The request holds one of the five web threads for as long as the query runs, the visitor sees a blank tab, and two visitors asking for the same thing run the same query twice. The page cache (PERF.104) already holds finished answers under the tags of the data they read; what is missing is a way to *produce* a slow one without holding a request.

## Design

**A page opts in; a miss that the cost rule calls slow becomes a job.** A route that expects slow answers runs inside `slowpages.allow_pending()`. Inside it, `apiclient` reads the page cache as always; on a miss it asks the cost rule (below) once, and if the request is slow it does not build the answer: it asks `pagejobs.ensure` for the job and raises `slowpages.Pending`, carrying the job's status. The page catches it, shows its loading state (UX.118) and polls `GET /api/page-jobs/<id>`. When the state is `succeeded` the page asks again and the answer is in the cache. Outside `allow_pending`, and in the job itself, a miss is built in the request exactly as before, so nothing changes for a page that does not opt in.

**The job builds through the page's own path.** `pagejobs.build_answer` is an RQ job on the interactive queue (the warm-up's queue, `api_jobs.submit`). It makes the request through the site's own app (`apiclient.fetch_target`), the way the warm-up does, so the answer is stored by the same `begin`/`put` as any cached answer: under the tags it read, in Redis (`planetgen:page:<db>:<digest>`), where every web process finds it. An edit that bumps one of those tags (an API write, the population pass, the orbital update, the content stamp) makes the stored answer stale and the next request queues a new job; nothing here keeps its own copy or its own expiry.

**One request, one job.** `planetgen:slowjob:<db>:<digest of the request>` holds the id of the job building that request, claimed with `SET NX`, so visitors in other web processes, asking in the same moment, attach to one job (`ensure` returns the same id to all). A finished job is replaced when the answer it stored has since been invalidated. A failed or cancelled job is shown to the visitors who ask for the same request for `SLOW_PAGE_FAILED_SECONDS` (30 s), so a failing request is not re-queued by every page load, then replaced; a page's Retry replaces it at once (`ensure(..., retry=True)`).

**The cost rule, in one place** (`slowpages.cost`). A request is slow when:

1. its last build took longer than `tuning.SLOW_PAGE_SECONDS` (3 s). Every cached build records its time in Redis (`planetgen:slowtime:...`, kept a week). A measured time beats any estimate, so a request that turned out quick is quick from then on, and one that turned out slow is slow whatever the tables hold; or
2. it has never been timed, and the tables its route reads in full (`slowpages.PATH_TABLES`: the phenomena list reads the scatter and every phenomenon table, the search the systems, bright stars and phenomena) hold more than `tuning.SLOW_PAGE_ROW_BUDGET` (2,000,000) rows by the database's own estimate (`GET /api/table-estimates`, from `information_schema`, instant, within some tens of percent; kept a minute per process).

The expected time of a request nothing has timed is its rows over `SLOW_PAGE_ROWS_PER_SECOND` (250,000); that is the first guess at the progress bar's ETA and is replaced by the measured time after the first build. A route that reads one page of an indexed table does not list its table, or does not opt in.

**Progress is an estimate, and says so.** The answer is one query, so there is no step to count. The status gives `done` of `total` rows with `estimated: true`: the elapsed share of the expected time, never past 95 %, so it cannot read as finished before the answer is stored; `eta_seconds` is the time left, or `null` once the expected time has passed ("still working"). At `succeeded`, `done == total` and the ETA is 0. UX.118 draws this as a progress bar whose text says it is an estimate. A real count would need the query split into steps (PERF.99 does that for the lists it can page through); then the same fields carry the real numbers.

**On the queue page and cancellable.** The job is a node of kind `slow-page` ("Slow page answer") on `/admin/queue` beside the warm-up's, with one result line (the request, the bytes, the seconds). An admin cancels it there, or with `POST /api/page-jobs/<id>/cancel` (admin only): a job that has not started is removed from the queue, one that has not begun its query stops, and a query already running finishes (MySQL cannot be asked to stop it from here) and its answer stays in the cache, because it is correct, but the job reports `cancelled`. A visitor's own Cancel only stops that page from waiting; every visitor shares the job, and one anonymous visitor must not cancel it for the rest (the API takes no anonymous write but sign-in).

**Redis down, or an update running: build it in the request, and say so.** `ensure` raises `NoQueue`; `slowpages.pending_for` returns nothing, `apiclient` builds the answer in the request as it did before, the log says why, and `slowpages.take_notice()` gives the page one sentence to show ("This page is slow to build and the work queue is not available, so it was built while you waited."). A page cache that cannot reach Redis is already a plain miss, so the page is slower, never wrong.

**With the early refresh (PERF.110).** A costly cached answer is fetched again shortly before it expires, by one reader at a time (`pagecache.ResponseCache._refresh_early`); that reader's cache read is a miss. On a page that opted in and a request that is slow, that miss queues the job like any other, so the early refresh of a slow answer is built on the queue, not in a request; the visitor of that moment sees the loading state for the length of the build. Every other reader is still served the old copy until the new one is stored.

## Redis keys (all with a TTL, PERF.86)

| Key | Holds | Kept |
|---|---|---|
| `planetgen:slowjob:<db>:<digest>` | the id of the job building that request | 1 h (30 s once failed or cancelled) |
| `planetgen:slowmeta:<id>` | the request, estimated rows and seconds, created, started, finished, cancel flag, queue job id | 1 h |
| `planetgen:slowtime:<db>:<digest>` | how long the last build took, seconds | 7 days |

## Not here

The loading state, the progress bar, the retry button and the swap-in of the finished list belong to UX.118; moving the phenomena lists and the search onto `allow_pending` belongs with them. Splitting a list into counted steps is PERF.99.
