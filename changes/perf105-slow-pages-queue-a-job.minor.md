**Headline:** A page request that is slow to build now queues one job that stores the answer in the page cache, with a status route and a progress estimate, instead of holding a web thread.

### Changed
- `apiclient` (`web/lib/apiclient.py`) can take a slow cache miss off the request. A page that opts in with `slowpages.allow_pending()` gets `slowpages.Pending`, carrying the job's status, for a request the cost rule calls slow; any other page, and the job itself, builds as before. The cost rule is in one place (`web/lib/slowpages.py`): a request is slow when its last timed build took more than `SLOW_PAGE_SECONDS` (3 s), or, never timed, when the tables its route reads in full hold more than `SLOW_PAGE_ROW_BUDGET` (2,000,000) rows by the database's estimate (both in `tuning.py`, with the rows per second and the times the jobs and timings are kept).
- `planetgen/queue/pagejobs.py` queues the job (`build_answer`, on the interactive queue beside the cache warm-up's): it makes the request through the site's own app, so the answer is stored under the tags it read and an edit invalidates it like any cached answer. Two visitors asking for the same request, even in the same moment and in different web processes, share one job (`SET NX`); a finished job whose answer was invalidated is replaced; a failed or cancelled one is shown for 30 s, then replaced, or at once on retry. Its record, progress and the last build time live in Redis, all with a TTL.
- The progress is an estimate and says so (`estimated: true`): the elapsed share of the expected time, never past 95 % before the answer is stored, with the time left. The answer is one query, so there is no step to count; PERF.99 splits the lists it can.
- The job is a node of kind "Slow page answer" on the admin queue page with a result line. An admin cancels it there or with `POST /api/page-jobs/<id>/cancel`: a job not yet started is removed, one not yet at its query stops, a running query finishes and its answer stays cached. A visitor's own Cancel only stops that page waiting. `api_jobs.cancel` removes a queued job.
- Without Redis, or while an update runs, the answer is built in the request as it is today, the log says why, and `slowpages.take_notice()` gives the page a sentence to show.
- `docs/design/slow-page-answers.md` is the design; `docs/api.md` documents the routes; the architecture document lists the module.

### Added
- `GET /api/table-estimates` and `GET /api/page-jobs/<id>` (public), `POST /api/page-jobs/<id>/cancel` (admin); the API route snapshot is updated (additions only, no version change).
- `tests/test_slowpages.py` (16 tests): the cost rule, one job per request across processes and threads, the answer stored and found by another process without an API call, an edit invalidating it, a failed job and its retry, the status route, the progress estimate, the cancels, the admin node, and the fallback with no queue or an update running.

TODO items resolved: PERF.105 / Added: none
