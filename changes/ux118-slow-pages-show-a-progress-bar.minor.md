**Headline:** The Phenomena list and the Search page show a progress bar while a slow answer is being built, then swap the finished page in without a reload.

### Changed
- `/phenomena` and `/search` are marked slow with `@slow_route` (`web/lib/slowroute.py`). When the answer is being built by a queued job (PERF.105) the visitor gets the page's own frame with a panel in place of the content: a progress bar, one line of status ("About 30% (an estimate: 300 of 1,000 rows), about 5 s left."), and Stop waiting. The answer is `202`, `Cache-Control: no-store`, `X-Robots-Tag: noindex`.
- `static/slowload.js` polls `GET /api/page-jobs/<id>` every second and draws the estimate in the site's own units. A job still queued leaves the bar without a value and says it is waiting; one past its expected time says it is still working. When the job has succeeded it fetches the same address again and swaps the finished content in without a page load (title updated, the page's own scripts loaded); a page that waits for a second answer (a search asks for several) carries on with that job.
- A failed job shows its error and Try again, which asks for the page with `X-Slow-Retry: 1` so a new job replaces it; a job an administrator cancelled says so. Stop waiting only stops this page polling (the job is shared), and Check again picks it up.
- The bar is a `<progress>` named by the panel heading, the status is a polite live region and the error an alert; the panel wraps at phone width; a bar with no value is hidden for a reader who asked for less motion and the words stand alone. Without JavaScript the page refreshes itself every five seconds.
- A slow answer built in the request because the queue could not take it now says so at the top of the page.
- The cost rule gained a middle step (`planetgen:slowpath:...`, `SLOW_PAGE_PATH_SECONDS`, ten minutes): a request not timed itself takes the slowest build of any request on its route lately, so a search for a new name does not queue a job just because the tables are big.
- `docs/design/slow-page-answers.md`, `docs/html-interface.md` and the architecture document describe it.

### Added
- `tests/js/slowload.test.mjs` (10 tests), `tests/test_slow_route.py` (12), `tests/test_web_browser_slow_pages.py` (34 runs: loading, running, queued, failed with Try again, cancelled, stop and check again, a second job, no scripts, reduced motion, light and dark, desktop and phone width), and one more cost-rule test in `tests/test_slowpages.py`.
- Follow-up for TODO: the table's own scrolling and sorting still fetch rows from a data route that is not marked slow (it builds in the request, as before); marking it needs `datatable.js` to wait on a job.

TODO items resolved: UX.118 / Added: the data-route follow-up above (a new item for TODO to file)
