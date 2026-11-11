**Headline:** The Queue run page lists each running unit with its state, last beat, counter, phase and rate, the ones that need a look first.

### Changed
- **Running units on the run page (OPS.139).** `/admin/queue/<run>` of a live run gets a "Running units" table drawn from `GET /api/admin/work/<id>/units` (OPS.138, through `apiclient.admin_work_units`): stage, unit, a state badge (working, slow, waiting, quiet, looping, dead; the reason is its tooltip), what it is doing and what it waits on, the counter (`done of estimate`, or "no counter"), the rate, the time left at that rate, the age of its last beat, how long it has run, its memory, and "to be stopped" on a unit two checks in a row have called quiet or dead. A summary line counts the states; the worst states are listed first and at most 200 rows are drawn. A run that is over, or one with no unit judged yet, has no table. Reused, not rebuilt: the verdicts and their view (`unitwatch.view`), the route, `_duration`, `counted`, the page's panel/table markup and `.badge`. No migration.

generation-output: unchanged
generation-affects: none
