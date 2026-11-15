**Headline:** The Queue page now shows every job of a run with its attempts and why it waits, and a history of runs, read from the job ledger.

### Added

- OPS.224 (job engine stage 3, L10): `GET /api/admin/work/<id>/jobs` and `/jobs/<job id>` serve a run's jobs and one job's attempts from the ledger (`queue/jobview.py`), and `GET /api/admin/work/history` the runs the ledger holds, filtered by state, kind or galaxy. The run page gets a Jobs panel; `/admin/queue/<id>/jobs` lists the jobs under a state menu, each with its attempts, how its newest attempt ended and, when it is not running, why it waits; `/admin/queue/<id>/job/<job>` shows every attempt with its cost, error text and stack file; the queue page ends with a Run history table.
- A job that a worker holds at the memory gate now says why in the ledger (`engine_jobs.waited`, the gate's own words from `memorygate.decide`, written when the reason changes and at most every 15 s), a retry hold shows how long is left, and the end of a stage names what it still needs. Reused: `runview`, `planrows.control_connection`, `datatable`, the run page's panel and table markup.

generation-output: unchanged
generation-affects: none
