**Headline:** Generation runs and web jobs no longer write the old MySQL job tree; the Redis work queue carries them, and a job keeps how long each step took in its jobs folder.

### Changed
- **The old job tree is no longer written (PERF.98).** `planetgen galaxy` and the other generate commands no longer open a root node in `work_jobs`, the population and territories passes and the page-job answer builder no longer wrap themselves in `job_node`, and the web job runner no longer records a `web-job` root with a node for each step. A web job estimates the time of each step from `step-times.json` in the jobs folder instead of finished step rows in the tree, and a generate command still refuses to start while an update is running (the check moved from `work.open_node` into `generate.main`). The tree functions in `queue/work.py` and the `work_*` tables stay until PERF.100 deletes them. No generated value, epoch or golden changes.

TODO items resolved: PERF.98 (step 4 of the work queue: sectors, backfills and the remaining stages run on Redis only; the old tree has no writer left). Added: none.
