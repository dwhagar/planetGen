**Headline:** The old queue executor is removed; nothing in a run used it any more.

### Changed
- **The old work-queue executor is gone (OPS.66, Boss).** `WorkQueue`, `RQExecutor`, `RQFuture`, `Channel`, `_rq_task` and the `--debug FILE` plumbing of that path had no caller in a run since runs became units queued by the wrapper; nothing is kept for compatibility. The job tree (`job_node`, the admin Queue page and API), `worker_count` and `task_seed` stay. The tests of the removed executor are removed with it, and the ones that only needed a queue node to look at use a plain job node.
- **A worker never writes the progress file.** The worker itself now drops `PLANETGEN_PROGRESS_FILE`, which the removed executor did for its tasks, so a worker started from a Generate page job does not overwrite the job's progress.
- **The debug docs say what the detached helpers do (OPS.67).** The background `planetgen warm` and an auto-started worker's own start and exit write the debug log file only, since they have no stdout.
