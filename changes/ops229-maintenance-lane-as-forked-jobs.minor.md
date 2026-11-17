**Headline:** The cache warm-up and the database check now run as jobs the supervisor hands out from a warm template, the way a generation run's units do, and the maintenance lane's round-robin worker is gone.

### Changed
- **The maintenance and check lanes move to the launcher (OPS.229, D4 of the job engine).** A cache warm-up or a database check that has a ledger writes its rows with the class `maintenance` (`planrows.job_class`), and `wrapper.forked` gives the units of every such run to the supervisor's launcher instead of RQ's queues. The launcher claims them on the maintenance lane, holds them to the share the setting `queue.maintenance_workers` gives upkeep while a generation run is live (none otherwise: the memory gate and the stage's own cap are the limit), and forks each from a template that holds the site's app built (`jobhost.WARM_SETS` `app`, `warmup.warm_app`), which a worker of the lane built for itself before (PERF.245).
- **A worker serves one queue.** RQ's `RoundRobinWorker` and the lane's queue-in-turn rule are deleted (`cli/worker.py`, `workers.lane_of`), so RQ's reliable dequeue holds for every queue that is left (closes OPS.189). The lane's two queues serve only a run with no ledger, which OPS.271 moves.
- Deleted: `cli/worker.APP_QUEUES` and the app a worker built once; `worker_class(round_robin=)`.

generation-affects: none

generation-output: unchanged
