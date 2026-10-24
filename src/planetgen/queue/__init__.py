# planetgen/queue/__init__.py

"""
The work queue: runs of generation and the API's long jobs, on Redis with RQ.

`redisclient.py` makes the one Redis client; `redisqueue.py` holds the
queues and the RQ worker class, `workers.py` starts the workers of the shared
queues, counts them and stops a running job, and `planetgen.cli.worker` is
the worker process. A run lives only in Redis: `runrecord.py` keeps its
record and counters, `runeta.py` works out the time left from them,
`progress_rate.py` the rate behind every bar's estimate, and
`progress_file.py` the progress a web-started run writes for the page.
`units.py` runs a stage of a run on that spine (PERF.97: the units queued as
JSON, the done-set barrier, the starter following the run; so far the
massive-star scatter). `work.py` is the generation work queue (units of work run on a pool of
workers), `api_jobs.py` queues the API's long operations, `load.py` reads
the server's load for the queue page and `redishealth.py` what the admin
Stats page shows about Redis. The design is
`docs/design/operations/redis-and-the-work-queue.md` (PERF.86).
"""
