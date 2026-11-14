**Headline:** Only the supervisor starts workers now; a submit, a status poll or a web thread leaves it a request.

### Changed
- **OPS.209, output unchanged.** `workers.ensure_workers` no longer counts or starts anything: it leaves a request for the queue (`supervisor.nudge`: the latest ask for each queue, with its worker count) and starts the supervisor if none is alive, so it is cheap to call from any web thread or enqueuer. The supervisor takes the requests at each turn (`supervisor.WorkerStarts`) and acts on each with `workers.start_workers`, which is the old function with the same rules (caps, lanes, memory limit, shedding the bulk lane for the maintenance one, lending). With one process deciding, the host-wide lock between callers (PERF.181) is gone. The count of workers started and not yet registered stays, because a worker clears it as it registers. A supervisor with a request waiting does not end for being idle.
generation-output: unchanged
generation-affects: none
