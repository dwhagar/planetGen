**Headline:** A supervisor process, one per host, now starts when work arrives and shows on Redis that it is alive; it does nothing else yet.

### Changed
- **OPS.207, output unchanged.** `planetgen.queue.supervisor` and `python -m planetgen.cli.supervisor` are the first step of the job engine's stage 2: one process per host by a kernel file lock (`supervisor.lock` in the state folder), a sign of life in Redis (`pid`, `since`, last turn, turns), a turn every 2 seconds, and an exit after 5 minutes with nothing to watch. A submit (`workers.ensure_workers`), a status poll (`GET /api/admin/work/workers`) and the cache warm-up's start each start it when none is alive. Its loop takes the clock, the store and the process table as arguments, so the engine simulator (OPS.212) can drive it. Watching, worker starts and the reconciler move into it in OPS.208 to OPS.211.
generation-output: unchanged
generation-affects: none
