**Headline:** Every job of the work queue now takes a claim on the host's memory and load before it starts, so a generation run no longer takes the machine down.

### Added
- **One memory and load gate for all jobs (PERF.254, `queue/memorygate.py`, `queue/claims.py`).** A worker takes a claim (memory in MB, cores, lane) for each job before it prepares it and waits, with its heartbeat going, until the host has room: the free memory from the kernel's MemAvailable (and a cgroup limit when there is one), less a reserve for MySQL, Redis and the web server, less the places kept for the interactive and long lanes, less the part of every other claim not yet in memory (read from the proportional memory of each worker and its horse, `psutil` 7.2.2, a new dependency in both lock files). Cores are claimed the same way and an ordinary job waits while the one-minute load per core is over 1.25. The wait is bounded (300 s) and then the old table cap decides, so a wrong gate cannot stall a run. Every kind of job has a claim (a run's units, the API's and the pages' jobs, the Generate page's, the node prefetch, the cut of a stage) and `tests/test_memory_gate.py` fails for a kind that has none. The test runner takes a claim for each phase.
- **Claims are learned.** The worker's heartbeat samples what the job's processes hold; the largest peak seen for a kind is claimed, with 10 percent on top, in place of the table's figure when it is larger.
- **A job killed for memory is retried once with a larger claim,** then the run stops with the unit (or the opening) named and what it claimed and held.

### Changed
- **The cut of a stage is a job of its own** (`units.open_stage`, claimed alone on the bulk lane) and no longer runs inside the work horse of the last unit of the stage before; a failure or a kill there is filed under its own stage. Cancel, Clear and the update guard know the job.
- `scripts/test.py` runs each phase under a claim on the gate when a Redis server answers.
- Generated output is unchanged.
