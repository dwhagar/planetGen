**Headline:** The Generate page's jobs, the API's jobs and the map's tile prefetch now give a sign of life while they run, and the Queue page lists them, so a dead one shows as dead instead of as nothing (OPS.142, heartbeat part 9).

### Added
- **`queue/jobbeat.py`.** `jobbeat.beating(kind, job_id, label)` runs the same beat thread as a unit (`unitbeat.UnitBeat`, Redis server clock) for a job that is not a unit of a run: an API job (`api_jobs.execute`), a web job (`web/job_runner.run`, with the step it is on as the phase) and a tile prefetch chunk (`tileprefetch.build_nodes`). The records sit in one hash, `planetgen:jobbeats`, removed when the job ends; a job that is killed leaves its record, which reads `dead` for an hour and is then dropped. These jobs have no unit of work to count, so the judge says only `working` or `dead`; nothing is ended on its state (RQ's time limit still is). Where Redis does not answer, the job runs without a beat.
- **`GET /api/admin/work/jobs`** (admin; additive, API version unchanged) and a **Jobs running now** panel on the Queue page: kind, label, state, what it is doing, how long it has run, and when it was last heard from.

### Changed
- `UnitBeat` can write to a hash and field of its own (`hash_key`, `extra`), which is how the job beat reuses it.
- No star, planet or stored number changes.

generation-output: unchanged
generation-affects: none
