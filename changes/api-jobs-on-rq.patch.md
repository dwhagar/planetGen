### Changed

- `POST /api/sectors/<id>/generate-neighborhood` now queues the run on Redis and answers `202` with a job id, instead of holding a web worker for the whole run (PERF.24, step 4a). The math check, the unknown-sector and missing-skeleton errors (404, 409) and the disk-space refusal (507) still answer before anything is queued; `estimate_only` still answers at once.
- Added `GET /api/jobs/<id>`, the state, result and error of a queued API job.
