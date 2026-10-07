### Changed

- `POST /api/sectors/<id>/regenerate` now queues the delete-and-regenerate on Redis and answers `202` with a job id, instead of holding a web worker (PERF.24, step 4b). Unknown (404) and off-grid or unplanned (409) sectors still answer before anything is deleted. The admin page's Regenerate button says the job is queued.
