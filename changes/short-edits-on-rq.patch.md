### Changed

- Creating a system, regenerating a system, a planet, a moon, a belt or a phenomenon, changing a planet's or moon's class, and changing a system's star now run on the Redis queue (PERF.24, step 4c). The request waits up to eight seconds for the job, so a quick edit still answers in the same response; a slower one answers `202` with a job id for `GET /api/jobs/<id>`. Deletes and plain renames still happen in the request. Refusals (404, 409, 400) come back unchanged, and `503` means no Redis server answered.
- A queued job that the work refuses (an `ApiError`) reports `error_status` with its HTTP status in `GET /api/jobs/<id>`.
