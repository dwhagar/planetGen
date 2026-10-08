### Changed

- Wiki uploads (`POST /api/systems/<id>/wiki`, `POST /api/sectors/<id>/wiki`) now run on the Redis queue and wait up to eight seconds, so a quick upload still answers `201` with the page; a slower one answers `202` with a job id (PERF.24, step 4e). The worker reads the wiki's login details from the environment and `config.json` itself, so they never pass through Redis.
