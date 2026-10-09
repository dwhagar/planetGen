### Changed
- Values that were JSON blobs now live in columns: a generation run's arguments are rows in `generation_run_arguments`, a work job's argv is rows in `work_job_args`, and the unused task `result` column is gone (schema v62, control schema v10). Adds indexes on the star-system quadrant and the star, planet and moon radius and star type columns.
