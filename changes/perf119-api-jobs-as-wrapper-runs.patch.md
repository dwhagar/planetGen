**Headline:** The API's generate-neighborhood, regenerate-sector and settle jobs run through the work queue's run record, so they show on the Queue page and wait their turn.

### Changed
- **The API's sector jobs are runs (PERF.119, second slice).** Generating a neighborhood, regenerating a sector and saving sector paths each start runs through the wrapper (the sectors, the bright-star backfill around the requested sector, and the paths), like the stages of `planetgen galaxy`, instead of generating inside the job's own process. A job that finds another generation run going waits in line for its turn, as every generation run does. The job's id and `GET /api/jobs/<id>` answer as before.

TODO items resolved: none (PERF.119 slice 2 of 3). Added: none.
