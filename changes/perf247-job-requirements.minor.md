**Headline:** A queued job can declare that it requires another job, so no worker is handed one that cannot run yet.

### Changed
- A job can require other jobs (`requires.enqueue`). It is held out of every queue until they are done, so the round robin of a worker that serves several queues cannot pick it early, and the worker that finishes the last one queues it. A job that declares no requirement is queued at once and runs beside the others, as before.
- Held, released and dropped jobs are logged (and in the run's own log for a unit of a run). A required job that fails or is stopped fails what waits for it, and the jobs that wait for those, with the line "required job X failed", instead of leaving them waiting for ever; Stop All drops the held jobs with the queued ones.
- The stages of a run that need another stage (the phenomena index after the phenomena scatter, the territories after the population pass, the link and the sector paths after the sectors) are declared; a list of steps in the wrong order is refused where it is made, a stage logs when it is released, and a failed run logs which stages it left undone and why.

generation-output: unchanged
