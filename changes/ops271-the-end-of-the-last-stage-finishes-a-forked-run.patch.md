**Headline:** A run whose last stage has units ends as done when those units end, as it does for a run in workers.
### Fixed
- A run whose units are handed out by the supervisor's launcher (OPS.228) and whose last stage held units tried to queue an opening for a stage that does not exist, failed six times, and ended the run as failed. The end of the last stage now finishes the run in place (OPS.271).
- The launcher read a job's row from the snapshot of an earlier look, so a unit that had finished in its own process could still read as running when its process ended, and was handed out a second time (it ended as skipped). It looks again before it settles a job.
