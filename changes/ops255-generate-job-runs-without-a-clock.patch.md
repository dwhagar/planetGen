**Headline:** a Generate job that runs on the queue is no longer ended by RQ's time limit (ten times what earlier, smaller galaxies recorded), which would stop the starter of a run that is still working.

### Fixed

- OPS.255: Such a job takes the 30-day ceiling, since the run's watch judges its units by their progress. When RQ or the runner ends any job, the job's log and records now say so as the failure of the stage it was at.
