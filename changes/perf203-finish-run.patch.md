**Headline:** "Finish the last run" now carries on a stopped Generate job as the one run it is, including the stages that never started.

### Changed

- PERF.203 step 5: a Generate job is one run (step 4), so "Finish the last run" plans what is left from that run's own stage list and no longer adds the job's later steps as separate runs. The task page no longer lists "steps of that job that never started"; the stages that run again include the ones the stopped run never reached. A scatter whose stage had not opened when the run stopped still works out its limits when its stage opens (`scatter-settings` is left deferred while the stopped run is planned), so a finished job stores what an uninterrupted one does. The queue page's Retry is unchanged. Nothing stored changes, no epoch.
