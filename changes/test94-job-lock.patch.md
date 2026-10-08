### Fixed
- The job tests treat a job as finished only once its runner has written a final status and released the job lock, so starting the next job right after one no longer meets a busy lock under load (TEST.94).
