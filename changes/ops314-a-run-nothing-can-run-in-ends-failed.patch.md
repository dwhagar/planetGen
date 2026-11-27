**Headline:** The supervisor ends a live run in which nothing can run any more, as failed with the cause said, instead of leaving it Running for days with no worker alive.

### Fixed
- OPS.314: Boss's Run 57 (Generation) stood Running from 11-25 12:29 AM EST at 21 of 233 jobs with 1 failed, no worker alive, while the Workers and Runs live cards showed 0. A live run whose ledger rows have no job running or ready and a job failed (or only jobs waiting) now ends `failed` after 10 minutes (30 minutes with no failed job): the run's error names the failed job, its stage and the first line of its error, how long nothing could run (HMS) and how many jobs wait on it; its rows are ended and the supervisor's acts say so.

generation-output: unchanged
