**Headline:** A reset that is stopped says what stays wiped, and a job whose row a reset wiped ends cleanly.

### Fixed

- OPS.295: the handler of a reset stopped by a signal called a log level that does not exist (`log.info`), so it died with an `AttributeError` before it could say that the tables wiped so far stay wiped; it uses `log.normal` and the reset exits 143 as meant.
- OPS.296: the job runner read the job's state from its `engine_jobs` row without checking the row was there, and failed with `'NoneType' object is not subscriptable` (stopping the step it followed, a reset among them) on a run of `test_web_generate` on a slow machine. A row that is not there is now no Cancel: the job runs on. What removed the row there is not found (the control schema is a database of its own that a reset does not touch, and the same test passed alone in 6m 37s with the row always there). The first time a job step finds its row gone the log now says so, with the job, run and step, to find the cause.

generation-output: unchanged
