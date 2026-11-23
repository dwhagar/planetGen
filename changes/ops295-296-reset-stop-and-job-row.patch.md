**Headline:** A reset that is stopped says what stays wiped, and a job whose row a reset wiped ends cleanly.

### Fixed

- OPS.295: the handler of a reset stopped by a signal called a log level that does not exist (`log.info`), so it died with an `AttributeError` before it could say that the tables wiped so far stay wiped; it uses `log.normal` and the reset exits 143 as meant.
- OPS.296: the job runner read the job's state from its `engine_jobs` row without checking the row was there; a reset wipes the control tables, its own job's row too, and the runner failed with `'NoneType' object is not subscriptable`. A row that is gone now ends the job as a Cancel does.

generation-output: unchanged
