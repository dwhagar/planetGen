**Headline:** A log line gives a duration under a second in whole milliseconds, so `test_job_logs` is green again.

### Fixed

- `logrecord.duration_text` read a quick run as `5.56 ms` (the time ladder's three significant figures), which turned the two `test_job_logs` tests red on main. From 0.5 ms to just under a second a log line now says whole milliseconds (`6 ms`, `340 ms`), 0.9995 s and up reads `1s`, and the pages and the ladder itself are unchanged (TEST.814).

generation-output: unchanged
