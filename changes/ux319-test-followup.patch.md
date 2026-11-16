**Headline:** Two more tests that still expected a decimal time now expect the hours, minutes and seconds form.

### Fixed
- `test_job_view` (a retry hold and an attempt's CPU time) and `test_web_sector_nav` (the stay at each stop) read times as `2m`, `1m 5s` and `30m`.

generation-output: unchanged
