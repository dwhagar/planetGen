### Added
- **Bright stars are drawn several layers at a time (PERF.7).**
  `generate.py plan` hands each layer of the galaxy to a worker process,
  like sector generation; three workers placed the same stars about
  three times faster (313 s down to 105 s at 20,000 L_sun on a 4-core
  machine). `--workers` sets the count. Each layer has its own random
  stream, so the same seed gives the same stars on any number of
  workers.
- **The Generate page shows the time left** ("about 4 m 10 s left") for
  a running job, the same estimate the terminal shows.

### Changed
- **Steadier time-remaining estimates (PERF.7).** Every generation
  progress bar's remaining time now comes from a decaying average of
  sectors (or layers) finished per second, weighted toward the last
  minute, instead of rich's own short-window estimate, so it holds
  steady while several workers report at once and stays on screen until
  the bar is done. A run stopped early by `--limit` now ends its bar
  finished rather than part way.
