**Headline:** a held unit's error now says where it failed, so a one-off failure under load can be traced from the Queue page and the log.

### Changed
- OPS.275: the error a held unit keeps (the held hash, the run's `unit.held` line, the job log and the reason the run is held with) ends with ` (at file:line in function, ...)` for the last three frames of the failure's traceback, innermost first, files named from the `planetgen` package (`held.where`). The whole error text, exception message included, is passed through the secrets mask before it is stored or said.

generation-output: unchanged
generation-affects: none
