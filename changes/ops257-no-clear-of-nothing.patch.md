**Headline:** a galaxy run that starts with a database wipe no longer clears scatters that hold nothing.

### Changed

- OPS.257: A galaxy run that starts with a database wipe no longer clears the scatters: a clearing stage whose target holds nothing (the phenomena and the nebula field) is left out of the plan and listed as skipped, by every route (the web's New galaxy, the CLI scatter and galaxy, run-job). The bright-star clear stays, since it stores the scatter's limits, and finds nothing to clear in a few reads.
