**Headline:** Three guard tests that went red with the template slices are green again.

### Fixed
- The runner of a unit forked from a template runs inside the unit's own memory claim, so it is named among the functions that take none of their own (`claims.CALLBACKS`).
- The way to try the no-kernel-handle path is the template's `--no-pidfd` argument (`JobHost(pidfd=False)`), not an environment variable the settings model does not name.
- The two handlers that turn a job's failure into an exit code now note it (`log.swallowed`).
