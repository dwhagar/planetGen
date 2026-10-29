**Headline:** The reset, the math check and the data check now run as recorded runs on the work queue, like every other command.

### Changed

- `planetgen reset` (and the Generate page's reset step) wipes the database as one `reset-galaxy` unit of a recorded run, titled "Reset the galaxy". The typed confirmation, the dry run and the progress bar behave as before. A cold worker takes a while to start, so the first reset after a quiet spell may wait a little longer before the wipe begins.
- `planetgen check-math` runs the math check as a `check-math` unit of a check run, and exits 1 on a failure as before.
- `planetgen datasets status` runs the data check as a `check-datasets` unit of a check run. `datasets fetch` is unchanged.
- The command inventory and the work queue design page list the three new commands.

TODO items resolved: PERF.119. Added: none.
