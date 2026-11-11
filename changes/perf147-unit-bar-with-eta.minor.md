**Headline:** Every running unit now has its own progress bar and time left, with an estimate from the first second, and the unit list in the work API carries it.

### Changed

- `GET /api/admin/work/<id>/units` gives each running unit a `progress` record (`queue/taskprogress.py`): `show` (true once the unit is expected to take, or has run, longer than 15 s), `percent`, `eta_seconds`, and `source` (`stored`, `blended` or `live`). It is worked out from the unit's heartbeat counter (OPS.136 to OPS.138), so the page, the API and the log read the same numbers.
- The time left moves from the rate recorded for the command to the running unit's own as it proceeds: the running rate's weight is the share of the unit done over a quarter, and the running rate is held back for its first 20 s (`progress_rate.HOLD_SECONDS`).
- New: the server records how fast each command's units work (`queue/unitrates.py`, items a second over the seconds worked, not counting waits on the gate or a lock). A run keeps the rate of its finished units, a run that ends done folds it into the server's, and a unit's bar uses the run's rate first, then the server's.
- The shell scatter's second bar now appears when tasks take longer than 15 s (it was 30 s) and goes away under 10 s (it was 20 s), the threshold every bar shares.
- A bar macro for the unit rows of the Queue page (`partials/task_bar.html`) is added for the unit rows that OPS.139 builds.

generation-output: unchanged
