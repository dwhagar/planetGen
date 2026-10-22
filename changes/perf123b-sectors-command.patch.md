**Headline:** `planetgen sectors` generates sectors in every mode of `galaxy` without the stages after them.

### Changed
- **The sectors of a galaxy run are a command of their own (PERF.123, second slice).** `planetgen sectors` takes the modes and sector options of `planetgen galaxy` (shell, band, wedge, span, course, neighborhood, random start, limit, estimate) and only generates the sectors, saved unlinked, then says the run's start time for `planetgen link`, `backfill` and `settle`. `planetgen galaxy` is unchanged for now and shares the code.

TODO items resolved: none (PERF.123 slice 2 of 4). Added: none.
