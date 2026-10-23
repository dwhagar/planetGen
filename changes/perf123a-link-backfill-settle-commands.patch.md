**Headline:** `planetgen link`, `backfill` and `settle` run the last stages of a galaxy run on their own.

### Changed
- **Three stages of `galaxy` are commands (PERF.123, first slice).** `planetgen link --since TIME` links the sectors saved since `TIME` to their neighbours, `planetgen backfill --since TIME` draws the bright stars out from them, and `planetgen settle --since TIME` saves their paths and those of the sectors around them. Each is one run on the work queue and does nothing when no sector was saved since `TIME`. A `galaxy` run now says its start time on the database's clock, the value these commands take. `galaxy` itself is unchanged for now.

TODO items resolved: none (PERF.123 slice 1 of 4). Added: none.
