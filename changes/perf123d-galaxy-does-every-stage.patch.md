**Headline:** `planetgen galaxy` does every stage with no way to leave one out, and the Generate page builds its jobs from the single-stage commands.

### Changed
- **`galaxy` is the whole list (PERF.123, fourth slice).** The options `--then-scatter`, `--population`, `--no-settle` and `--backfill-from` are gone, with no replacement flag: a `galaxy` run generates the sectors, links them, scatters the stars and phenomena, backfills, runs the population pass and saves the paths. Whoever wants fewer stages runs the commands in turn (`planetgen sectors`, `link`, `scatter`, `backfill`, `population`, `settle`). `planetgen link`, `backfill` and `settle` take `--since` as before, and without it work on the sectors saved since the start of the newest `sectors` or `galaxy` run. `planetgen sector` keeps `--population` until PERF.124.
- **The Generate page's jobs are lists of commands.** New galaxy plans, then runs `sectors` around a random start, `link`, the three scatters, `backfill` and `settle` as steps of one job (the scatters are left out when "skip the bright-star scatter" is ticked); Generate sectors runs `sectors`, `link`, `backfill` and `settle`. The estimate is taken from the `sectors` step. The Queue page's retry of a failed sector runs `sectors` for it and then `link`, `backfill` and `settle`.

TODO items resolved: none (PERF.123 slice 4 of 4). Added: none.
