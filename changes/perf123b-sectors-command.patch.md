**Headline:** `planetgen sectors` generates sectors in every mode of `galaxy` without the stages after them.

### Changed
- **The sectors of a galaxy run are a command of their own (PERF.123, second slice).** `planetgen sectors` takes the modes and sector options of `planetgen galaxy` (shell, band, wedge, span, course, neighborhood, random start, limit, estimate) and only generates the sectors, saved unlinked, then says the run's start time for `planetgen link`, `backfill` and `settle`. `planetgen galaxy` is unchanged for now and shares the code.
- **A shell, band, span, course or neighborhood run hands its stage the arguments that name its sectors, not the sectors (PERF.123).** The stage cuts them itself when it opens (`planned_items`: the same addresses less those that exist now, up to `--limit`), so a run of a million sectors no longer carries three million numbers in its record. `run_shell`, `run_span` and `run_course` are one `run_addresses`; the neighborhood is cut the same way, in the fill order, and each sector's log line carries its distance from the center; the single-address and random-start modes still hand over their lists.
- **`planetgen sectors` leaves the density to the plan like `galaxy` does.** It no longer falls back to a flat 10 systems a sector when no density option is given.

TODO items resolved: none (PERF.123 slice 2 of 4). Added: none.
