**Headline:** Two object timing tests follow the Seconds per object table of the Stats page again.

### Changed
- **Test fix, no code change (TEST.546).** PERF.197 replaced `run_totals` and `type_means` with `kind_summaries`, which lists the kinds with the one that took the most time first. The unit test expected the system first, but a unit splits its time by each object's own time, so the stars (built inside the systems and alone) take more than the systems' own time, and the table puts the stars first. The test now expects that order and that the totals fall down the list. The galaxy run test called the removed `run_totals`; it now adds up the saved rows of the run. The Stats table and the saved rows are unchanged.
generation-output: unchanged
