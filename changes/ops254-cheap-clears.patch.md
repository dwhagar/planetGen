**Headline:** A New galaxy run on an empty database no longer spends minutes clearing what is not there.

### Fixed
- **OPS.254 (bug), output unchanged.** The reset step of a New galaxy (`planetgen reset`, the Generate page) truncated all 111 tables, 36 seconds on Boss's run, when one held a row. It now leaves a table alone that holds no row and whose id counter was never used, and still truncates a table that holds rows or whose counter has moved, so the next galaxy numbers from 1 as before. "Clear the earlier bright stars" (129 seconds on that run, with no bright star to clear) now reads five things and leaves when `bright_stars` is empty and unused and the scatter's threshold, seed, mass limit and epochs are not recorded (`store.bright_stars_hold_nothing`); a clear that has something to clear runs every step as before, with the `TRUNCATE` of `bright_stars` last so a clear that was stopped half way is done again in full.
generation-output: unchanged
generation-affects: none
