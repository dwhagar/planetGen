**Headline:** The cache warm-up of sector, system and phenomenon pages is cut into short units, and emptying a scatter no longer runs one statement over millions of rows.

### Changed
- **PERF.180: the three object families of the cache warm-up are slices of their listing.** `warm-sectors`, `warm-systems` and `warm-phenomena` were one unit each, minutes on a large galaxy once the cap is raised; each is now units of 250 objects (a whole multiple of that past 4,000 units), declared at 75 s against the 300 s unit ceiling. The stage leaves one line in the activity log for all its units.
- **Emptying the bright-star or phenomenon scatter removes its member and node rows in a way that cannot run for minutes in one statement.** The two tables are truncated when the other scatter has none, else the rows go in batches of 100,000 that each commit. No stored value changes; the timing on a large galaxy still has to be measured.
