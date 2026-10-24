**Headline:** The size and time estimate for a generation run is no longer inflated by the scatter's tables, the progress bars keep their elapsed and remaining times in a job log, and the wipe line no longer quotes a stale row count (PERF.178, PERF.179, UX.164).

### Changed

- The measured bytes per star system now leave out the scatter's own tables (black holes and other phenomena, their companions, the galaxy blocks, components and containers). Boss's log estimated 33 GB for 105 sectors because the 24 million black holes were counted as if the systems had made them.
- The default time per system when no speed is measured yet rises from 0.2 s to 0.35 s, from the log (2m33s estimated, 4m28s taken for 105 sectors).
- The progress bars in a job log (output that is not a terminal) are drawn 160 columns wide instead of rich's 80, so "elapsed" and "remaining" are no longer cut off.
- `planetgen reset` prints "Wiped N table(s)" without the row count, which was the database's stale statistic and read "about 1 row" on a database that had held millions.
