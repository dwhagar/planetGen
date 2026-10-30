**Headline:** The block light totals that a scatter or `fill-block-totals` rebuilds at the end are written about a fifth faster in the part that fills the table, and the rows are the same.

### Changed

- **The block light rebuild writes its rows as plain INSERT text (PERF.190).** On a galaxy of 3 million pre-placed stars and 1.5 million scattered phenomena (19.6 million rows) the shell scan is about 1.5 percent of the rebuild; nearly all of it is the 83 units, whose rows were turned into SQL one value at a time by the driver. They are now formatted as the batched inserts are, so a unit takes 10.7 s against 13.3 s with the same stored rows. The rest of a unit is the per-shell grouped reads, the sums rolled up in Python and the server's own insert time.

generation-output: unchanged

TODO items resolved: PERF.190.
