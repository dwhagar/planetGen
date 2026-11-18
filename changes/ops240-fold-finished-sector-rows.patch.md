**Headline:** a galaxy run of many thousands of sectors no longer keeps a kilobyte of history for every sector it made.

### Changed
- OPS.240: when a stage that has more than 10,000 unit rows ends, the unit rows that finished at their first attempt lose their stored arguments and their attempt row (`planrows.fold_finished`), except the 100 slowest. Every row stays, so the stage's count, its status lines and the Queue page's totals are as real as before; a unit that was retried or failed keeps everything. What was folded (rows, attempts, their seconds, the rows kept) is written to the stage's `release` row. A failed fold never fails the stage's end.

generation-output: unchanged (no golden, no epoch, no lock item)
