**Headline:** The progress lines of the long loops (block light totals, sector stages, downloads) work again: the ticker called a time formatter that the message text standard had removed, so its lines would have failed.

### Changed
- **`util/ticker.py` writes its time with `logrecord.duration_text` (OPS.173).** `joblog.clock_text` is gone; a ticker line now reads `Things: 1,240 of 4,000 done (31%), 41.3 per s, 30.0 s elapsed.` and the stamp `Adding up the levels (4:07 elapsed).` keeps its form from a minute on.

generation-output: unchanged
generation-affects: none
