**Headline:** The cache warm-up now finishes on a large galaxy, in about 40 percent less time.

### Fixed
- **`planetgen warm` stopped after its fourth phase on a large galaxy (PERF.156).** A slow page answers `202` with a progress bar while a job builds it, and the warm-up counted that as a page that did not answer. The warm-up now builds a slow page in its own request, so the answer is cached when it ends.
- **A tile of the Galaxy Map lost its bright stars when its query ran past three seconds (PERF.157).** The bright-star piece chose how to read by the box's share of the disc's area, not its height, and walked a whole population index for boxes far above the plane; it chooses by the box's share of the stars now. Of 660 sampled tiles at the default galaxy size 58 were served without their bright stars before, 1 after.

### Changed
- **The tile builder reads far fewer rows (PERF.157).** The scattered-point query lists a range per shell instead of one range over all of them, which for a 256 pc tile at the centre was a scan of the 3.6 million row table (22.4 million rows over 83 sampled tiles, 1.6 million now). The warm-up at the default size takes 6.5 minutes instead of 10.8.
- **Scattered multi-star systems set their member counts in one statement per batch (PERF.158).** The band scatter at the default size takes 536 seconds instead of 617.
