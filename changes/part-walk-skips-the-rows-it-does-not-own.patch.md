**Headline:** The bright-star backfill no longer works out sky rows another part owns, so each part runs 7 to 35 percent faster with identical stars.

### Changed
- PERF.280: each part of the bright-star backfill (OPS.253) now skips the window arithmetic of the sphere rows the other parts own, instead of working every row out and dropping 15 of 16. The cells and stars are identical (the cell lists of two parts of three sky tiers hash the same before and after); a part of the deep tiers walks 7 to 35 percent faster.
