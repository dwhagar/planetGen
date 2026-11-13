**Headline:** Two leftover back-compatibility paths for results and rows saved before 8.0 are gone.

### Changed
- **TEST.736.** `radiation.total_dose_msv_yr` now requires the belt part (every stored body has one, 0.0 outside a belt) instead of reading a missing one as zero, and the benchmark report no longer tolerates a result saved before the buffer pool, table sizes and worker start-up were recorded. Their two tests that kept the old behaviour alive are dropped and the belt-part assertions stay. Nothing generated changes.
generation-output: unchanged
