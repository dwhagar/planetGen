**Headline:** The memory-reading test of the worker cap no longer fails under gate load.

### Changed

- `test_peak_counts_the_worker_beside_its_horse` (TEST.757) held `peak_resident_mb` to a reading of the worker and its horse taken afterwards; the system divides a shared page by the processes sharing it, so the later reading can be higher when another test worker ends in between. The test now takes the reading before and after the call and holds the figure to the lower of the two.
