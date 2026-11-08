**Headline:** A run killed while it was cutting the next stage into units now says so, names that stage, and logs what the cut is doing and how much memory it holds; the starting-sector and scatter workers are sized from the full-galaxy run that went over.

### Changed
- **The opening of a stage is logged and marked (PERF.252).** The unit that completes a stage opens the next one in its own work horse, and the next stage's start line is written only after its cut into units returns, so a slow or killed cut was silent and a kill read as a unit of the stage before it. The cut now logs `Opening stage N: ...` before it starts and after the cluster field, the preparation and the cut, each with the memory the process holds, and the run record names the stage being opened until the cut ends. A kill or a failure then reads `unit 0 of stage X was killed while it was opening stage Y`.
- **Worker sizes (PERF.252).** `SECTOR_WORKER_MEMORY_MB` 2,600 to 3,050 (the starting sector held 2,762 MB) and `BULK_WORKER_MEMORY_MB` 1,700 to 2,100 (the bright stars held 1,879 MB), each 10 percent over the measure. On a 16 GB host a scatter stage now runs 3 workers, a sector stage 2.

Output unchanged.
