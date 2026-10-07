### Changed
- **The work queue and job runner move into `planetgen.queue` and `planetgen.cli.job` (OPS.24, step 9 of 14).** `workQueue`, `progressFile`, `progressRate` and `systemLoad` are now `planetgen.queue.work`, `progress_file`, `progress_rate` and `load`. `src/jobRunner.py` is gone; the Generate page starts `python -m planetgen.cli.job` from the checkout's `src/`. Every caller moved too.
