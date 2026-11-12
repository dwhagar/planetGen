**Headline:** The unused whole-job shell rate reader is retired; the per-unit ETA already does its job.

### Changed
- **PERF.147 closed on the reader side (GEN.1559).** `generation.stats.job_shells_per_second`, built to start the ETA of a task over 15 seconds, is deleted: the ETA now comes from `queue/taskprogress.py` and the per-command rates of `queue/unitrates.py`, and the stage times of `queue/runeta.py`, so nothing read the whole-job rate. The whole-job history row is still stored. Its register line is removed.
generation-output: unchanged
