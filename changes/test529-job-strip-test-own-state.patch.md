**Headline:** The admin job strip test no longer depends on whether the machine it runs on has Redis, a writable jobs folder or jobs of its own.

### Changed
- **Test fix, no code change.** The strip is hidden when Redis or the jobs folder cannot be read, so on a machine without a Redis server the test saw it hidden while a job ran and failed. The test now sets up its own job state (an empty jobs folder, nothing waiting, no running lock) and still fails when the strip stops showing a running job.
