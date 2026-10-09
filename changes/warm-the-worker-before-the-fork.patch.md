### Changed
- The queue worker loads the generation code and the bright-star sampling tables once, before it forks a work horse per job, so each job no longer spends about 2.6 s importing and rebuilding them (PERF.42). A worker also restarts itself when an update changes the release.
