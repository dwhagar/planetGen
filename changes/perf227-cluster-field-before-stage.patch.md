**Headline:** A scatter or generation stage now makes the galaxy's cluster field once before its workers start, instead of every worker building it for itself.

### Changed
- Before it is cut into units, a stage whose workers read the galaxy's shape first thing (the massive-star, bright-star and phenomena scatters and the other stages that warm the shape) makes the cluster field a file in the state folder (`galaxy.field_job.ensure_file`, called from `queue.wrapper._expand`). The warm-up that normally makes that file refuses to run while a generation run is live, and a new or re-planned galaxy has none, so every worker used to build the field itself: a minute of work and several hundred megabytes in each, ten at a time.

### Fixed
- Starting a scatter on a galaxy with many clusters no longer leaves the run at "Massive stars, 2 of 456 units" for tens of minutes with ten workers holding about 7 GB between them, which also kept every job started after it waiting at position 1.

generation-output: unchanged
TODO items resolved: PERF.227
