**Headline:** warming the map no longer goes quiet for a minute while it counts how much disk the tile cache uses, and it counts less often.

### Fixed

- PERF.274: `cache_is_full` walked every folder of the tile cache before each chunk of nodes. On a cold disk with 181,000 files that walk alone took 58 s of silence. The walk now tells the unit beat at each folder it enters (`tilecache.ON_DIRECTORY`, set by the warm-up to `unitbeat.moving`), and a walk that took at least a second is trusted for ten times as long as it took before the next one. Small caches are still read fresh every time, so what the cache says when it is nearly full is unchanged.
