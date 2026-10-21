**Headline:** One setting, `cache.max_mb` (1 GB by default), sizes the stored caches, and when it is full the warm-up stops and only tiles nobody asked for lately are dropped.

### Changed
- **One cache size key, `cache.max_mb` (PERF.139).** It replaces `tile_cache.max_mb` and `PLANETGEN_TILE_CACHE_MAX_MB` with no compatibility shim (an old `tile_cache.max_mb` in `config.json` is ignored), defaults to 1024 MB instead of 200, and `0` turns the disk caches off. The Galaxy Map's tile cache is the first thing under it; the object index and the stored lists will share it. `update.sh`'s cache directory step and `docs/config.md` follow the new key.
- **What happens at the budget (PERF.133).** The tiles `planetgen warm`, the background build and the rebuild after a change build ahead are marked, and a visit sets the modification time of the file it reads. When the cache is over the budget `prune` deletes only the files a visit asked for, least recently visited first, down to 80%; a tile built ahead is never deleted for size. `planetgen warm` then stops building tiles and says how many it left out (`N tiles left out, the tile cache is full`).

TODO items resolved / Added: PERF.139 (key and budget behaviour; the slider is ADM.97), PERF.133 (second part) / Added: none
