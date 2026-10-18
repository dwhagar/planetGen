**Headline:** `planetgen warm` now builds the Galaxy Map's tiles ahead of the first visitor, the first levels in full and the zoom run in to the centre.

### Added
- **`warm-galaxy-tiles`, a new cache family and command (PERF.130).** `planetgen warm` (and so every update) now builds every tile of the first `tile_cache.warm_levels` levels of the map's tile octree (default 8: 4,713 tiles on the scaled test galaxy), then the tiles of a zoom run from the opening view in to the galaxy's centre along the opening camera's line, down to level `tile_cache.warm_zoom_levels` (default 14, the finest the map's cut uses; 7,112 tiles in all). The work is cut into chunks of 128 tiles that the queue's workers build, coarsest level first, so a stopped run has built what the most views use; a second run finds every tile stored and builds none. On the scaled test galaxy the defaults took 12 minutes 3 seconds in one process and 40 MB, and the time for a fresh visitor to get every tile of a view along the zoom run fell from 111 seconds over eight views (53.5 s at 3,000 pc) to 3.2 seconds. A level count that would build more than 150,000 tiles is refused with its count. Two new settings, `tile_cache.warm_levels` and `tile_cache.warm_zoom_levels`, in `config.json`, the example, the schema and `docs/config.md`.

### Changed
- `planetgen warm` lists seven families instead of six (`galaxy-tiles` follows `galaxy-map`).

TODO items resolved: PERF.130 / Added: none
