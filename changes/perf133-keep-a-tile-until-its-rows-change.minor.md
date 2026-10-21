**Headline:** A tile of the Galaxy Map is built again only when something it holds has changed, and a release alone no longer throws the stored tiles away.

### Changed
- **The stored tiles no longer follow the release number (PERF.133).** They are keyed by a new tile stamp (`db.query.galaxy_tile_stamp`, in `GET /api/galaxy/changes` as `tile_stamp` and `tile_full`) made of the galaxy shape, the bright-star seed, the naming key, a tile-format number (`db.query.TILE_FORMAT`) and the generator epoch, with the sector change counters. A release that changes neither the tile format nor the epoch keeps every stored tile, so the `planetgen warm` that `update.sh` starts after it finds them stored and builds none; the pages and counts cached by the old stamp are still dropped by a release. A change to what a tile holds or how it is packed raises `TILE_FORMAT`, and `tests/test_tile_format.py` fails until it does.
- **A change rebuilds what it touched at once.** When a sector is edited, added or given a system, only the tiles holding it (one per level) are deleted, and those that were stored are queued on the bulk queue for a rebuild straight away, coarsest first (at most 2,048 a change), so the next visit finds them built. It is off when `tile_cache.prefetch_per_minute` is 0, during an update and while the cache is full.

TODO items resolved / Added: PERF.133 / Added: none
