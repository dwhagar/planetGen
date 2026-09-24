### Changed
- **`install.sh` and `update.sh` now create the Galaxy Map's tile cache
  directory** (`/var/cache/planetgen/tiles`, or `tile_cache.dir` from
  `config.json`) and give it to Apache's user, so the cache doesn't fall
  back to a private `/tmp` folder that's cleared on every restart. The
  step is also available on its own as
  `sudo examples/apache/create-cache-dir.sh [dir]`.
