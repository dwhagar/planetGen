### Fixed
- **The Galaxy Map could take the whole site down (Apache OOM-killed on
  2026-09-24).** One view request listed every not-yet-generated sector
  slot within 200 pc of the camera target (~770,000 of them) before
  keeping the nearest 4,000: ~400 MB and up to 140 s per request, several
  at once across the API's threads. Slot search now skips slots outside
  the view's azimuth (same results, ~100x faster), keeps only the nearest
  results in memory, and stops at 40 pc.

### Added
- **The Galaxy Map loads space one cube ("tile") at a time, like a web
  map.** New `GET /api/galaxy/tiles` returns fixed cubes of an octree
  (placed sectors capped per cube, not-yet-generated slots only for
  16 pc cubes) and `GET /api/galaxy/stamp` a token that changes when the
  galaxy's contents do. No request can scan an unbounded region.
- **Tiles are cached on the server's disk and in the browser.** The web
  layer (`html/lib/tilecache.py`, behind the new `galaxy_tiles.py`) only
  asks the API for tiles it hasn't cached, and the map keeps tiles in
  `localStorage`, so revisiting or reloading doesn't refetch them. New
  sectors show up within a minute. Configure with `tile_cache.dir`/
  `tile_cache.max_mb` or `PLANETGEN_TILE_CACHE_DIR`/
  `PLANETGEN_TILE_CACHE_MAX_MB`.

### Changed
- The example Apache config runs the API daemon with
  `request-timeout=60`, so a runaway request restarts the daemon instead
  of piling up until the OOM killer stops Apache.
- The planned-slot ball no longer appears around the galactic core at
  full-galaxy zoom; planned slots show within 20 pc of the camera target
  once zoomed in.
