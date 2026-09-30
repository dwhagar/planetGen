### Changed
- **The Galaxy Map moved to the Flask app at `/galaxy`.** It is a plain,
  bookmarkable GET URL with no database name in it
  (`/galaxy?quadrant=II&page=2` for one Quadrant's sector list), under
  the new header with the Galaxy section marked and breadcrumbs. The
  map's "View sector" button and every table row are real `<a href>`
  links. The map's script fetches its tiles from `/galaxy/tiles` (JSON,
  no database in the URL). Tiles are still cached on the server's disk,
  now by the WSGI process, and in the browser's `localStorage` under the
  same keys as before, so nothing already cached is lost. The map's
  viewport is larger on wide screens.
- `galaxy.py` and `galaxy_tiles.py` now answer 301 to `/galaxy` and
  `/galaxy/tiles`, keeping their parameters, so old links, bookmarks and
  open tabs still work. A tile cache directory set only with `SetEnv
  PLANETGEN_TILE_CACHE_DIR` in the vhost no longer applies (the WSGI
  daemon doesn't see `SetEnv`); set `tile_cache.dir` in `config.json`
  instead.
