**Headline:** The Galaxy Map now draws the sector-block nodes by their size on screen and the cube tiles are gone, so a zoom no longer swaps one set of tiles for another.

### Changed
- **The map draws nodes, not cube tiles (MAP.261).** `static/galaxynodecut.js` (twin of `planetgen.galaxy.nodecut`) picks the blocks to fetch and draw by projected size: a block opens into its children as it grows on screen, between 70 and 140 px, inside the cone the camera looks along and a budget of nodes, nearest first. Each node has a weight that its ancestors' hand-over leaves it, a block whose children are not yet loaded keeps all of its weight, and `starfade.js` shows a star at the weight of the node that lists it, so a parent and its children always add up to one and a star never dips while the children arrive. The page fetches the children of what it shows before they are needed (still camera, at most 12 nodes, header `X-Map-Prefetch`) and holds nodes in IndexedDB `planetgen-nodes` and in local storage under `planetgen:node:`.
- **Warm-up settings.** `tile_cache.warm_zoom_floor` (default 9, the finest block level of the zoom run) replaces `tile_cache.warm_levels` and `tile_cache.warm_zoom_levels`; the admin settings page has one stepper for it. `planetgen warm` has the families `galaxy-nodes` and `galaxy-deep-nodes` and the phases `built-levels` and `deep-nodes`; the `warm-galaxy-tiles` and `warm-galaxy-deep-tiles` commands are gone.
- **`TILE_FORMAT` is 3.** Every stored map file is built again once, as the cache folder fills. The folder, the `tile_cache` settings section and the packed form keep their names.

### Removed
- **Cube tiles (MAP.262).** `GET /galaxy/tiles`, `GET /api/galaxy/tiles`, `viewport.distance_cut_tiles` and its twin `galaxytilecut.js`, the tile read in `db/query.py`, the tile-keyed cache files and the tile entries of the change feed (`/api/galaxy/changes` no longer has `tiles`; it has `nodes`). The API version is 8. The rate limit `galaxy_tiles` ("Map data limit") now covers `/galaxy/nodes`.

### Fixed
- **A block wholly above or below a thin shell no longer errors.** The bright-star band read took the arc cosine of a ratio that passes 1 for such a block; the opening request of `/galaxy` answered 400 and the page 502. The ratio is clamped.
