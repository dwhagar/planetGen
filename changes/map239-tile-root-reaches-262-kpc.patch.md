**Headline:** The Galaxy Map's tile grid now covers the whole 262 kpc reach, so a tile exists for every address the grid can hold, in the halo and beyond the disc.

### Changed
- `TILE_ROOT_EDGE_PC` in `galaxy/viewport.py` is twice the reach (`cellkey.REACH_PC`), 524,288 pc, where it was 65,536 pc (nothing was tiled beyond 32.8 kpc). The finest tile is still 16 pc, now level 15 (`TILE_MAX_LEVEL`, derived from the two). The browser takes the root and the finest level from the page data; its fallbacks were changed to match.
- A tile of a given size is three levels deeper in its key than it was. The per-level budgets in `db/query.py` (generated stars, bright stars, point phenomena) and `POINT_PHENOMENON_MIN_LEVEL` are shifted by `TILE_ROOT_EXTRA_LEVELS`, so each size keeps its budget; the three new coarsest levels take the coarsest level's bright-star budget and list no generated stars. The tile cache starts a new generation with the release, as it does for any release, so no old key is read as a new one.
- `docs/api.md`, the architecture document, the precision audit (item 5.6) and the zoom-star ladder say so.

### Added
- `test_tiles_cover_the_reach`: every corner of the reach, to its last representable value, is in one tile of every level, and nothing beyond it is.
- The tests that named a tile level now name it through `TILE_MAX_LEVEL` and `TILE_ROOT_EXTRA_LEVELS`.

TODO items resolved: MAP.239 / Added: none
