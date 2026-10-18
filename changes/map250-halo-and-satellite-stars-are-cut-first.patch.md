**Headline:** When the Galaxy Map has room for only some of the stars, the stars of the halo and of the satellite galaxies are always the ones left out first.

### Changed
- Every bright-star list ranks by a tier before the MAP.154 brightness key: the disc, the bulge, the nuclei and the globular clusters first, then the stellar halo, then the satellite galaxies (`LOWEST_PRIORITY_TIER`, `_top_by_key` in `planetgen/db/query.py`). Before, a halo or satellite star was ranked by brightness like any other, so a tile over a loose group could list its stars ahead of the disc's. The order is the same at every tile level, so a child tile still lists every star of its parent inside it.
- It covers the tile lists (all three ways `galaxy_bright_stars_in_box` reads), the galaxy-wide sample that answers the biggest tiles, the finest tiles' detail lists and the stars round a nebula.
- `docs/design/zoom-star-visibility.md` says what is built, and what is not covered and why (generated stars carry no population and their sectors lie inside the outline; a sector page cuts no star).

### Added
- `test_halo_and_satellite_stars_are_listed_last_whatever_their_brightness` (`tests/test_db_persistence.py`): halo and satellite stars brighter than every disc and bulge star are listed after them, the satellite galaxies' last, and no disc or bulge star is left out while any of them is listed, in a tile that reads its addresses, one that walks the luminosity index and the sample.

TODO items resolved: MAP.250 / Added: none
