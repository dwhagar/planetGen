### Added

- The Sectors table (home page and `/sectors`), the All Systems table
  (`/systems`) and the Standalone Systems table are data tables like the
  Phenomena list (UX.41): click a header to sort, filter from the menus
  above (sectors by Quadrant; systems by where they are, single or binary,
  and octant) and scroll through every row, 50 fetched at a time. Each table
  on a page keeps its own sort, filters and place in the address bar
  (`sectors_sort`, `systems_sort`, `standalone_sort`, ...). Sectors are
  still nearest the core first until you sort them.
- `GET /api/sectors` and `GET /api/systems` take `sort`, `order`, filters and
  `facets=1`.

### Changed

- The Sectors table shows its density and distance in plain text with
  Unicode superscripts instead of HTML ones.
