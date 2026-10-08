### Added

- A sector's Contents and a Galaxy Map Quadrant's sectors are data tables
  like the others (UX.41): sort from the headers, filter Contents by Type and
  Octant and a Quadrant by Zone, and scroll through every row. Many rogue
  planets still read as one folded row; it links to the Contents filtered to
  Rogue Planet, which lists them one by one.
- The Contents "Show on map" buttons keep working as rows scroll in and out.
- `GET /table/sector-contents?sector=<id>` and
  `GET /table/galaxy-quadrant?quadrant=<I-IV>` serve their rows to the page.

### Changed

- The Contents table's folded rogue-planet group is no longer an expandable
  `<details>` row; see above.
