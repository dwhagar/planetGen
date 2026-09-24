### Changed
- **The browse page's sector list now runs outward from the galactic
  core.** Sectors are sorted by distance from the core (nearest first),
  with sectors never placed in a galaxy listed after them by name, and a
  new "Distance from core" column shows each one's distance. Paging walks
  the same order.
- **A sector's systems table now runs outward from the sector's center.**
  Systems are sorted by distance from the sector's center (nearest first),
  with a new "From center" column. `GET /api/sectors` gains
  `galactic_radius_pc`/`galactic_radius_ly` and `GET /api/sectors/<id>`'s
  systems gain `center_distance_ly`.
